#include "SlopeSolver.h"

#include "TPZStructMatrix.h"
#include "pzstepsolver.h"

#include <algorithm>
#include <cmath>

namespace slope {

const REAL kTol = 1.e-7; // |R| / |fbar| nas equacoes livres (a tolerancia do TPZUnsatAnalysis)

REAL PathProblem::Dot(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) const
{
    REAL s = 0.;
    for (int64_t i : free) s += a.GetVal(i, 0) * b.GetVal(i, 0);
    return s;
}

const char *Name(FactorMethod m)
{
    switch (m) {
    case FactorMethod::GIMArcLength: return "GIM arc-length";
    case FactorMethod::GIMBisection: return "GIM bisseccao";
    case FactorMethod::SRMArcLength: return "SRM arc-length";
    case FactorMethod::SRMBisection: return "SRM bisseccao";
    }
    return "";
}

void SetFreeEquations(TPZLinearAnalysis &an, const std::vector<int64_t> &free)
{
    TPZStack<int64_t> act;
    for (int64_t i : free) act.Push(i);
    an.StructMatrix()->EquationFilter().Reset();
    an.StructMatrix()->EquationFilter().SetActiveEquations(act);
}

TPZFMatrix<STATE> SolveK(TPZLinearAnalysis &an, const TPZFMatrix<STATE> &b)
{
    TPZEquationFilter &filtro = an.StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> br(filtro.NActiveEquations(), 1, 0.), yr, y(b.Rows(), 1, 0.);
    filtro.Gather(b, br);
    an.MatrixSolver<STATE>().Solve(br, yr); // Solve do TPZLinearAnalysis carregaria y na malha
    filtro.Scatter(yr, y);
    return y;
}

bool SolveFree(TPZLinearAnalysis &an, const std::vector<int64_t> &free, TPZFMatrix<STATE> &x)
{
    SetFreeEquations(an, free);
    REAL r0 = 0.;
    for (int it = 0; it < 5; it++) {
        {
            QuietCout q;
            an.LoadSolution(x);
            an.Assemble();
        }
        const TPZFMatrix<STATE> rhs = an.Rhs();
        REAL r = 0.;
        for (int64_t i : free) r += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
        r = std::sqrt(r);
        if (it == 0) r0 = r;
        if (r <= 1.e-10 * r0) return true;
        x += SolveK(an, rhs);
    }
    return false;
}

// x com as fixas em lambda x0 (scaleFixed)
static TPZFMatrix<STATE> Loaded(const PathProblem &p, const TPZFMatrix<STATE> &x)
{
    TPZFMatrix<STATE> xl(x);
    if (p.scaleFixed)
        for (int64_t i : p.fixed) xl(i, 0) = p.lambda * p.x0.GetVal(i, 0);
    return xl;
}

// -R(x) nas equacoes livres e a tangente no solver
static TPZFMatrix<STATE> AssembleAt(PathProblem &p, const TPZFMatrix<STATE> &x)
{
    QuietCout q;
    p.an->LoadSolution(Loaded(p, x));
    p.an->Assemble();
    return p.an->Rhs();
}

// So -R(x) (a tangente fatorada fica)
static TPZFMatrix<STATE> ResidualAt(PathProblem &p, const TPZFMatrix<STATE> &x)
{
    QuietCout q;
    p.an->LoadSolution(Loaded(p, x));
    p.an->AssembleResidual();
    return p.an->Rhs();
}

// Newton com os parametros ja ajustados, a partir de x
static bool Newton(PathProblem &p, TPZFMatrix<STATE> &x, REAL fref, int maxit, int &it)
{
    for (int k = 1; k <= maxit; k++) {
        it = k;
        const TPZFMatrix<STATE> rhs = AssembleAt(p, x);
        const REAL e = std::sqrt(p.Dot(rhs, rhs)) / fref;
        if (!std::isfinite(e)) return false;
        if (k > 1 && e < kTol) return true;
        try {
            x += SolveK(*p.an, rhs);
        } catch (...) { // tangente singular
            return false;
        }
    }
    return false;
}

// Parametro do caminho: set(v, preditor) e a carga tangente fbar = -dR/dv em x (rhs = -R em x)
struct PathParam {
    std::function<void(REAL, bool)> set;
    std::function<TPZFMatrix<STATE>(const TPZFMatrix<STATE> &, const TPZFMatrix<STATE> &)> load;
};

struct ArcResult {
    bool conv = false, vertical = false;
    int iter = 0;
    REAL value = 0.;
    TPZFMatrix<STATE> du, x; // incremento e estado final
};

// Um incremento do arc-length cilindrico (Box 4.4), como TPZUnsatAnalysis::ArcLengthIterate: preditor secante (o
// incremento anterior com comprimento l) e raiz (4.118) de menor angulo
static ArcResult ArcLengthStep(PathProblem &p, const PathParam &pp, const TPZFMatrix<STATE> &X, REAL fref, REAL vn,
                               REAL l, const TPZFMatrix<STATE> &dun, REAL dvn, int maxit)
{
    ArcResult res;
    res.du.Redim(X.Rows(), 1);
    TPZFMatrix<STATE> x(X), &Du = res.du;
    REAL v = vn;
    const bool hasPrev = dun.Rows() == X.Rows();
    int kfirst = 1;
    if (hasPrev) {
        const REAL prev2 = p.Dot(dun, dun);
        if (prev2 > 0. && l * l < 1.e4 * prev2) {
            const REAL s = l / std::sqrt(prev2);
            for (int64_t i : p.free) x(i, 0) += s * dun.GetVal(i, 0), Du(i, 0) = s * dun.GetVal(i, 0);
            v += s * dvn;
            kfirst = 2;
        }
    }
    for (int k = kfirst; k <= maxit; k++) {
        res.iter = k - kfirst + 1;
        if (!(v > 0.)) break;
        pp.set(v, k == 1);
        const TPZFMatrix<STATE> rhs = AssembleAt(p, x);
        const REAL e = std::sqrt(p.Dot(rhs, rhs)) / fref;
        if (!std::isfinite(e) || e > 1.e6) break;
        if (k > 1 && e < kTol) {
            res.conv = true;
            break;
        }
        const TPZFMatrix<STATE> fbar = pp.load(x, rhs);
        TPZFMatrix<STATE> dustar, dubar;
        try {
            dustar = SolveK(*p.an, rhs);
            dubar = SolveK(*p.an, fbar);
        } catch (...) {
            break;
        }
        const REAL a = p.Dot(dubar, dubar);
        if (!(a > 0.)) {
            res.vertical = a == 0.;
            break;
        }
        REAL dv;
        if (k == 1) { // preditor tangente, com o sinal do criterio secante (4.123)
            dv = (hasPrev && p.Dot(dun, dubar) < 0. ? -1. : 1.) * l / std::sqrt(a);
        } else {      // restricao |Du| = l (4.116)-(4.118)
            TPZFMatrix<STATE> w(Du);
            w += dustar;
            const REAL b = 2. * p.Dot(w, dubar), c = p.Dot(w, w) - l * l, disc = b * b - 4. * a * c;
            if (disc < 0.) break;
            const REAL q = -0.5 * (b + (b >= 0. ? 1. : -1.) * std::sqrt(disc));
            const REAL r1 = q / a, r2 = q != 0. ? c / q : r1, dir = p.Dot(dubar, Du);
            dv = r1 * dir >= r2 * dir ? r1 : r2;
        }
        v += dv;
        for (int64_t i : p.free) {
            const REAL di = dustar.GetVal(i, 0) + dv * dubar.GetVal(i, 0);
            x(i, 0) += di;
            Du(i, 0) += di;
        }
    }
    res.value = v;
    res.x = x;
    return res;
}

static void Row(int passo, REAL v, const char *dv, int it, REAL l, Clock::time_point t0, bool ok)
{
    char arco[16] = "", num[16] = "fixo";
    if (l > 0.) std::snprintf(arco, sizeof arco, "%11.3e", l);
    if (passo > 0) std::snprintf(num, sizeof num, "%d", passo);
    Print("    %-6s %10.5f %10s %4d %11s %9.1f%s", num, v, dv, it, arco, Seconds(t0), ok ? "" : "  falhou");
}

static std::string Delta(REAL d)
{
    char s[16];
    std::snprintf(s, sizeof s, "%+10.5f", d);
    return s;
}

// Continuacao com bisseccao (Continuation do SlopeAnalysis de referencia): passo dv, x 1,5 se convergiu em ate 8
// iteracoes, metade se falhou, ate dv < 0,002 v ou v = vmax
static FactorResult Continuation(PathProblem &p, TPZFMatrix<STATE> &X, REAL fref, const std::function<void(REAL)> &set,
                                 REAL v, REAL dv, REAL vmax, Clock::time_point t0, int &steps)
{
    TPZFMatrix<STATE> ultimo;
    while (v < vmax && dv > 2.e-3 * std::max(v, REAL(1.e-2))) {
        const REAL vt = std::min(v + dv, vmax);
        set(vt);
        TPZFMatrix<STATE> x(X);
        int it = 0;
        const bool ok = Newton(p, x, fref, 25, it);
        Row(++steps, vt, ok ? Delta(vt - v).c_str() : "", it, 0., t0, ok);
        if (ok) {
            p.accept();
            ultimo = x;
            ultimo -= X;
            X = x;
            v = vt;
            if (it <= 8) dv *= 1.5;
        } else {
            dv *= 0.5;
        }
    }
    set(v);
    return {v, steps, v >= vmax ? "limite" : "bisseccao", ultimo};
}

// Arc-length com os parametros do FactorGIMArco do SlopeRandomField: l ajustado para ~10 iteracoes, ate d v < 0,001 (no
// l pedido), v >= 10, l < 1e-4 |u| ou 100 incrementos; o fator e o maior v do caminho. Caminho vertical (nenhum ponto
// plastico, so no SRM): um passo de Newton com v 10% maior
static FactorResult ArcLengthPath(PathProblem &p, const PathParam &pp, TPZFMatrix<STATE> &X, REAL fref, REAL v,
                                  Clock::time_point t0, int &steps)
{
    const REAL tol = 0.001, l0 = 0.25, lMaxR = 2., lMinR = 1.e-4, vLimite = 10.;
    const int alvo = 10, maxIt = 12, maxPassos = 100;
    const REAL uRef = std::sqrt(p.Dot(X, X)), lMax = lMaxR * uRef, lMin = lMinR * uRef;
    REAL l = std::min(lMax, std::max(lMin, l0 * uRef)), lPedido = l, vMax = v, dvn = 0.;
    TPZFMatrix<STATE> dun, ultimo; // ultimo: o incremento do ultimo estado aceito
    std::string fim = "maxPassos";
    for (int passo = 0, verticais = 0; passo < maxPassos; passo++) {
        const ArcResult r = ArcLengthStep(p, pp, X, fref, v, l, dun, dvn, maxIt);
        Row(++steps, r.value, r.conv ? Delta(r.value - v).c_str() : "", r.iter, l, t0, r.conv);
        if (!r.conv && r.vertical) {
            const REAL vt = 1.1 * v;
            pp.set(vt, false);
            TPZFMatrix<STATE> x(X);
            int it = 0;
            if (++verticais > 20 || !Newton(p, x, fref, 25, it)) {
                fim = "vertical";
                break;
            }
            Row(++steps, vt, Delta(vt - v).c_str(), it, 0., t0, true);
            p.accept();
            ultimo = x;
            ultimo -= X;
            X = x;
            v = vt, vMax = std::max(vMax, v), dun = TPZFMatrix<STATE>();
            continue;
        }
        if (!r.conv) {
            l *= 0.5;
            if (l < lMin) {
                fim = "lMin";
                break;
            }
            continue;
        }
        p.accept();
        X = r.x;
        dvn = r.value - v;
        v = r.value;
        dun = r.du;
        ultimo = r.du;
        vMax = std::max(vMax, v);
        if (dvn * lPedido / l < tol) {
            fim = "patamar";
            break;
        }
        if (v >= vLimite) {
            fim = "limite";
            break;
        }
        l = std::min(lMax, l * std::min(2., std::max(0.5, std::sqrt(REAL(alvo) / r.iter))));
        lPedido = l;
    }
    return {vMax, steps, fim, ultimo};
}

FactorResult ComputeFactor(PathProblem &p, int64_t neq, FactorMethod method)
{
    if (p.x0.Rows() != neq) p.x0.Redim(neq, 1); // seco: fixas em 0
    std::vector<bool> livre(neq, false);
    for (int64_t i : p.free) livre[i] = true;
    p.fixed.clear();
    for (int64_t i = 0; i < neq; i++)
        if (!livre[i]) p.fixed.push_back(i);
    TPZFMatrix<STATE> X(p.x0); // estado aceito: sem tensoes, com as pressoes fixas
    p.SetReduction(1.);
    p.SetLambda(1.);
    TPZFMatrix<STATE> fbar = AssembleAt(p, X); // carga por unidade de lambda: -R(lambda = 1) + R(lambda = 0)
    p.SetLambda(0.);
    fbar -= AssembleAt(p, X);
    const REAL fref = std::sqrt(p.Dot(fbar, fbar));
    const auto t0 = Clock::now();
    const bool gim = method == FactorMethod::GIMArcLength || method == FactorMethod::GIMBisection;
    auto header = [](const char *v) {
        Print("    %-6s %10s %10s %4s %11s %9s", "passo", v, "d", "it", "arco l", "tempo (s)");
    };
    auto setLambda = [&p](REAL v) { p.SetLambda(v); };
    int steps = 0;
    FactorResult res;

    if (method == FactorMethod::GIMBisection) {
        header("lambda");
        res = Continuation(p, X, fref, setLambda, 0., 0.5, 10., t0, steps);
    } else if (method == FactorMethod::GIMArcLength) {
        header("lambda");
        REAL lam = 0.5; // primeiro passo com lambda prescrito
        for (;;) {
            p.SetLambda(lam);
            TPZFMatrix<STATE> x(X);
            int it = 0;
            const bool ok = Newton(p, x, fref, 12, it);
            Row(0, lam, "", it, 0., t0, ok);
            if (ok) {
                p.accept();
                X = x;
                break;
            }
            lam *= 0.5;
            if (lam < 0.01) return {0., steps, "sem equilibrio inicial"};
        }
        PathParam pp; // lambda; no preditor a superficie encolhe um pouco (tangente plastica nos pontos sobre ela)
        pp.set = [&p](REAL v, bool pred) { p.SetLambda(v), p.SetReduction(pred ? 1. + 1.e-6 : 1.); };
        pp.load = [&fbar](const TPZFMatrix<STATE> &, const TPZFMatrix<STATE> &) { return fbar; };
        res = ArcLengthPath(p, pp, X, fref, lam, t0, steps);
    } else {
        header("lambda"); // SRM: primeiro o peso total com a resistencia real
        const FactorResult g = Continuation(p, X, fref, setLambda, 0., 0.5, 1., t0, steps);
        if (g.value < 1.) {
            p.SetLambda(1.);
            return {0., steps, "sem equilibrio com o peso total (FS < 1)"};
        }
        header("F");
        if (method == FactorMethod::SRMBisection) {
            res = Continuation(p, X, fref, [&p](REAL F) { p.SetReduction(F); }, 1., 0.25, 10., t0, steps);
        } else {
            PathParam pp; // F; fbar = -dR/dF por diferenca finita (so os pontos plasticos contribuem)
            pp.set = [&p](REAL F, bool pred) { p.SetReduction(pred ? F * (1. + 1.e-6) : F); };
            pp.load = [&p](const TPZFMatrix<STATE> &x, const TPZFMatrix<STATE> &rhs) {
                const REAL F = p.reduction, h = 1.e-6 * F;
                p.SetReduction(F + h);
                TPZFMatrix<STATE> f = ResidualAt(p, x);
                p.SetReduction(F);
                f -= rhs;
                f *= 1. / h;
                return f;
            };
            res = ArcLengthPath(p, pp, X, fref, 1., t0, steps);
        }
    }
    p.SetReduction(1.);
    p.SetLambda(1.);
    return res;
}

} // namespace slope
