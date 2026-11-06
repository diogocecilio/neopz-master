#pragma once
// V2_TaludesCho: os dois taludes classicos de Cho (2010), como em Liao et al. (2025), secos: FS por reducao de
// resistencia (SRM) e P_f por Monte Carlo com campos aleatorios lognormais gerados por Karhunen-Loeve.
// Referencia: Liao, Zhao, Wu, Miao, Pan e Beer (2025), Efficient reliability analysis of slopes in spatially variable
// soils with active learning-assisted bootstrap PCE, Comput. Geotech. 179:107022; Cho (2010), JGGE 136(7):975-984.
//
//  - cho1 (Liao, sec. 4.1, Fig. 2, Tabela 1): argila saturada nao drenada, c_u medio 23 kPa, phi = 0 (Tresca, o limite
//    phi -> 0 do Mohr-Coulomb associado PV2), gamma = 20 kN/m3, E = 100 MPa, nu = 0,3. FS = 1,369 (Liao, SIGMA/W +
//    SLOPE/W) e 1,366 (Jiang et al. 2014). Campo lognormal de c_u: CoV 0,3, h_x = 25 m, h_z = 2,5 m; P_f = 10,12%
//    (Monte Carlo simples, 10 000 realizacoes, Tabela 2).
//  - cho2 (Liao, sec. 4.2, Fig. 6, Tabela 3): talude c-phi de 10 m a 1V:1H, c = 10 kPa, phi = 30 graus, gamma = 20,
//    E = 100 MPa, nu = 0,3. FS = 1,205 (Bishop), 1,204 (Cho 2010) e 1,284 (elementos finitos). Campos lognormais de c
//    (CoV 0,3) e phi (CoV 0,2), h_x = 20 m, h_z = 2 m, correlacao cruzada r = -0,5 nos coeficientes de KL;
//    P_f = 0,12% (MEF, 100 000 realizacoes, Tabela 4).
//
// Analise seca: poropressoes nulas e congeladas (so o equilibrio e resolvido; com p_w = p_a = 0 a tensao efetiva e a
// total). Geometria lida das Figs. 2 e 6 (coordenadas no README.md): malha estruturada de triangulos (u quadratico,
// T6) com espacamento h (0,5 m, o da grade de Liao). Apoios: base fixa nas duas direcoes, laterais fixas na
// horizontal (Liao, secs. 4.1 e 4.2). O FS e o do SRM do projeto (c/F e tan(phi)/F, fluxo associado), com a nova
// tentativa do F que falhou do projeto SlopeRandomField (branch UnsaturatedSoil).
#include "RahardjoDriver.h" // UnsatParams, UnsatMaterial, CreateSlopeCompMesh, DeleteCompMesh, ids dos contornos
#include "KLGauss.h"        // KL separavel com o nucleo gaussiano (Nystrom) e o exponencial (forma fechada)

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>

namespace v2 {

using Relogio = std::chrono::steady_clock;
inline double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }

/// 0: silencioso; 1: uma linha por tentativa do SRM (o executavel)
inline int &Verboso() { static int v = 0; return v; }

// ------------------------------------------------------------------------------------------------------------- casos

/// Um dos dois exemplos (os valores sao os de Liao et al. 2025; ver a tabela de dados do README.md)
struct Caso {
    std::string nome;
    // geometria (Figs. 2 e 6): base z = 0 em 0 <= x <= L; crista z = zc para x <= xc; face reta ate o pe (xp, zp);
    // z = zp para x >= xp
    double L = 30., xc = 10., zc = 10., xp = 20., zp = 5.;
    // solo (medias)
    double c = 23., phiGraus = 0., gamma = 20., E = 1.e5, nu = 0.3;
    // campos aleatorios: CoV, distancias de autocorrelacao h_x, h_z da eq. (1) e correlacao cruzada r (c x phi)
    double covC = 0.3, covPhi = 0., hx = 25., hz = 2.5, rho = 0.;
    bool phiAleatorio = false;
    // valores publicados
    double FSref = 1.369, FSref2 = 1.366, PfRef = 0.1012;
    int NmcRef = 10000;

    /// Cota da superficie do terreno em x
    double Superficie(double x) const {
        if (x <= xc) return zc;
        if (x >= xp) return zp;
        return zc + (zp - zc) * (x - xc) / (xp - xc);
    }
    /// Area da secao (m2)
    double AreaSecao() const { return xc * zc + 0.5 * (zc + zp) * (xp - xc) + (L - xp) * zp; }
};

/// Exemplo 1 (Liao, sec. 4.1): talude de 5 m a 1V:2H em argila nao drenada (Fig. 2: (0,0), (30,0), (30,5), (20,5),
/// (10,10), (0,10))
inline Caso Cho1() {
    Caso c;
    c.nome = "cho1";
    c.zc = 10., c.zp = 5.;
    c.c = 23., c.phiGraus = 0.;
    c.covC = 0.3, c.covPhi = 0., c.hx = 25., c.hz = 2.5, c.rho = 0., c.phiAleatorio = false;
    c.FSref = 1.369, c.FSref2 = 1.366, c.PfRef = 0.1012, c.NmcRef = 10000;
    return c;
}

/// Exemplo 2 (Liao, sec. 4.2): talude c-phi de 10 m a 1V:1H (Fig. 6: (0,0), (30,0), (30,5), (20,5), (10,15), (0,15))
inline Caso Cho2() {
    Caso c;
    c.nome = "cho2";
    c.zc = 15., c.zp = 5.;
    c.c = 10., c.phiGraus = 30.;
    c.covC = 0.3, c.covPhi = 0.2, c.hx = 20., c.hz = 2., c.rho = -0.5, c.phiAleatorio = true;
    c.FSref = 1.284, c.FSref2 = 1.205, c.PfRef = 0.0012, c.NmcRef = 100000; // FSref2: Bishop
    return c;
}

inline Caso CasoPorNome(const std::string &n) {
    if (n == "cho1") return Cho1();
    if (n == "cho2") return Cho2();
    throw std::invalid_argument("CASO deve ser cho1 ou cho2 (recebido: " + n + ")");
}

// ------------------------------------------------------------------------------------------------------------- malha

/// n intervalos iguais em [a, b], n = round((b - a) / h) >= 1
inline std::vector<double> Uniforme(double a, double b, double h) {
    const int n = std::max(1, (int)std::lround((b - a) / h));
    std::vector<double> x(n + 1);
    for (int i = 0; i <= n; ++i) x[i] = a + (b - a) * i / n;
    return x;
}

/**
 * Malha estruturada de triangulos do talude: linhas verticais com espacamento ~h (x = 0, xc, xp e L sao linhas da
 * malha, entao as quinas da superficie sao nos); em cada linha, round(z_s / h) intervalos iguais de z = 0 ate a
 * superficie z_s(x). A faixa entre duas linhas e triangulada avancando pelo lado cujo proximo no e mais baixo; com nos
 * na mesma cota (fora da face) cada quadrado vira dois triangulos pela diagonal paralela a face. Com h = 0,5 m e a face
 * a 45 graus (cho2) os nos sao os da grade de Liao (1281 nos; 2400 triangulos no lugar de 1210 quadrilateros e
 * triangulos). Ids: solo matSolo; contornos bcBase (z = 0), bcDireita (x = L), bcSuperficie e bcEsquerda (x = 0).
 */
inline TPZGeoMesh *MalhaTalude(const Caso &cs, double h) {
    if (!(h > 0.)) throw std::invalid_argument("MalhaTalude: h deve ser > 0");
    std::vector<double> xs = Uniforme(0., cs.xc, h);
    for (double v : Uniforme(cs.xc, cs.xp, h)) if (v > xs.back() + 1.e-9) xs.push_back(v);
    for (double v : Uniforme(cs.xp, cs.L, h)) if (v > xs.back() + 1.e-9) xs.push_back(v);
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const size_t nl = xs.size();
    std::vector<int64_t> inicio(nl);
    std::vector<int> nint(nl);
    int64_t nn = 0;
    for (size_t i = 0; i < nl; ++i) {
        nint[i] = std::max(1, (int)std::lround(cs.Superficie(xs[i]) / h));
        inicio[i] = nn;
        nn += nint[i] + 1;
    }
    g->NodeVec().Resize(nn);
    for (size_t i = 0; i < nl; ++i)
        for (int k = 0; k <= nint[i]; ++k) {
            TPZManVector<REAL, 3> x = {xs[i], cs.Superficie(xs[i]) * k / nint[i], 0.};
            g->NodeVec()[inicio[i] + k] = TPZGeoNode(inicio[i] + k, x, *g);
        }
    auto X = [&](int64_t no, int d) { return g->NodeVec()[no].Coord(d); };
    auto tri = [&](int64_t a, int64_t b, int64_t c) { // anti-horario (det J > 0)
        const double A = (X(b, 0) - X(a, 0)) * (X(c, 1) - X(a, 1)) - (X(c, 0) - X(a, 0)) * (X(b, 1) - X(a, 1));
        if (!(std::fabs(A) > 1.e-12)) throw std::runtime_error("MalhaTalude: triangulo degenerado");
        TPZManVector<int64_t, 3> v = {a, A > 0. ? b : c, A > 0. ? c : b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g);
    };
    for (size_t i = 0; i + 1 < nl; ++i) {
        int a = 0, b = 0;
        while (a < nint[i] || b < nint[i + 1]) {
            const int64_t La = inicio[i] + a, Rb = inicio[i + 1] + b;
            bool esquerda;
            if (a == nint[i]) esquerda = false;
            else if (b == nint[i + 1]) esquerda = true;
            else esquerda = X(La + 1, 1) <= X(Rb + 1, 1) + 1.e-9;
            if (esquerda) { tri(La, Rb, La + 1); ++a; }
            else { tri(La, Rb, Rb + 1); ++b; }
        }
    }
    auto lin = [&](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, id, *g);
    };
    const size_t N = nl - 1;
    for (size_t i = 0; i < N; ++i) lin(inicio[i], inicio[i + 1], bcBase);
    for (int k = 0; k < nint[N]; ++k) lin(inicio[N] + k, inicio[N] + k + 1, bcDireita);
    for (size_t i = N; i > 0; --i) lin(inicio[i] + nint[i], inicio[i - 1] + nint[i - 1], bcSuperficie);
    for (int k = nint[0]; k > 0; --k) lin(inicio[0] + k, inicio[0] + k - 1, bcEsquerda);
    g->BuildConnectivity();
    return g;
}

// ------------------------------------------------------------------------------------------------------------ modelo

/// Malha, material e analise de um caso; pontos de integracao do solo (indice de memoria, coordenadas e area)
class Modelo {
public:
    Caso cs;
    double h = 0.5;
    UnsatParams p;
    TPZGeoMesh *gmesh = nullptr;
    TPZMultiphysicsCompMesh *cmesh = nullptr;
    TPZUnsatAnalysis *an = nullptr;
    TUnsat *mat = nullptr;
    std::vector<int64_t> ips;
    std::vector<double> xs, ys, area;

    Modelo(const Caso &c, double h_, int maxitFrozen = 25, int porder = 2) : cs(c), h(h_) {
        p.E = cs.E, p.nu = cs.nu, p.c = cs.c, p.phi = cs.phiGraus * M_PI / 180., p.gamma = cs.gamma;
        p.mixture = false; // peso especifico constante (gamma_sat da Tabela 1 / 3)
        gmesh = MalhaTalude(cs, h);
        cmesh = CreateSlopeCompMesh(gmesh, porder, p);
        an = new TPZUnsatAnalysis(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
        an->SetTolerances(1.e-7, 1.e-6, 40, maxitFrozen);
        mat = dynamic_cast<TUnsat *>(cmesh->FindMaterial(matSolo));
        if (!mat) throw std::runtime_error("Modelo: material do solo nao encontrado");
        for (TPZCompEl *cel : cmesh->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
            TPZManVector<int64_t, 10> idx;
            cel->GetMemoryIndices(idx);
            const TPZIntPoints &rule = cel->GetIntegrationRule();
            if (idx.size() != rule.NPoints()) throw std::runtime_error("Modelo: memoria e regra de integracao diferem");
            TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
            TPZFNMatrix<9, REAL> jac, axes, jacinv;
            REAL w, detjac;
            for (int i = 0; i < rule.NPoints(); ++i) {
                rule.Point(i, qsi, w);
                cel->Reference()->X(qsi, x);
                cel->Reference()->Jacobian(qsi, jac, axes, detjac, jacinv);
                ips.push_back(idx[i]);
                xs.push_back(x[0]);
                ys.push_back(x[1]);
                area.push_back(w * std::fabs(detjac));
            }
        }
    }
    ~Modelo() {
        delete an;
        if (cmesh) DeleteCompMesh(cmesh);
        delete gmesh;
    }
    Modelo(const Modelo &) = delete;
    Modelo &operator=(const Modelo &) = delete;

    int64_t NElementos() const {
        int64_t n = 0;
        for (TPZGeoEl *gel : gmesh->ElementVec()) if (gel && gel->Dimension() == 2) ++n;
        return n;
    }
    double AreaPontos() const { double s = 0.; for (double a : area) s += a; return s; }

    /// c (kPa) e phi (rad) em cada ponto de integracao (psi = phi: associado)
    void Propriedades(const std::vector<double> &c, const std::vector<double> &phi) {
        auto mem = mat->GetMemory();
        for (size_t i = 0; i < ips.size(); ++i) (*mem)[ips[i]].m_elastoplastic_state.fmatprop = {c[i], phi[i], phi[i]};
    }
    /// As medias do caso em todos os pontos (o caso deterministico)
    void PropriedadesMedias() {
        Propriedades(std::vector<double>(ips.size(), p.c), std::vector<double>(ips.size(), p.phi));
    }
};

// --------------------------------------------------------------------------------------------------------------- SRM

/// Uma solucao da busca do SRM: um dos 4 passos de carga do peso com c/F0, tan(phi)/F0 ("gravidade") ou um F
/// tentado depois ("srm")
struct TentativaSRM {
    std::string fase;
    double g = 1., F = 1.;
    bool conv = false, retry = false;
    int iter = 0;
    double seg = 0.;
};

struct ResultadoSRM {
    double FS = -1.;      ///< ultimo F convergido (-1: sem equilibrio nem com F0)
    double Fparada = -1.; ///< onde o FactorSRM de Comum pararia (sem a nova tentativa)
    double F0 = 0.;       ///< fator do estado geostatico usado
    std::vector<TentativaSRM> tent;
    int newton = 0;
    double seg = 0.;
};

/// Lista de F do SRM: 0,3 a 3,0 de 0,1 em 0,1, depois 3,5, 4, 5 e 6 (o FS e localizado entre dois valores da lista e
/// refinado por bissecao ate a tolerancia)
inline std::vector<double> ListaF(double Fmax = 0.) {
    std::vector<double> f;
    for (int k = 3; k <= 30; ++k) f.push_back(k / 10.);
    for (double v : {3.5, 4., 5., 6.}) f.push_back(v);
    if (Fmax > 0.) { // busca censurada em Fmax: a lista termina em Fmax (FS >= Fmax se Fmax convergir)
        while (!f.empty() && f.back() >= Fmax - 1.e-12) f.pop_back();
        f.push_back(Fmax);
    }
    return f;
}

/**
 * FS por reducao de resistencia, seco (pressoes nulas e congeladas): o FactorSRM de Comum/RahardjoDriver.h com a nova
 * tentativa do F que falhou de FactorSRMRetry (SlopeRandomField, branch UnsaturatedSoil). Parte do estado sem tensoes
 * (a memoria volta ao padrao, mantendo c e phi de cada ponto); estado geostatico com c/F0, tan(phi)/F0 (peso em 4
 * incrementos); F cresce pela lista Fs a partir do ultimo estado convergido; no primeiro F sem convergencia, bissecao
 * ate o intervalo ficar menor que tol; entao o F que falhou e tentado de novo a partir do ultimo convergido (incremento
 * < tol): se convergir, a busca continua. Retorna o ultimo F convergido.
 */
inline ResultadoSRM FatorSRM(TPZUnsatAnalysis &an, double F0, double tol, const std::vector<double> &Fs) {
    const auto inicio = Relogio::now();
    ResultadoSRM res;
    res.F0 = F0;
    an.ResetMemory(); // mantem as propriedades de cada ponto (os campos)
    TPZFMatrix<STATE> &X = an.X();
    X.Redim(an.NEq(), 1); // deslocamentos e pressoes nulos (seco)
    an.Restore();
    const auto fixed = an.FrozenFixed(X);
    an.SetStrengthReductionFactor(F0);
    an.SetGravityScalesSuction(true);
    auto resolve = [&](const char *fase, double g, double F, bool retry) {
        const auto t0 = Relogio::now();
        const auto r = an.Solve(0., fixed, {}, true);
        res.tent.push_back({fase, g, F, r.conv, retry, r.iter, Segundos(t0)});
        res.newton += r.iter;
        if (Verboso())
            std::cout << "SRM: " << fase << " g = " << g << " F = " << F << (r.conv ? " convergiu" : " nao convergiu")
                      << " (" << r.iter << " iteracoes)" << (retry ? " [nova tentativa]" : "") << std::endl;
        return r.conv;
    };
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        if (!resolve("gravidade", REAL(k) / n, F0, false)) {
            an.Restore();
            an.SetGravityFactor(1.);
            an.SetStrengthReductionFactor(1.);
            res.seg = Segundos(inicio);
            return res; // FS < F0
        }
        an.Accept();
    }
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) {
                if (res.Fparada < 0.) res.Fparada = Fok;
                F = Ffail;
                retry = true;
            } else F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        if (resolve("srm", 1., F, retry)) {
            an.Accept();
            Fok = F;
            if (retry) Ffail = -1.;
        } else {
            an.Restore();
            if (retry) break;
            Ffail = F;
        }
    }
    if (res.Fparada < 0.) res.Fparada = Fok;
    an.SetStrengthReductionFactor(1.);
    res.FS = Fok;
    res.seg = Segundos(inicio);
    return res;
}

/// FatorSRM com F0 = F0, F0/2 e F0/5 ate o estado geostatico convergir (FS muito baixo em alguma realizacao).
/// Fmax > 0: busca censurada (para P_f basta saber se FS < 1: com Fmax = 1 as tentativas ate F = 1 sao as mesmas da
/// busca completa, entao a indicadora FS < 1 e identica; o FS informado e min(FS, Fmax))
inline ResultadoSRM FatorSRMSeco(TPZUnsatAnalysis &an, double F0, double tol, double Fmax = 0.) {
    ResultadoSRM r;
    int newton = 0;
    double seg = 0.;
    for (double f : {F0, F0 / 2., F0 / 5.}) {
        r = FatorSRM(an, f, tol, ListaF(Fmax));
        newton += r.newton;
        seg += r.seg;
        if (r.FS > 0.) break;
    }
    r.newton = newton;
    r.seg = seg;
    return r;
}

// ------------------------------------------------------------------------------------------------------- configuracao

inline double EnvD(const char *k, double d) { const char *v = std::getenv(k); return v && *v ? std::atof(v) : d; }
inline int EnvI(const char *k, int d) { const char *v = std::getenv(k); return v && *v ? std::atoi(v) : d; }
inline std::string EnvS(const char *k, const std::string &d) { const char *v = std::getenv(k); return v && *v ? std::string(v) : d; }

/// Parametros de uma execucao (os padroes sao os do artigo; o executavel le as variaveis de ambiente do README)
struct Config {
    std::string caso = "cho1";
    std::string nucleo = "gauss"; ///< "gauss" (rho2, o do artigo) ou "exp" (rho1, o de Comum/KLSeparavel.h)
    int real = 0, campos = 1;
    uint64_t semente = 20241220ULL;
    double covC = 0.3, covPhi = 0., hx = 25., hz = 2.5, rho = 0.;
    double klFrac = 0.95; ///< fracao de variancia da KL (o criterio epsilon >= 95% de Liao, eq. 4)
    int klM = 2000;       ///< maximo de termos da KL
    double h = 0.5;       ///< espacamento da malha (m)
    int porder = 2;       ///< ordem de u (2: T6)
    double fsTol = 0.01, F0 = 0.5;
    double fsMax = 0.;    ///< > 0: busca do FS censurada nesse valor (FS_MAX=1 basta para P_f; padrao: sem limite)
    int maxitFrozen = 25; ///< iteracoes de Newton antes de declarar a falta de convergencia num F
    double klMargem = 0.; ///< o retangulo da KL e [-m, L + m] x [-m, zc + m] (m = 0,25 m: o de Liao, ver README)
    double c = 23., phiGraus = 0., E = 1.e5, nu = 0.3, gamma = 20.;
    int camposCsv = 0;
    int soCampos = 0; ///< 1: so grava campos.csv (sem SRM nem resultado.txt), para ver os campos
    std::string versao = "-";

    /// Os padroes do caso (Tabelas 1 e 3 de Liao)
    static Config Padrao(const std::string &nome) {
        const Caso cs = CasoPorNome(nome);
        Config cf;
        cf.caso = nome;
        cf.covC = cs.covC, cf.covPhi = cs.covPhi, cf.hx = cs.hx, cf.hz = cs.hz, cf.rho = cs.rho;
        cf.c = cs.c, cf.phiGraus = cs.phiGraus, cf.E = cs.E, cf.nu = cs.nu, cf.gamma = cs.gamma;
        return cf;
    }
    /// Os padroes do caso CASO com as trocas das variaveis de ambiente
    static Config DoAmbiente() {
        Config cf = Padrao(EnvS("CASO", "cho1"));
        cf.nucleo = EnvS("NUCLEO", cf.nucleo);
        cf.real = EnvI("REAL", cf.real);
        cf.campos = EnvI("CAMPOS", cf.campos);
        if (const char *s = std::getenv("SEMENTE"); s && *s) cf.semente = std::strtoull(s, nullptr, 10);
        cf.covC = EnvD("COV_C", cf.covC);
        cf.covPhi = EnvD("COV_PHI", cf.covPhi);
        cf.hx = EnvD("LX", cf.hx);
        cf.hz = EnvD("LY", cf.hz);
        cf.rho = EnvD("RHO", cf.rho);
        cf.klFrac = EnvD("KL_FRAC", cf.klFrac);
        cf.klM = EnvI("KL_M", cf.klM);
        cf.h = EnvD("H_MALHA", cf.h);
        cf.porder = EnvI("ORDEM", cf.porder);
        cf.fsTol = EnvD("FS_TOL", cf.fsTol);
        cf.fsMax = EnvD("FS_MAX", cf.fsMax);
        cf.klMargem = EnvD("KL_MARGEM", cf.klMargem);
        cf.F0 = EnvD("FS_F0", cf.F0);
        cf.maxitFrozen = EnvI("MAXIT_FROZEN", cf.maxitFrozen);
        cf.c = EnvD("COES", cf.c);
        cf.phiGraus = EnvD("PHI_DEG", cf.phiGraus);
        cf.E = EnvD("E", cf.E);
        cf.nu = EnvD("NU", cf.nu);
        cf.gamma = EnvD("GAMMA", cf.gamma);
        cf.camposCsv = EnvI("CAMPOS_CSV", cf.camposCsv);
        cf.soCampos = EnvI("SO_CAMPOS", cf.soCampos);
        cf.versao = EnvS("V2_VERSAO", cf.versao);
        return cf;
    }
    /// Verifica os valores (mensagem vazia: tudo certo)
    std::string Erro() const {
        if (caso != "cho1" && caso != "cho2") return "CASO deve ser cho1 ou cho2";
        if (nucleo != "exp" && nucleo != "gauss") return "NUCLEO deve ser exp ou gauss";
        if (campos != 0 && campos != 1) return "CAMPOS deve ser 0 ou 1";
        if (real < 0) return "REAL deve ser >= 0";
        if (!(covC >= 0.) || !(covPhi >= 0.)) return "COV_C e COV_PHI devem ser >= 0";
        if (!(hx > 0.) || !(hz > 0.)) return "LX e LY devem ser > 0";
        if (!(rho > -1. && rho < 1.)) return "RHO deve estar em (-1, 1)";
        if (!(klFrac > 0. && klFrac <= 1.) || klM < 1) return "KL_FRAC em (0, 1] e KL_M >= 1";
        if (!(h > 0.) || h > 5.) return "H_MALHA deve estar em (0, 5]";
        if (porder < 2) return "ORDEM deve ser >= 2";
        if (!(fsTol > 0.) || !(F0 > 0.) || maxitFrozen < 2) return "FS_TOL > 0, FS_F0 > 0, MAXIT_FROZEN >= 2";
        if (fsMax < 0. || (fsMax > 0. && fsMax <= F0)) return "FS_MAX deve ser 0 (sem limite) ou maior que FS_F0";
        if (!(klMargem >= 0.)) return "KL_MARGEM deve ser >= 0";
        if (!(c > 0.) || !(phiGraus >= 0. && phiGraus < 89.) || !(E > 0.) || !(nu > -1. && nu < 0.5) || !(gamma > 0.))
            return "COES > 0, 0 <= PHI_DEG < 89, E > 0, -1 < NU < 0,5 e GAMMA > 0";
        return "";
    }
    /// O caso com o solo e a estatistica desta configuracao
    Caso CasoEfetivo() const {
        Caso cs = CasoPorNome(caso);
        cs.c = c, cs.phiGraus = phiGraus, cs.E = E, cs.nu = nu, cs.gamma = gamma;
        cs.covC = covC, cs.covPhi = covPhi, cs.hx = hx, cs.hz = hz, cs.rho = rho;
        return cs;
    }
};

// ------------------------------------------------------------------------------------------------- campos aleatorios

/// Media, desvio, minimo e maximo ponderados pela area
struct Estat {
    double media = 0., dp = 0., min = 0., max = 0.;
};
inline Estat Estatistica(const std::vector<double> &v, const std::vector<double> &w) {
    Estat e;
    double sw = 0., s1 = 0., s2 = 0.;
    e.min = 1.e300, e.max = -1.e300;
    for (size_t i = 0; i < v.size(); ++i) {
        sw += w[i], s1 += w[i] * v[i], s2 += w[i] * v[i] * v[i];
        e.min = std::min(e.min, v[i]), e.max = std::max(e.max, v[i]);
    }
    e.media = s1 / sw;
    e.dp = std::sqrt(std::max(0., s2 / sw - e.media * e.media));
    return e;
}

/// Gerador dos campos de um caso nos pontos de integracao: KL no retangulo [-m, L + m] x [-m, zc + m] que contem o
/// talude (m = KL_MARGEM; com m = 0,25 m os autovalores coincidem com os das Figs. 3 e 8 de Liao, ver README),
/// lognormal pelas eqs. (5)-(7) e correlacao cruzada pela eq. (12) (Cholesky de [[1, r], [r, 1]] nos coeficientes)
class GeradorCampos {
public:
    std::unique_ptr<v2kl::KLSep2D> kl;
    v2kl::KLSep2D::Pontos P;
    Caso cs;
    std::string nucleo;
    uint64_t semente;

    GeradorCampos(const Caso &c, const std::string &nuc, double fracao, int Mmax, uint64_t sem,
                  const std::vector<double> &x, const std::vector<double> &y, double margem = 0.)
        : cs(c), nucleo(nuc), semente(sem) {
        kl = std::make_unique<v2kl::KLSep2D>(v2kl::CriaKL(nuc, -margem, -margem, cs.L + 2. * margem, cs.zc + 2. * margem,
                                                           cs.hx, cs.hz, fracao, Mmax));
        P = kl->PreparaPontos(x, y);
    }
    int M() const { return kl->M(); }

    /// Normais padrao do campo "campo" (0: c, 1: phi) da realizacao real: sequencia propria por (SEMENTE, REAL, campo)
    std::vector<double> Xi(int real, int campo) const {
        return klsep::Normal(semente * 0x9E3779B97F4A7C15ULL + 3ULL * (uint64_t)real + (uint64_t)campo).Vetor(M());
    }
    /// Campos gaussianos padrao de c e de phi (z_phi vazio se phi nao e aleatorio)
    void Gaussianos(int real, std::vector<double> &zc, std::vector<double> &zphi) const {
        const auto xc = Xi(real, 0);
        zc = kl->Campo(P, xc);
        zphi.clear();
        if (!cs.phiAleatorio) return;
        auto xp = Xi(real, 1);
        const double r = cs.rho, s = std::sqrt(1. - r * r);
        for (size_t k = 0; k < xp.size(); ++k) xp[k] = r * xc[k] + s * xp[k]; // [xi_c, xi_phi] L^T, eq. (12)
        zphi = kl->Campo(P, xp);
    }
    /// c (kPa) e phi (rad) nos pontos
    void Gera(int real, std::vector<double> &c, std::vector<double> &phi) const {
        std::vector<double> zc, zphi;
        Gaussianos(real, zc, zphi);
        const klsep::Lognormal lc(cs.c, cs.covC);
        c.resize(zc.size());
        for (size_t i = 0; i < zc.size(); ++i) c[i] = lc(zc[i]);
        phi.assign(zc.size(), cs.phiGraus * M_PI / 180.);
        if (cs.phiAleatorio && cs.covPhi > 0.) {
            const klsep::Lognormal lp(cs.phiGraus, cs.covPhi);
            for (size_t i = 0; i < zphi.size(); ++i) phi[i] = lp(zphi[i]) * M_PI / 180.;
        }
    }
};

// ------------------------------------------------------------------------------------------------------- realizacoes

/// Resultado de uma realizacao (ou do caso deterministico)
struct Realizacao {
    int real = 0;
    bool campos = false;
    ResultadoSRM srm;
    Estat c, phi; ///< phi em graus
    std::string situacao = "ok";
    double seg = 0.;
};

/// Modelo e gerador de campos montados uma vez; Roda(real) calcula o FS de uma realizacao
class Simulador {
public:
    Config cfg;
    Caso cs;
    std::unique_ptr<Modelo> mod;
    std::unique_ptr<GeradorCampos> ger;
    std::vector<double> c, phi; ///< as propriedades da ultima realizacao

    explicit Simulador(const Config &cf) : cfg(cf), cs(cf.CasoEfetivo()) {
        const std::string e = cfg.Erro();
        if (!e.empty()) throw std::invalid_argument(e);
        mod = std::make_unique<Modelo>(cs, cfg.h, cfg.maxitFrozen, cfg.porder);
    }
    GeradorCampos &Gerador() {
        if (!ger)
            ger = std::make_unique<GeradorCampos>(cs, cfg.nucleo, cfg.klFrac, cfg.klM, cfg.semente, mod->xs, mod->ys,
                                                  cfg.klMargem);
        return *ger;
    }
    Realizacao Roda(int real, bool campos) {
        const auto t0 = Relogio::now();
        Realizacao r;
        r.real = real;
        r.campos = campos;
        if (campos) Gerador().Gera(real, c, phi);
        else {
            c.assign(mod->ips.size(), cs.c);
            phi.assign(mod->ips.size(), cs.phiGraus * M_PI / 180.);
        }
        mod->Propriedades(c, phi);
        std::vector<double> phiG(phi.size());
        for (size_t i = 0; i < phi.size(); ++i) phiG[i] = phi[i] * 180. / M_PI;
        r.c = Estatistica(c, mod->area);
        r.phi = Estatistica(phiG, mod->area);
        r.srm = FatorSRMSeco(*mod->an, cfg.F0, cfg.fsTol, cfg.fsMax);
        if (!(r.srm.FS > 0.)) r.situacao = "sem_equilibrio_geostatico";
        r.seg = Segundos(t0);
        return r;
    }
};

// ------------------------------------------------------------------------------------------- execucao (o contrato)

inline std::string Num(double v) { std::ostringstream o; o << std::setprecision(10) << v; return o.str(); }

/// Grava o conteudo em arq.tmp e renomeia (o arquivo nunca fica pela metade)
inline bool GravaAtomico(const std::filesystem::path &arq, const std::string &conteudo) {
    const std::filesystem::path tmp = arq.string() + ".tmp";
    {
        std::ofstream o(tmp);
        if (!o) return false;
        o << conteudo;
        if (!o) return false;
    }
    std::error_code ec;
    std::filesystem::rename(tmp, arq, ec);
    return !ec;
}

/**
 * Uma execucao do contrato de Monte Carlo na pasta: CAMPOS = 1 roda a realizacao REAL, CAMPOS = 0 o caso
 * deterministico (as medias). Grava fs.csv (t_h = 0, FS_SRM, ...), fs_tentativas.csv (cada solucao da busca do SRM),
 * campos.csv (se CAMPOS_CSV = 1), resumo.csv (cabecalho e uma linha) e, por ultimo, resultado.txt ("chave=valor" com
 * real, FS_min, situacao e seg). Retorna o codigo de saida (0: sucesso).
 */
inline int RodaExecucao(const Config &cfg, const std::filesystem::path &pasta = ".") {
    const auto inicio = Relogio::now();
    const std::string erro = cfg.Erro();
    if (!erro.empty()) {
        std::cerr << "V2_TaludesCho: " << erro << std::endl;
        return 2;
    }
    Simulador sim(cfg);
    int Musado = 0;
    double frac = 0.;
    if (cfg.campos) {
        Musado = sim.Gerador().M();
        frac = sim.Gerador().kl->fracao;
    }
    auto gravaCampos = [&]() {
        const Modelo &m = *sim.mod;
        std::ofstream f(pasta / "campos.csv");
        f << "ip,x,y,area,c,phi_deg\n" << std::setprecision(7);
        for (size_t i = 0; i < m.ips.size(); ++i)
            f << m.ips[i] << ',' << m.xs[i] << ',' << m.ys[i] << ',' << m.area[i] << ',' << sim.c[i] << ','
              << sim.phi[i] * 180. / M_PI << '\n';
    };
    if (cfg.soCampos) { // so os campos: a mesma SEMENTE e REAL refazem o campo da realizacao
        if (cfg.campos) sim.Gerador().Gera(cfg.real, sim.c, sim.phi);
        else {
            sim.c.assign(sim.mod->ips.size(), sim.cs.c);
            sim.phi.assign(sim.mod->ips.size(), sim.cs.phiGraus * M_PI / 180.);
        }
        gravaCampos();
        std::cout << "V2 " << cfg.caso << " real " << cfg.real << ": campos.csv (KL_M = " << Musado << ", fracao "
                  << frac << ")" << std::endl;
        return 0;
    }
    const Realizacao r = sim.Roda(cfg.real, cfg.campos != 0);
    const double seg = Segundos(inicio);
    const Modelo &m = *sim.mod;
    {
        std::ofstream f(pasta / "fs.csv");
        f << "t_h,FS_SRM,FS_parada,F0,tentativas,newton,seg\n" << std::setprecision(10);
        f << 0 << ',' << r.srm.FS << ',' << r.srm.Fparada << ',' << r.srm.F0 << ',' << r.srm.tent.size() << ','
          << r.srm.newton << ',' << r.srm.seg << '\n';
    }
    {
        std::ofstream f(pasta / "fs_tentativas.csv");
        f << "t_h,ordem,fase,g,F,convergiu,newton,nova_tentativa,seg\n" << std::setprecision(10);
        for (size_t k = 0; k < r.srm.tent.size(); ++k) {
            const auto &q = r.srm.tent[k];
            f << 0 << ',' << k + 1 << ',' << q.fase << ',' << q.g << ',' << q.F << ',' << q.conv << ',' << q.iter << ','
              << q.retry << ',' << q.seg << '\n';
        }
    }
    if (cfg.camposCsv) gravaCampos();
    std::string v = cfg.versao;
    std::replace(v.begin(), v.end(), ',', ';');
    std::replace(v.begin(), v.end(), ' ', '_');
    const std::vector<std::pair<std::string, std::string>> res = {
        {"real", std::to_string(cfg.real)}, {"caso", cfg.caso}, {"nucleo", cfg.nucleo},
        {"campos", std::to_string(cfg.campos)}, {"semente", std::to_string(cfg.semente)},
        {"KL_M", std::to_string(Musado)}, {"KL_variancia", cfg.campos ? Num(frac) : ""},
        {"cov_c", Num(cfg.covC)}, {"cov_phi", Num(cfg.covPhi)}, {"Lx", Num(cfg.hx)}, {"Ly", Num(cfg.hz)},
        {"rho", Num(cfg.rho)}, {"h_malha", Num(cfg.h)}, {"elementos", std::to_string(m.NElementos())},
        {"neq", std::to_string(m.an->NEq())}, {"pontos", std::to_string(m.ips.size())},
        {"c_med", Num(r.c.media)}, {"c_dp", Num(r.c.dp)}, {"c_min", Num(r.c.min)}, {"c_max", Num(r.c.max)},
        {"phi_med", Num(r.phi.media)}, {"phi_dp", Num(r.phi.dp)}, {"phi_min", Num(r.phi.min)},
        {"phi_max", Num(r.phi.max)}, {"FS_min", Num(r.srm.FS)}, {"FS_parada", Num(r.srm.Fparada)},
        {"F0", Num(r.srm.F0)}, {"FS_max", Num(cfg.fsMax)}, {"tentativas", std::to_string(r.srm.tent.size())},
        {"newton", std::to_string(r.srm.newton)}, {"situacao", r.situacao}, {"seg", Num(seg)}, {"versao", v}};
    std::ostringstream cab, lin;
    for (size_t i = 0; i < res.size(); ++i) {
        cab << (i ? "," : "") << res[i].first;
        lin << (i ? "," : "") << res[i].second;
    }
    bool ok = GravaAtomico(pasta / "resumo.csv", cab.str() + "\n" + lin.str() + "\n");
    std::ostringstream txt;
    txt << "real=" << cfg.real << " FS_min=" << Num(r.srm.FS) << " situacao=" << r.situacao << " seg=" << Num(seg)
        << " caso=" << cfg.caso << " nucleo=" << cfg.nucleo << " campos=" << cfg.campos << " KL_M=" << Musado
        << " newton=" << r.srm.newton << '\n';
    ok = GravaAtomico(pasta / "resultado.txt", txt.str()) && ok; // por ultimo
    std::cout << "V2 " << cfg.caso << " real " << cfg.real << " (campos " << cfg.campos << ", " << cfg.nucleo
              << "): FS = " << r.srm.FS << ", " << r.situacao << ", " << seg << " s" << std::endl;
    return ok && r.situacao == "ok" ? 0 : 1;
}

// ------------------------------------------------------------------------------------- probabilidade de ruptura

/// Intervalo de Wilson (95%, z = 1,96) para k sucessos em n
inline std::pair<double, double> Wilson(int k, int n, double z = 1.96) {
    if (n <= 0) return {0., 1.};
    const double p = double(k) / n, z2 = z * z, d = 1. + z2 / n;
    const double c = (p + z2 / (2. * n)) / d, m = z * std::sqrt(p * (1. - p) / n + z2 / (4. * n * n)) / d;
    return {std::max(0., c - m), std::min(1., c + m)};
}

} // namespace v2
