#pragma once
// V3_SuperficieLivre: percolacao permanente com superficie livre e face de percolacao (variante iii: esqueleto rigido,
// fluxo de Richards so na agua, ar a pressao atmosferica), contra os exemplos de Sharma, Fujisawa & Murakami (2021),
// Space-Time Finite Element Method for Transient and Unconfined Seepage Flow Analysis, Finite Elem. Anal. Des. 197:103632:
//   4.1 barragem retangular 16 x 24 m com agua a jusante (montante 24 m, jusante 4 m): saida 8,9 m acima do nivel de
//       jusante (y = 12,9 m); a vazao e exatamente a de Dupuit-Charny q = K (H1^2 - H2^2) / (2 L);
//   4.2 barragem de face a jusante inclinada (base 7 m, crista 2 m, montante 5 m, jusante 1 m): saida em y = 3,4 m,
//       face de percolacao de 2,4 m;
//   4.3 barragem retangular 5 x 10 m em dois blocos (k esquerdo = k direito / 10), montante 10 m, jusante 2 m: a linha
//       freatica cruza a interface a 8,05 m e sai 0,46 m acima do nivel de jusante.
// O artigo resolve so a regiao saturada com malha movel (superficie livre nitida). Aqui a malha e fixa e o solo e nao
// saturado: a superficie livre e a isolinha p_w = 0, com uma curva de retencao de Gardner muito ingreme e k_rw = S_e^p
// (pouca agua e pouco fluxo acima da linha freatica), e o regime permanente e alcancado por um transiente (Euler
// implicito com passos crescentes) a partir do estado hidrostatico (o "modelo A" do artigo). A face de percolacao
// (p_w <= 0 com saida de agua) e a logica do passo() do SlopeRandomField (branch UnsaturatedSoil): os nos da face a
// jusante acima do nivel de jusante ficam com p_w = 0 enquanto a agua sai por eles e voltam a fluxo nulo quando a agua
// entraria; um no com fluxo nulo e p_w > 0 passa a p_w = 0. O passo so e aceito quando o conjunto nao muda.
#include "UnsatAnalysis.h"
#include "UnsatProblem.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace v3 {

constexpr REAL gw = 9.81; ///< peso especifico da agua (kN/m3)

/// Um exemplo do artigo: geometria (quadrilatero A(0,0), B(xB,0), D(xD,yD), E(0,yE)), niveis e permeabilidade
struct Caso {
    std::string nome;   ///< retangular (4.1), inclinada (4.2), dois_blocos (4.3)
    std::string secao;  ///< secao do artigo
    REAL xB = 16., xD = 16., yD = 24., yE = 24.;
    REAL H1 = 24., H2 = 4.; ///< niveis de montante (x = 0) e de jusante (m)
    /// condutividade saturada: k = 1e-10 m2 com mu_w = 1e-3 Pa s, rho_w = 1000 kg/m3, g = 9.81 m/s2 (artigo, Sec. 4):
    /// K = k rho_w g / mu_w = 9.81e-4 m/s (conversao nossa; a posicao da superficie livre nao depende de K)
    REAL K = 9.81e-4;
    REAL xInterface = -1.; ///< caso 4.3: interface vertical; a esquerda K * fatorEsq
    REAL fatorEsq = 1.;
    // referencias do artigo (texto)
    REAL ySaidaArtigo = 12.9; ///< cota do ponto de saida na face a jusante (m)
    REAL yInterfaceArtigo = -1.; ///< caso 4.3: cota da linha freatica na interface (m)
    std::string figura;   ///< figura com a linha freatica digitalizada
    std::string csvLinha; ///< arquivo de dados/ com a linha digitalizada
    REAL L() const { return xB; }
    /// a face a jusante B-D na cota y (x da face)
    REAL XFace(REAL y) const { return xB + (xD - xB) * y / yD; }
};

inline Caso CasoRetangular() {
    Caso c;
    c.nome = "retangular", c.secao = "4.1";
    c.xB = 16., c.xD = 16., c.yD = 24., c.yE = 24., c.H1 = 24., c.H2 = 4.;
    c.ySaidaArtigo = 4. + 8.9; // "intersects the downstream face at an elevation of ys = 8.9 m above the tailwater"
    c.figura = "Fig. 7", c.csvLinha = "sharma2021_fig7_retangular.csv";
    return c;
}

inline Caso CasoInclinada() {
    Caso c;
    c.nome = "inclinada", c.secao = "4.2";
    // Fig. 13a: base 7 m, face a montante vertical de 5 m, crista de 2 m, face a jusante reta de D(2,5) a B(7,0)
    c.xB = 7., c.xD = 2., c.yD = 5., c.yE = 5., c.H1 = 5., c.H2 = 1.;
    c.ySaidaArtigo = 3.4; // "intersects the downstream face at y = 3.4 m, and forms a seepage face of height 2.4 m"
    c.figura = "Fig. 14", c.csvLinha = "sharma2021_fig14_inclinada.csv";
    return c;
}

inline Caso CasoDoisBlocos() {
    Caso c;
    c.nome = "dois_blocos", c.secao = "4.3";
    c.xB = 5., c.xD = 5., c.yD = 10., c.yE = 10., c.H1 = 10., c.H2 = 2.;
    c.xInterface = 2.5, c.fatorEsq = 0.1; // k = 1e-11 m2 a esquerda, 1e-10 m2 a direita
    c.ySaidaArtigo = 2. + 0.46;           // Fig. 17: saida 0,46 m acima do nivel de jusante
    c.yInterfaceArtigo = 8.05;            // a linha freatica encontra a interface a 8,05 m da base
    c.figura = "Fig. 17", c.csvLinha = "sharma2021_fig17_dois_blocos.csv";
    return c;
}

inline Caso CasoPorNome(const std::string &n) {
    if (n == "retangular") return CasoRetangular();
    if (n == "inclinada") return CasoInclinada();
    if (n == "dois_blocos") return CasoDoisBlocos();
    std::cerr << "CASO desconhecido: " << n << " (retangular, inclinada ou dois_blocos)\n";
    return Caso{};
}

/// Parametros numericos e do modelo nao saturado (nossos: o artigo nao tem zona nao saturada)
struct Numerica {
    int nx = 32, ny = 48;     ///< divisoes do quadrilatero (ny: a face a jusante tem ny + 1 nos)
    REAL rx = 1.;             ///< razao entre o ultimo e o primeiro intervalo em x (< 1 refina perto da face a jusante)
    // curva de retencao: van Genuchten S_e = [1 + (a p_c)^b]^(-1 + 1/b) ("vg", derivada nula em p_c = 0) ou Gardner
    // S_e = exp(-a p_c) ("gardner", com um bico em p_c = 0); k_rw = S_e^pot. A franja capilar tem espessura da ordem
    // de 1 / (a gamma_w) (0,1 m com a = 1/kPa)
    std::string ret = "vg";
    REAL aG = 1.;             ///< a (1/kPa)
    REAL bVG = 3.;            ///< expoente b do van Genuchten
    REAL pot = 3.;
    REAL krwMin = 1.e-9;      ///< limite inferior de k_rw (condicionamento da zona seca)
    REAL n = 0.3;             ///< porosidade (so muda o pseudo-transiente)
    REAL Kw = 2.2e6;          ///< modulo volumetrico da agua (kPa)
    // pseudo-transiente: passos crescentes ate o regime permanente
    REAL dt0 = 1., fator = 1.5, dtMax = 1.e5, tMax = 1.e7; ///< s
    REAL tolEst = 1.e-4;      ///< regime permanente: |q_entra - q_sai| / q_entra < tolEst
    int passosEstaveis = 3;   ///< ... em tantos passos seguidos, com a mesma face de percolacao
    REAL tolP = 1.e-6;        ///< kPa: no livre da face com p_w > tolP passa a p_w = 0
    REAL tolQ = 1.e-12;       ///< m/s: no com p_w = 0 e entrada de agua > tolQ volta a fluxo nulo
    int maxItFace = 40, maxCortes = 8;
    REAL tolNewtonP = 1.e-6;  ///< maior correcao de p_w no Newton (kPa)
    int maxNewton = 30;
    bool topoFace = true;     ///< o topo (crista) tambem pode ser face de percolacao
    std::string ci = "hidrostatica"; ///< estado inicial: hidrostatica (p_w = gamma_w (H1 - y)) ou linear
    int verbose = 0;
};

/// Malha estruturada de triangulos no quadrilatero A(0,0), B(xB,0), D(xD,yD), E(0,yE) pelo mapa bilinear do quadrado
/// unitario (xi, eta): x = xi (1 - eta) xB + xi eta xD, y = xi eta yD + (1 - xi) eta yE. As linhas eta = cte sao
/// horizontais quando yD = yE (os tres casos). Cada celula vira dois triangulos, com a diagonal alternada (malha
/// "union jack", sem direcao preferencial). Contornos: base AB (EBase), face a jusante BD (ERight), topo DE (ETop),
/// face a montante EA (ELeft). Dominio: id 1, ou 2 a direita da interface do caso 4.3 (alinhada com uma coluna).
inline TPZGeoMesh *MalhaBarragem(const Caso &c, const Numerica &num) {
    std::vector<REAL> xi;
    if (c.xInterface > 0.) {
        const REAL f = c.xInterface / c.xB;
        const int n1 = std::max(1, int(std::lround(num.nx * f)));
        xi = Join({Graded(0., f, n1), Graded(f, 1., std::max(1, num.nx - n1), num.rx)});
    } else
        xi = Graded(0., 1., num.nx, num.rx);
    const std::vector<REAL> eta = Graded(0., 1., num.ny);
    const int nx = int(xi.size()) - 1, ny = num.ny;
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    g->NodeVec().Resize(int64_t(nx + 1) * (ny + 1));
    auto node = [&](int i, int j) { return int64_t(i) * (ny + 1) + j; };
    for (int i = 0; i <= nx; ++i)
        for (int j = 0; j <= ny; ++j) {
            const REAL s = xi[i], e = eta[j];
            TPZManVector<REAL, 3> x = {s * (1. - e) * c.xB + s * e * c.xD, s * e * c.yD + (1. - s) * e * c.yE, 0.};
            g->NodeVec()[node(i, j)] = TPZGeoNode(node(i, j), x, *g);
        }
    auto X = [&](int64_t k, int d) { return g->NodeVec()[k].Coord(d); };
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            const int64_t a = node(i, j), b = node(i + 1, j), cc = node(i + 1, j + 1), d = node(i, j + 1);
            std::array<std::array<int64_t, 3>, 2> t;
            if ((i + j) % 2 == 0) t = {{{a, b, cc}, {a, cc, d}}};
            else t = {{{a, b, d}, {b, cc, d}}};
            for (auto &tri : t) {
                const REAL xc = (X(tri[0], 0) + X(tri[1], 0) + X(tri[2], 0)) / 3.;
                const int id = (c.xInterface > 0. && xc > c.xInterface) ? 2 : 1;
                TPZManVector<int64_t, 3> nn = {tri[0], tri[1], tri[2]};
                new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nn, id, *g);
            }
        }
    auto line = [&](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> nn = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nn, id, *g);
    };
    for (int i = 0; i < nx; ++i) line(node(i, 0), node(i + 1, 0), EBase);
    for (int j = 0; j < ny; ++j) line(node(nx, j), node(nx, j + 1), ERight);
    for (int i = nx; i > 0; --i) line(node(i, ny), node(i - 1, ny), ETop);
    for (int j = ny; j > 0; --j) line(node(0, j), node(0, j - 1), ELeft);
    g->BuildConnectivity();
    return g;
}

/// Parametros do solo da variante iii: esqueleto rigido (nao acoplado), retencao ingreme com k_rw = S_e^pot; os parametros
/// mecanicos nao entram (todos os deslocamentos ficam presos em zero)
inline UnsatParams Solo(const Caso &c, const Numerica &num, REAL fatorK) {
    UnsatParams p;
    p.ks = c.K * fatorK;
    p.n = num.n, p.Kw = num.Kw, p.gw = gw;
    p.coupled = false;
    p.lumped = true;
    if (num.ret == "gardner") {
        p.hydro.fRetention = TPZUnsaturatedHydraulics::EGardner;
        p.hydro.fAG = num.aG;
    } else {
        p.hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        p.hydro.fAVG = num.aG;
        p.hydro.fBVG = num.bVG;
    }
    p.hydro.fSwr = 0.;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
    p.hydro.fPower = num.pot;
    p.hydro.fKrwMin = num.krwMin;
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
    return p;
}

/// Um passo aceito do pseudo-transiente (passos.csv)
struct Passo {
    int passo = 0;
    REAL t = 0., dt = 0.;
    REAL qEsq = 0., qJus = 0., qFace = 0., qTopo = 0.; ///< entrada a montante e saidas (m3/s por metro)
    REAL dVdt = 0.;      ///< variacao do armazenamento (m3/s por metro): soma de R_w / dt em todos os nos
    REAL est = 0.;       ///< |q_entra - q_sai| / q_entra
    int nFace = 0;       ///< nos da face a jusante com p_w = 0 (acima do nivel de jusante)
    REAL ySaida = 0.;    ///< no mais alto da face de percolacao
    int itsFace = 0, newton = 0;
    double seg = 0.;
};

/// No da face a jusante no estado final (face.csv)
struct NoFace {
    REAL x, y, pw, q; ///< q: saida de agua por metro de face (m/s), positiva saindo
    int tipo;         ///< 0: fluxo nulo; 1: face de percolacao (p_w = 0); 2: abaixo do nivel de jusante
};

struct Resultado {
    bool ok = false;
    std::string msg;
    int nos = 0, triangulos = 0, neq = 0;
    int passos = 0, cortes = 0, newton = 0;
    REAL t = 0.;
    bool permanente = false;
    Passo ultimo;
    REAL ySaida = 0., yAcima = 0.; ///< no mais alto da face de percolacao e o no seguinte acima dele
    REAL ySaidaEst() const { return 0.5 * (ySaida + yAcima); }
    REAL qDC = 0.;     ///< Dupuit-Charny K (H1^2 - H2^2) / (2 L) (faces verticais; em serie no caso 4.3)
    REAL qDCcorr = 0.; ///< o mesmo com o termo da zona nao saturada da face a jusante (identidade de Charny com k(p))
    /// limite a priori do acrescimo de vazao da franja capilar na identidade de Charny: (topo - H2) |Phi~_inf| / R, com
    /// |Phi~_inf| = int_0^inf k_rw dp_c / gamma_w (a "altura" da franja) e R = sum(L_i / K_i)
    REAL qFranjaMax = 0.;
    /// caso 4.3: cota em que a linha freatica do bloco esquerdo chega a interface, extrapolada linearmente das duas
    /// ultimas colunas de nos do bloco (x = xI - 2h e xI - h); yInterface0 e a isolinha p_w = 0 na propria interface
    /// (no nosso modelo p_w fica na franja capilar, -1 a 0 kPa, ao longo do trecho em que a linha desce pela interface)
    REAL yInterface = -1., yInterface0 = -1.;
    std::vector<NoFace> face;
    std::vector<std::array<REAL, 3>> nosPw; ///< x, y, p_w de todos os vertices
    std::vector<std::array<int64_t, 3>> tri; ///< vertices dos triangulos
    std::vector<Passo> historico;
    double seg = 0.;
    /// cota da linha freatica (p_w = 0) na vertical x (a mais alta; -1 se nao cruza), com p_w linear em cada triangulo
    REAL Freatica(REAL x) const {
        REAL best = -1.;
        for (const auto &t : tri) {
            const auto &a = nosPw[t[0]], &b = nosPw[t[1]], &c = nosPw[t[2]];
            const REAL xmin = std::min({a[0], b[0], c[0]}), xmax = std::max({a[0], b[0], c[0]});
            if (x < xmin - 1.e-12 || x > xmax + 1.e-12) continue;
            // plano p = p0 + px (x - xa) + py (y - ya)
            const REAL det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
            if (std::fabs(det) < 1.e-14) continue;
            const REAL px = ((b[2] - a[2]) * (c[1] - a[1]) - (c[2] - a[2]) * (b[1] - a[1])) / det;
            const REAL py = ((c[2] - a[2]) * (b[0] - a[0]) - (b[2] - a[2]) * (c[0] - a[0])) / det;
            if (py >= -1.e-12) continue; // a linha freatica tem p_w decrescendo para cima
            const REAL y = a[1] - (a[2] + px * (x - a[0])) / py;
            // dentro do triangulo (coordenadas baricentricas)
            const REAL l1 = ((b[0] - x) * (c[1] - y) - (c[0] - x) * (b[1] - y)) / det;
            const REAL l2 = ((c[0] - x) * (a[1] - y) - (a[0] - x) * (c[1] - y)) / det;
            const REAL l3 = 1. - l1 - l2;
            const REAL tol = 1.e-9;
            if (l1 < -tol || l2 < -tol || l3 < -tol) continue;
            best = std::max(best, y);
        }
        return best;
    }
};

using Relogio = std::chrono::steady_clock;
inline double Segundos(Relogio::time_point t0) {
    return std::chrono::duration<double>(Relogio::now() - t0).count();
}

/**
 * Roda um caso ate o regime permanente. Contornos da agua: x = 0 com p_w = gamma_w (H1 - y); face a jusante com
 * p_w = gamma_w (H2 - y) ate o nivel de jusante e face de percolacao acima dele (e no topo, se topoFace); base
 * impermeavel. Ar: p_a = 0 em todos os nos (variante iii); deslocamentos: todos nulos (esqueleto rigido). Cada passo
 * aceito vai para passos (se nao nulo).
 */
inline Resultado Executa(const Caso &c, const Numerica &num, std::ostream *passos = nullptr) {
    const auto t0 = Relogio::now();
    Resultado res;
    std::map<int, UnsatParams> params = {{1, Solo(c, num, c.xInterface > 0. ? c.fatorEsq : 1.)}};
    if (c.xInterface > 0.) params[2] = Solo(c, num, 1.);
    UnsatProblem pb(MalhaBarragem(c, num), params, {});
    UnsatAnalysis &an = pb.A();
    an.SetTolerances(1.e-7, num.tolNewtonP, num.maxNewton, 25);
    an.SetVerbose(num.verbose > 1);
    res.neq = int(an.NEq());
    res.nos = int(an.WaterEq().size());
    for (int64_t i = 0; i < pb.gmesh->NElements(); ++i) {
        TPZGeoEl *gel = pb.gmesh->Element(i);
        if (gel && gel->Dimension() == 2) {
            res.triangulos++;
            res.tri.push_back({gel->NodeIndex(0), gel->NodeIndex(1), gel->NodeIndex(2)});
        }
    }
    const REAL tolY = 1.e-9 * c.yE;

    // ---- contornos fixos: deslocamentos (esqueleto rigido), ar e as cargas de montante e de jusante
    auto contorno = an.AllDisplacements(an.X());
    for (const auto &e : an.AirEq()) contorno[e.second] = 0.;
    std::set<int64_t> eqEsq, eqJus;
    for (const auto &b : an.BoundaryNodes(ELeft)) {
        const REAL y = an.Coord(b.first, 1);
        if (y <= c.H1 + tolY) {
            contorno[an.WaterEq().at(b.first)] = gw * (c.H1 - y);
            eqEsq.insert(an.WaterEq().at(b.first));
        }
    }
    const auto bordaJus = an.BoundaryNodes(ERight);
    for (const auto &b : bordaJus) {
        const REAL y = an.Coord(b.first, 1);
        if (y <= c.H2 + tolY && !eqEsq.count(an.WaterEq().at(b.first))) {
            contorno[an.WaterEq().at(b.first)] = gw * (c.H2 - y);
            eqJus.insert(an.WaterEq().at(b.first));
        }
    }
    // candidatos a face de percolacao (no -> comprimento tributario): face a jusante acima do nivel de jusante e topo
    std::map<int64_t, REAL> cand;
    std::set<int64_t> nosFace;
    for (const auto &b : bordaJus)
        if (!contorno.count(an.WaterEq().at(b.first))) cand[b.first] += b.second, nosFace.insert(b.first);
    if (num.topoFace)
        for (const auto &b : an.BoundaryNodes(ETop))
            if (!contorno.count(an.WaterEq().at(b.first))) cand[b.first] += b.second;

    // ---- estado inicial: hidrostatico com o nivel de montante em todo o dominio (modelo A do artigo) ou linear
    TPZFMatrix<STATE> x0(an.NEq(), 1, 0.);
    for (const auto &e : an.WaterEq()) {
        const REAL xx = an.Coord(e.first, 0), y = an.Coord(e.first, 1);
        const REAL h = num.ci == "linear" ? c.H1 + (c.H2 - c.H1) * xx / c.L() : c.H1;
        x0(e.second, 0) = gw * (h - y);
    }
    for (const auto &f : contorno) x0(f.first, 0) = f.second;
    an.X() = x0;
    an.Restore();

    // ---- um passo de dt com a face de percolacao (adaptado do passo() do SlopeRandomField, sem chuva)
    std::set<int64_t> saida;
    int nPassos = 0;
    auto vazoes = [&](const TPZFMatrix<STATE> &R, REAL dt, const std::set<int64_t> &face, Passo &ps) {
        ps.qEsq = ps.qJus = ps.qFace = ps.qTopo = ps.dVdt = 0.;
        for (auto e : eqEsq) ps.qEsq += R.GetVal(e, 0) / dt;
        for (auto e : eqJus) ps.qJus -= R.GetVal(e, 0) / dt;
        for (auto n : face) (nosFace.count(n) ? ps.qFace : ps.qTopo) -= R.GetVal(an.WaterEq().at(n), 0) / dt;
        for (const auto &e : an.WaterEq()) ps.dVdt += R.GetVal(e.second, 0) / dt;
        ps.est = std::fabs(ps.qEsq - ps.qJus - ps.qFace - ps.qTopo) / std::max(std::fabs(ps.qEsq), REAL(1.e-30));
        ps.nFace = 0;
        ps.ySaida = c.H2;
        for (auto n : face)
            if (nosFace.count(n)) ps.nFace++, ps.ySaida = std::max(ps.ySaida, an.Coord(n, 1));
    };
    Passo ultimo;
    auto passo = [&](REAL t, REAL dt) -> bool {
        const auto r0 = Relogio::now();
        auto tentativa = saida;
        for (const auto &b : cand)
            if (an.X().GetVal(an.WaterEq().at(b.first), 0) > num.tolP) tentativa.insert(b.first);
        int newton = 0;
        for (int it = 0; it < num.maxItFace; ++it) {
            auto fixed = contorno;
            for (auto nd : tentativa) fixed[an.WaterEq().at(nd)] = 0.;
            const auto sr = an.Solve(dt, fixed, {}, false);
            newton += sr.iter;
            res.newton += sr.iter;
            if (!sr.conv) break;
            const auto R = an.FullResidual(); // R_w / dt: entrada de agua nos nos com p_w imposto
            const TPZFMatrix<STATE> sol = pb.mphys->Solution();
            auto nova = tentativa;
            for (const auto &b : cand) {
                const int64_t eq = an.WaterEq().at(b.first);
                if (tentativa.count(b.first)) {
                    if (R.GetVal(eq, 0) / (dt * b.second) > num.tolQ) nova.erase(b.first);
                } else if (sol.GetVal(eq, 0) > num.tolP)
                    nova.insert(b.first);
            }
            if (nova == tentativa) {
                Passo ps;
                ps.passo = ++nPassos;
                ps.t = t + dt, ps.dt = dt;
                vazoes(R, dt, tentativa, ps);
                ps.itsFace = it + 1, ps.newton = newton, ps.seg = Segundos(r0);
                an.Accept();
                saida = tentativa;
                ultimo = ps;
                res.historico.push_back(ps);
                if (passos)
                    *passos << ps.passo << ',' << ps.t << ',' << ps.dt << ',' << ps.qEsq << ',' << ps.qJus << ','
                            << ps.qFace << ',' << ps.qTopo << ',' << ps.dVdt << ',' << ps.est << ',' << ps.nFace << ','
                            << ps.ySaida << ',' << ps.itsFace << ',' << ps.newton << ',' << ps.seg << std::endl;
                if (num.verbose)
                    std::cout << "passo " << ps.passo << " t = " << ps.t << " s dt = " << dt << " q_esq = " << ps.qEsq
                              << " q_sai = " << ps.qJus + ps.qFace + ps.qTopo << " est = " << ps.est
                              << " face = " << ps.nFace << " y_saida = " << ps.ySaida << " (" << ps.itsFace
                              << " tentativas, " << newton << " Newton)" << std::endl;
                return true;
            }
            tentativa = nova;
        }
        an.Restore();
        return false;
    };
    std::function<bool(REAL, REAL, int)> avanca = [&](REAL t, REAL dt, int nivel) -> bool {
        if (passo(t, dt)) return true;
        if (nivel >= num.maxCortes) return false;
        res.cortes++;
        return avanca(t, 0.5 * dt, nivel + 1) && avanca(t + 0.5 * dt, 0.5 * dt, nivel + 1);
    };

    // ---- pseudo-transiente ate o regime permanente
    REAL t = 0., dt = num.dt0;
    int estaveis = 0;
    std::set<int64_t> faceAnt;
    while (t < num.tMax - 1.e-9) {
        dt = std::min(dt, num.tMax - t);
        const int cortes0 = res.cortes;
        if (!avanca(t, dt, 0)) {
            res.msg = "passo sem convergencia em t = " + std::to_string(t) + " s";
            break;
        }
        t += dt;
        estaveis = (ultimo.est < num.tolEst && saida == faceAnt) ? estaveis + 1 : 0;
        faceAnt = saida;
        if (estaveis >= num.passosEstaveis) {
            res.permanente = true;
            break;
        }
        // depois de um corte, o passo volta a crescer a partir do subpasso que convergiu
        if (res.cortes > cortes0) dt = res.historico.back().dt;
        dt = std::min(dt * num.fator, num.dtMax);
    }
    res.t = t;
    res.passos = nPassos;
    res.ultimo = ultimo;
    res.ok = res.permanente;
    if (!res.permanente && res.msg.empty()) res.msg = "tMax sem regime permanente";

    // ---- estado final: p_w nos vertices, face a jusante, ponto de saida, Dupuit-Charny
    const TPZFMatrix<STATE> &X = an.X();
    res.nosPw.assign(pb.gmesh->NNodes(), {0., 0., 0.});
    for (const auto &e : an.WaterEq())
        res.nosPw[e.first] = {an.Coord(e.first, 0), an.Coord(e.first, 1), X.GetVal(e.second, 0)};
    {
        // vazao de cada no da face (residuo do estado aceito, com o ultimo dt)
        an.LoadSolution(X);
        const auto R = an.FullResidual();
        const REAL dtf = ultimo.dt > 0. ? ultimo.dt : 1.;
        for (const auto &b : bordaJus) {
            const int64_t eq = an.WaterEq().at(b.first);
            NoFace nf;
            nf.x = an.Coord(b.first, 0), nf.y = an.Coord(b.first, 1), nf.pw = X.GetVal(eq, 0);
            nf.tipo = eqJus.count(eq) ? 2 : (saida.count(b.first) ? 1 : 0);
            nf.q = nf.tipo ? -R.GetVal(eq, 0) / (dtf * b.second) : 0.;
            res.face.push_back(nf);
        }
        std::sort(res.face.begin(), res.face.end(), [](const NoFace &a, const NoFace &b) { return a.y < b.y; });
    }
    res.ySaida = c.H2;
    for (const auto &f : res.face)
        if (f.tipo >= 1) res.ySaida = std::max(res.ySaida, f.y);
    res.yAcima = res.ySaida;
    for (const auto &f : res.face)
        if (f.y > res.ySaida + tolY) {
            res.yAcima = f.y;
            break;
        }
    if (c.xD == c.xB) {
        // blocos em serie (caso 4.3): em cada bloco q = -K_i dF/dx com F(x) = int Phi~ dy, Phi~ = int_0^P k_rw ds; F e
        // continua na interface, entao q sum(L_i / K_i) = F(0) - F(L) e q = (H1^2 - H2^2) / (2 sum(L_i / K_i))
        const REAL RL = c.xInterface > 0. ? c.xInterface / (c.K * c.fatorEsq) + (c.L() - c.xInterface) / c.K
                                          : c.L() / c.K;
        res.qDC = (c.H1 * c.H1 - c.H2 * c.H2) / (2. * RL);
        // identidade de Charny com a permeabilidade k(P) (P = p_w / gamma_w): q L = int Phi(P) dy na face a montante
        // menos a mesma integral na face a jusante, com o potencial de Kirchhoff Phi(P) = int_0^P k(s) ds (topo e base
        // sem fluxo). Montante: P = H1 - y >= 0 (Phi = K P); jusante: K H2^2 / 2 abaixo do nivel, 0 na face de
        // percolacao e Phi < 0 acima da saida (succao), que aumenta a vazao. Phi pela regra de Simpson com a curva do
        // material (o limite krwMin incluido); a integral na altura pela regra dos trapezios nos nos da face
        const TPZUnsaturatedHydraulics h = Solo(c, num, 1.).hydro;
        auto phi = [&](REAL pw) {
            if (pw >= 0.) return REAL(0.);
            const int m = 400; // p_c de 0 a -pw
            const REAL pc = -pw, d = pc / m;
            REAL soma = 0.;
            for (int k = 0; k <= m; ++k) {
                REAL S, G, kr, dkr;
                h.Retention(k * d, S, G);
                h.WaterPermeability(S, kr, dkr);
                soma += kr * ((k == 0 || k == m) ? 1. : (k % 2 ? 4. : 2.));
            }
            return -soma * d / 3. / gw; // Phi~ (m): dP = -dp_c / gamma_w
        };
        REAL extra = 0.;
        for (size_t i = 1; i < res.face.size(); ++i)
            extra -= 0.5 * (phi(res.face[i].pw) + phi(res.face[i - 1].pw)) * (res.face[i].y - res.face[i - 1].y);
        res.qDCcorr = res.qDC + extra / RL;
        res.qFranjaMax = (c.yE - c.H2) * (-phi(-50. / num.aG)) / RL;
    }
    if (c.xInterface > 0.) {
        const int n1 = std::max(1, int(std::lround(num.nx * c.xInterface / c.xB)));
        const REAL h = c.xInterface / n1;
        res.yInterface0 = res.Freatica(c.xInterface);
        res.yInterface = 2. * res.Freatica(c.xInterface - h) - res.Freatica(c.xInterface - 2. * h);
    }
    res.seg = Segundos(t0);
    return res;
}

/// Leitura de uma linha digitalizada (dados/<arquivo>: colunas x_m,y_m; linhas com # sao comentarios)
inline std::vector<std::array<REAL, 2>> LeLinha(const std::string &arquivo) {
    std::vector<std::array<REAL, 2>> v;
    std::ifstream in(arquivo);
    std::string s;
    while (std::getline(in, s)) {
        if (s.empty() || s[0] == '#' || !(std::isdigit(static_cast<unsigned char>(s[0])) || s[0] == '-' || s[0] == '.'))
            continue;
        std::replace(s.begin(), s.end(), ',', ' ');
        std::istringstream ss(s);
        REAL x, y;
        if (ss >> x >> y) v.push_back({x, y});
    }
    return v;
}

/// Interpolacao linear numa linha (x crescente); fora do intervalo, -1
inline REAL Interpola(const std::vector<std::array<REAL, 2>> &v, REAL x) {
    for (size_t i = 1; i < v.size(); ++i)
        if (x >= v[i - 1][0] - 1.e-12 && x <= v[i][0] + 1.e-12) {
            const REAL f = (x - v[i - 1][0]) / (v[i][0] - v[i - 1][0]);
            return v[i - 1][1] + f * (v[i][1] - v[i - 1][1]);
        }
    return -1.;
}

} // namespace v3
