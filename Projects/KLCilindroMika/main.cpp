// KLCilindroMika - expansao de Karhunen-Loeve de Galerkin de um campo aleatorio com nucleo exponencial no meio cilindro
// vazado de Mika, Hughes, Schillinger, Wriggers e Hiemstra (2021, CMAME 379:113730), Secao 6.2, Exemplo 1-1:
// D = {8 <= rho <= 10, 0 <= theta <= pi, 0 <= z <= 15}, G(x,x') = sigma^2 exp(-|x - x'| / (b L)), b = 0.5, L = 10,
// sigma^2 = 1. O problema de autovalores C v = lambda B v usa o que o NeoPZ ja tem: o material TPZKarhunenLoeveMat
// (nucleo do tipo 4, exp(-dist/Lx)) e TPZInterpolationSpace::CalcStiffB (B_ij = int N_i N_j) e CalcStiffC
// (C_ij = int int N_i(x) G(x,x') N_j(x') dx dx', um par de elementos). B e C sao densas e montadas percorrendo todos os
// pares de elementos; os autovalores saem do Eigen (Cholesky de B e autovalores da matriz simetrica L^-1 C L^-T).
//
// Malha: hexaedros H1 de ordem P em theta x rho x z, com a geometria exata do cilindro (TPZGeoBlend com TPZArc3D em
// todas as arestas circunferenciais). Parametros (variaveis de ambiente):
//   P=2                 ordem dos elementos (1 ou 2)
//   MALHAS=8x1x4,...    malhas nt x nr x nz (elementos em theta, rho e z), resolvidas em sequencia
//   NG=0                pontos de Gauss por direcao nos pares de elementos sem no em comum (0: regra do elemento, P+1)
//   NG_VIZINHO=0        idem, nos pares de elementos distintos com no em comum (0: o mesmo de NG)
//   NG_PROPRIO=0        idem, no par de um elemento com ele mesmo (0: o mesmo de NG)
//   THREADS=2           threads na montagem de C
//   VALIDAR=0           1: confere a montagem com CalcStiffC par a par, o algoritmo anterior do CalcStiffC e o caso 2D
//   CSV=                arquivo CSV onde acrescentar uma linha por malha (autovalores, soma, tempos)
// Ver README.md.
#include <Eigen/Dense>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"
#include "TPZGeoCube.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "tpzarc3d.h"
#include "pzinterpolationspace.h"
#include "TPZElementMatrixT.h"
#include "TPZMaterialDataT.h"
#include "TPZKarhunenLoeveMat.h"
#include "pzquad.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

// ---- problema (Mika et al. 2021, Secao 6.2) ----
constexpr REAL kRaioInterno = 8., kRaioExterno = 10., kAltura = 15.;
constexpr REAL kB = 0.5, kL = 10.; // comprimento de correlacao b L = 5
// sigma^2 = 1: o nucleo do TPZKarhunenLoeveMat nao tem a variancia, e os autovalores sao os de sigma^2 = 1
const REAL kVolume = 0.5 * M_PI * (kRaioExterno * kRaioExterno - kRaioInterno * kRaioInterno) * kAltura; // 270 pi
constexpr int kMatVolume = 1, kMatArco = 2;
constexpr int kNModos = 20;

// Tabela 4 de Mika et al. (2021): Galerkin isogeometrico com o nucleo exato, NURBS p = (2,2,2), quadratura de Gauss
// (p+1)^3. Galerkin 2: 3800 elementos e 6642 graus de liberdade; Galerkin 1: 256 elementos e 1050.
const double kGalerkin2[kNModos] = {162.7965791, 91.42804062, 57.56447741, 51.08762278, 38.79483423, 27.90128438,
    25.05356161, 19.36659412, 16.15443088, 15.79273209, 15.14385016, 11.21090778, 10.17639037, 9.690982310,
    8.052352020, 7.576621410, 6.719925970, 6.443915210, 6.174771170, 5.761319370};
const double kGalerkin1[kNModos] = {162.7991539, 91.43063070, 57.56702901, 51.09017918, 38.79740931, 27.90386143,
    25.05611145, 19.36893419, 16.15700369, 15.79530798, 15.14622914, 11.21328690, 10.17896939, 9.693564320,
    8.054783080, 7.578085990, 6.722378350, 6.445576690, 6.177345410, 5.763786920};

int envi(const char *k, int d) { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }

using Relogio = std::chrono::steady_clock;
double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }

// pico da memoria residente do processo (VmHWM), em MB
double PicoMemoriaMB() {
    std::ifstream f("/proc/self/status");
    std::string l;
    while (std::getline(f, l))
        if (l.rfind("VmHWM:", 0) == 0) return std::atof(l.c_str() + 6) / 1024.;
    return -1.;
}

struct Malha { int nt = 0, nr = 0, nz = 0; };

// "8x1x4,12x2x6" -> malhas
std::vector<Malha> LeMalhas(const std::string &s) {
    std::vector<Malha> m;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) {
        if (t.empty()) continue;
        Malha a;
        if (std::sscanf(t.c_str(), "%dx%dx%d", &a.nt, &a.nr, &a.nz) != 3 || a.nt < 1 || a.nr < 1 || a.nz < 1) {
            std::cerr << "malha invalida em MALHAS: " << t << " (use nt x nr x nz, por exemplo 16x2x8)\n";
            std::exit(2);
        }
        m.push_back(a);
    }
    return m;
}

// regra de Gauss-Legendre do hexaedro com n pontos por direcao (a ordem 2n - 1 e a maior integrada exatamente)
TPZIntCube3D RegraCubo(int n) {
    TPZIntCube3D r(2 * n - 1);
    if (r.NPoints() != n * n * n) {
        std::cerr << "regra de " << n << " pontos por direcao indisponivel\n";
        std::exit(2);
    }
    return r;
}

// ---- malha geometrica ----
// No (i, j, k): theta_i = pi i / nt, rho_j = r + (R - r) j / nr, z_k = H k / nz. Hexaedro (i, j, k) com xi -> rho,
// eta -> theta e zeta -> z (jacobiano positivo). Cada aresta circunferencial e um TPZArc3D (com um no no meio do arco),
// vizinho dos hexaedros TPZGeoBlend que a contem: o mapa do hexaedro fica rho(xi) (cos theta(eta), sin theta(eta)) +
// z(zeta) e_z, o do cilindro, com rho, theta e z lineares nas coordenadas de referencia (conferido em VerificaGeometria).
TPZGeoMesh *MalhaMeioCilindro(const Malha &m, std::vector<std::array<int, 3>> &celula) {
    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    const REAL dr = (kRaioExterno - kRaioInterno) / m.nr, dt = M_PI / m.nt, dz = kAltura / m.nz;
    auto novoNo = [&](REAL rho, int metadesTheta, REAL z) -> int64_t {
        // theta = metadesTheta * dt / 2; nos extremos, cos e sin exatos (y = 0 nas faces planas theta = 0 e pi)
        const REAL th = 0.5 * dt * metadesTheta;
        REAL c = std::cos(th), s = std::sin(th);
        if (metadesTheta == 0) { c = 1.; s = 0.; }
        if (metadesTheta == 2 * m.nt) { c = -1.; s = 0.; }
        const int64_t id = gmesh->NodeVec().AllocateNewElement();
        TPZManVector<REAL, 3> co = {rho * c, rho * s, z};
        gmesh->NodeVec()[id].Initialize(co, *gmesh);
        return id;
    };
    auto rho = [&](int j) { return kRaioInterno + dr * j; };
    auto zz = [&](int k) { return dz * k; };
    std::vector<int64_t> no((size_t)(m.nt + 1) * (m.nr + 1) * (m.nz + 1));
    auto id = [&](int i, int j, int k) { return no[(size_t)i + (m.nt + 1) * ((size_t)j + (m.nr + 1) * (size_t)k)]; };
    for (int k = 0; k <= m.nz; k++)
        for (int j = 0; j <= m.nr; j++)
            for (int i = 0; i <= m.nt; i++)
                no[(size_t)i + (m.nt + 1) * ((size_t)j + (m.nr + 1) * (size_t)k)] = novoNo(rho(j), 2 * i, zz(k));

    using Hexa = TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>>;
    using Arco = TPZGeoElRefPattern<pzgeom::TPZArc3D>;
    celula.assign(0, {0, 0, 0});
    for (int k = 0; k < m.nz; k++)
        for (int j = 0; j < m.nr; j++)
            for (int i = 0; i < m.nt; i++) {
                TPZManVector<int64_t, 8> nos = {id(i, j, k), id(i, j + 1, k), id(i + 1, j + 1, k), id(i + 1, j, k),
                    id(i, j, k + 1), id(i, j + 1, k + 1), id(i + 1, j + 1, k + 1), id(i + 1, j, k + 1)};
                auto *gel = new Hexa(nos, kMatVolume, *gmesh);
                if ((int64_t)celula.size() <= gel->Index()) celula.resize(gel->Index() + 1, {-1, -1, -1});
                celula[gel->Index()] = {i, j, k};
            }
    // arcos: nos (inicio, fim, ponto no meio do arco)
    for (int k = 0; k <= m.nz; k++)
        for (int j = 0; j <= m.nr; j++)
            for (int i = 0; i < m.nt; i++) {
                TPZManVector<int64_t, 3> nos = {id(i, j, k), id(i + 1, j, k), novoNo(rho(j), 2 * i + 1, zz(k))};
                new Arco(nos, kMatArco, *gmesh);
            }
    gmesh->BuildConnectivity();
    return gmesh;
}

// ---- dados de cada elemento ----
struct ElementoKL {
    TPZInterpolationSpace *intel = nullptr;
    std::vector<int64_t> eq;                 // equacao global de cada funcao de forma (na ordem das linhas de phi)
    std::array<int, 3> celula = {0, 0, 0};    // (i, j, k)
    std::vector<int> toca;                    // elementos distintos com algum no em comum
    // dados nos pontos das regras dos pares sem no em comum, com no em comum e do elemento com ele mesmo; quando duas
    // regras sao iguais, os dados sao calculados uma vez so (pontosVizinho e pontosProprio ficam vazios)
    TPZKLPointData pontos, pontosVizinho, pontosProprio;
    int fonteVizinho = 0, fonteProprio = 0; // 0: pontos, 1: pontosVizinho, 2: pontosProprio
    const TPZKLPointData &Fonte(int f) const { return f == 0 ? pontos : (f == 1 ? pontosVizinho : pontosProprio); }
    const TPZKLPointData &Vizinho() const { return Fonte(fonteVizinho); }
    const TPZKLPointData &Proprio() const { return Fonte(fonteProprio); }
};

std::vector<int64_t> Equacoes(TPZCompMesh *cmesh, TPZCompEl *cel) {
    std::vector<int64_t> eq;
    for (int ic = 0; ic < cel->NConnects(); ic++) {
        TPZConnect &c = cel->Connect(ic);
        if (c.HasDependency() || c.IsCondensed()) {
            std::cerr << "conector com dependencia ou condensado: a malha devia ser conforme\n";
            std::exit(1);
        }
        const int64_t seq = c.SequenceNumber();
        const int64_t pos = cmesh->Block().Position(seq);
        const int sz = cmesh->Block().Size(seq);
        for (int k = 0; k < sz; k++) eq.push_back(pos + k);
    }
    return eq;
}

// C(eq1, eq2) += ce e, se transposta, C(eq2, eq1) += ce^T
void Espalha(const TPZFMatrix<STATE> &ce, const std::vector<int64_t> &eq1, const std::vector<int64_t> &eq2,
             bool transposta, Eigen::MatrixXd &C) {
    const int64_t n1 = (int64_t)eq1.size(), n2 = (int64_t)eq2.size();
    if (n1 == 0 || n2 == 0) return;
    if (ce.Rows() != n1 || ce.Cols() != n2) {
        std::cerr << "matriz do par com dimensoes diferentes das equacoes dos elementos\n";
        std::exit(1);
    }
    const STATE *v = &ce.g(0, 0); // por colunas
    for (int64_t j = 0; j < n2; j++) {
        const int64_t gj = eq2[j];
        double *col = C.data() + gj * C.rows();
        for (int64_t i = 0; i < n1; i++) col[eq1[i]] += v[i + j * n1];
    }
    if (!transposta) return;
    for (int64_t i = 0; i < n1; i++) {
        double *col = C.data() + eq1[i] * C.rows();
        for (int64_t j = 0; j < n2; j++) col[eq2[j]] += v[i + j * n1];
    }
}

// mapa exato da celula (i, j, k) no ponto de referencia qsi e |det J| correspondente
void MapaExato(const Malha &m, const std::array<int, 3> &c, const TPZVec<REAL> &qsi, REAL x[3], REAL &detjac) {
    const REAL dr = (kRaioExterno - kRaioInterno) / m.nr, dt = M_PI / m.nt, dz = kAltura / m.nz;
    const REAL rho = kRaioInterno + dr * (c[1] + 0.5 * (qsi[0] + 1.));
    const REAL th = dt * (c[0] + 0.5 * (qsi[1] + 1.));
    const REAL z = dz * (c[2] + 0.5 * (qsi[2] + 1.));
    x[0] = rho * std::cos(th);
    x[1] = rho * std::sin(th);
    x[2] = z;
    detjac = rho * dr * dt * dz / 8.;
}

// o CalcStiffC antes desta correcao: a regra do primeiro elemento tambem no segundo, e os dados do segundo elemento
// calculados de novo para cada par de pontos; usado so para comparar com a versao corrigida em elementos iguais
void CalcStiffCAnterior(TPZInterpolationSpace *el1, TPZInterpolationSpace *el2, TPZKarhunenLoeveMat *mat,
                        TPZFMatrix<STATE> &ce) {
    TPZMaterialDataT<REAL> data1, data2;
    el1->InitMaterialData(data1);
    el2->InitMaterialData(data2);
    data1.fNeedsSol = data2.fNeedsSol = false; // a solucao nao entra em C
    const int dim = el1->Dimension();
    TPZManVector<REAL, 3> pt1(dim, 0.), pt2(dim, 0.);
    REAL w1 = 0., w2 = 0.;
    TPZAutoPointer<TPZIntPoints> regra = el1->GetIntegrationRule().Clone();
    const int np = regra->NPoints();
    ce.Redim(el1->NShapeF(), el1->NShapeF());
    TPZFMatrix<STATE> cet;
    for (int ip = 0; ip < np; ip++) {
        regra->Point(ip, pt1, w1);
        el1->ComputeRequiredData(data1, pt1);
        w1 *= std::fabs(data1.detjac);
        for (int jp = 0; jp < np; jp++) {
            regra->Point(jp, pt2, w2);
            el2->ComputeRequiredData(data2, pt2);
            w2 *= std::fabs(data2.detjac);
            mat->ContributeC(data1, data2, w1, w2, cet);
            ce += cet;
        }
    }
}

// max |a - b| / max |b|
double DifRel(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) {
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) return 1.e300;
    double d = 0., n = 0.;
    for (int64_t i = 0; i < a.Rows(); i++)
        for (int64_t j = 0; j < a.Cols(); j++) {
            d = std::max(d, std::fabs(a.GetVal(i, j) - b.GetVal(i, j)));
            n = std::max(n, std::fabs(b.GetVal(i, j)));
        }
    return n > 0. ? d / n : d;
}

double DifRelTransposta(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) {
    TPZFMatrix<STATE> bt;
    b.Transpose(&bt);
    return DifRel(a, bt);
}

// Caso 2D (o uso anterior do CalcStiffC): malha de quadrilateros de ordem 2 e uma malha com um triangulo e um
// quadrilatero. Nos quadrilateros, o CalcStiffC corrigido da o mesmo C do algoritmo anterior. No par triangulo x
// quadrilatero, com Lx = 1e9 (G = 1 a menos de 1e-8), C_ij = (int N_i)(int N_j), conferido com a regra de cada
// elemento; o algoritmo anterior usava a regra do triangulo no quadrilatero.
void ValidaCaso2D() {
    for (int caso = 0; caso < 2; caso++) {
        auto *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(2);
        auto no = [&](REAL x, REAL y) {
            const int64_t id = gmesh->NodeVec().AllocateNewElement();
            TPZManVector<REAL, 3> co = {x, y, 0.};
            gmesh->NodeVec()[id].Initialize(co, *gmesh);
            return id;
        };
        if (caso == 0) { // 3 x 2 quadrilateros em [0,3] x [0,2], com um no deslocado (elementos nao afins)
            std::vector<int64_t> n;
            for (int j = 0; j <= 2; j++)
                for (int i = 0; i <= 3; i++) n.push_back(no(i + (i == 1 && j == 1 ? 0.2 : 0.), j + (i == 2 && j == 1 ? -0.15 : 0.)));
            for (int j = 0; j < 2; j++)
                for (int i = 0; i < 3; i++) {
                    TPZManVector<int64_t, 4> q = {n[i + 4 * j], n[i + 1 + 4 * j], n[i + 1 + 4 * (j + 1)], n[i + 4 * (j + 1)]};
                    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(q, kMatVolume, *gmesh);
                }
        } else { // um quadrilatero [0,1]^2 e um triangulo (1,0), (2,0), (1,1)
            const int64_t a = no(0, 0), b = no(1, 0), c = no(1, 1), d = no(0, 1), e = no(2, 0);
            TPZManVector<int64_t, 4> q = {a, b, c, d};
            TPZManVector<int64_t, 3> t = {b, e, c};
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(q, kMatVolume, *gmesh);
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(t, kMatVolume, *gmesh);
        }
        gmesh->BuildConnectivity();
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(2);
        cmesh->SetDimModel(2);
        cmesh->SetAllCreateFunctionsContinuous();
        const REAL lx = caso == 0 ? 1.5 : 1.e9;
        auto *mat = new TPZKarhunenLoeveMat(kMatVolume, lx, lx, lx, 2, 4);
        cmesh->InsertMaterialObject(mat);
        cmesh->AutoBuild();
        std::vector<TPZInterpolationSpace *> els;
        for (auto *cel : cmesh->ElementVec())
            if (cel) els.push_back(dynamic_cast<TPZInterpolationSpace *>(cel));
        double difAnterior = 0., difProduto = 0., difSim = 0.;
        for (auto *e1 : els)
            for (auto *e2 : els) {
                TPZElementMatrixT<REAL> ce, ce21;
                e1->CalcStiffC(e2, ce);
                e2->CalcStiffC(e1, ce21);
                difSim = std::max(difSim, DifRelTransposta(ce.fMat, ce21.fMat));
                if (caso == 0) {
                    TPZFNMatrix<81, STATE> ant;
                    CalcStiffCAnterior(e1, e2, mat, ant);
                    difAnterior = std::max(difAnterior, DifRel(ce.fMat, ant));
                } else {
                    TPZKLPointData d1, d2;
                    e1->KLPointData(e1->GetIntegrationRule(), d1);
                    e2->KLPointData(e2->GetIntegrationRule(), d2);
                    TPZFNMatrix<81, STATE> prod(d1.fPhi.Rows(), d2.fPhi.Rows(), 0.);
                    for (int64_t i = 0; i < d1.fPhi.Rows(); i++)
                        for (int64_t j = 0; j < d2.fPhi.Rows(); j++) {
                            double s1 = 0., s2 = 0.;
                            for (int64_t p = 0; p < d1.fPhi.Cols(); p++) s1 += d1.fW[p] * d1.fPhi(i, p);
                            for (int64_t q = 0; q < d2.fPhi.Cols(); q++) s2 += d2.fW[q] * d2.fPhi(j, q);
                            prod(i, j) = s1 * s2;
                        }
                    difProduto = std::max(difProduto, DifRel(ce.fMat, prod));
                }
            }
        if (caso == 0)
            std::cout << "  2D, 6 quadrilateros de ordem 2: CalcStiffC corrigido x anterior: " << difAnterior
                      << "; simetria C(e1,e2) = C(e2,e1)^T: " << difSim << "\n";
        else
            std::cout << "  2D, triangulo e quadrilatero de ordem 2, G = 1: C_ij x (int N_i)(int N_j): " << difProduto
                      << "; simetria: " << difSim << "\n";
        delete cmesh;
        delete gmesh;
    }
}

struct Resultado {
    Malha m;
    int64_t nel = 0, neq = 0;
    double h = 0.;
    std::vector<double> lambda; // os kNModos maiores, em ordem decrescente
    double soma = 0., traco = 0., menor = 0.;
    double tMalha = 0., tPontos = 0., tB = 0., tC = 0., tAuto = 0., tTotal = 0., memoria = 0.;
};

Resultado Resolve(const Malha &m, int p, int ng, int ngv, int ngp, int nthreads, bool validar) {
    Resultado res;
    res.m = m;
    const auto t0 = Relogio::now();

    // ---- malhas ----
    std::vector<std::array<int, 3>> celula;
    TPZGeoMesh *gmesh = MalhaMeioCilindro(m, celula);
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(p);
    cmesh->SetDimModel(3);
    cmesh->SetAllCreateFunctionsContinuous();
    // tipo 4: G = exp(-dist/Lx), com Lx = b L
    auto *mat = new TPZKarhunenLoeveMat(kMatVolume, kB * kL, kB * kL, kB * kL, 3, 4);
    cmesh->InsertMaterialObject(mat);
    std::set<int> ids = {kMatVolume};
    cmesh->AutoBuild(ids);
    const int64_t neq = cmesh->NEquations();
    std::vector<ElementoKL> el;
    for (auto *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 3) continue;
        ElementoKL e;
        e.intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        e.eq = Equacoes(cmesh, cel);
        e.celula = celula[cel->Reference()->Index()];
        if (!e.intel || (int64_t)e.eq.size() != e.intel->NShapeF()) {
            std::cerr << "elemento inesperado na malha computacional\n";
            std::exit(1);
        }
        el.push_back(std::move(e));
    }
    const int nel = (int)el.size();
    res.nel = nel;
    res.neq = neq;
    // elementos com no em comum e maior diagonal
    {
        std::vector<std::vector<int>> doNo(gmesh->NNodes());
        for (int e = 0; e < nel; e++) {
            TPZGeoEl *gel = el[e].intel->Reference();
            for (int in = 0; in < gel->NCornerNodes(); in++) doNo[gel->NodeIndex(in)].push_back(e);
            for (int a = 0; a < 8; a++)
                for (int b = a + 1; b < 8; b++) {
                    TPZManVector<REAL, 3> xa(3), xb(3);
                    gel->Node(a).GetCoordinates(xa);
                    gel->Node(b).GetCoordinates(xb);
                    res.h = std::max(res.h, std::hypot(xa[0] - xb[0], xa[1] - xb[1], xa[2] - xb[2]));
                }
        }
        for (int e = 0; e < nel; e++) {
            TPZGeoEl *gel = el[e].intel->Reference();
            for (int in = 0; in < gel->NCornerNodes(); in++)
                for (int o : doNo[gel->NodeIndex(in)])
                    if (o != e) el[e].toca.push_back(o);
            std::sort(el[e].toca.begin(), el[e].toca.end());
            el[e].toca.erase(std::unique(el[e].toca.begin(), el[e].toca.end()), el[e].toca.end());
        }
    }
    res.tMalha = Segundos(t0);

    // ---- dados nos pontos de integracao (coordenadas, pesos x |det J| e funcoes de forma), uma vez por elemento ----
    const auto t1 = Relogio::now();
    const int nElem = p + 1; // a regra do elemento (ordem 2p): p + 1 pontos por direcao
    if (ng <= 0) ng = nElem;
    if (ngv <= 0) ngv = ng;
    if (ngp <= 0) ngp = ng;
    const TPZIntCube3D regra = RegraCubo(ng), regraV = RegraCubo(ngv), regraP = RegraCubo(ngp);
    double erroX = 0., erroW = 0., volume = 0.;
    for (auto &e : el) {
        e.intel->KLPointData(regra, e.pontos);
        e.fonteVizinho = ngv == ng ? 0 : 1;
        e.fonteProprio = ngp == ng ? 0 : (ngp == ngv ? 1 : 2);
        if (e.fonteVizinho == 1) e.intel->KLPointData(regraV, e.pontosVizinho);
        if (e.fonteProprio == 2) e.intel->KLPointData(regraP, e.pontosProprio);
        // geometria: pontos e |det J| contra o mapa exato do cilindro
        TPZManVector<REAL, 3> qsi(3);
        REAL w;
        for (int ip = 0; ip < regra.NPoints(); ip++) {
            regra.Point(ip, qsi, w);
            REAL x[3], dj;
            MapaExato(m, e.celula, qsi, x, dj);
            for (int k = 0; k < 3; k++) erroX = std::max(erroX, std::fabs(e.pontos.fX(k, ip) - x[k]));
            erroW = std::max(erroW, std::fabs(e.pontos.fW[ip] - w * dj) / (w * dj));
            volume += e.pontos.fW[ip];
        }
    }
    res.tPontos = Segundos(t1);

    // ---- B: CalcStiffB de cada elemento ----
    const auto t2 = Relogio::now();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(neq, neq);
    for (auto &e : el) {
        TPZElementMatrixT<REAL> be;
        e.intel->CalcStiffB(be);
        Espalha(be.fMat, e.eq, e.eq, false, B);
    }
    res.tB = Segundos(t2);

    std::cout << "\nMalha " << m.nt << " x " << m.nr << " x " << m.nz << " (theta x rho x z), ordem " << p << ": " << nel
              << " elementos, " << neq << " equacoes, maior diagonal " << std::setprecision(4) << res.h << "\n";
    std::cout << std::setprecision(10) << "  volume (soma dos pesos x |det J|) = " << volume << ", 270 pi = " << kVolume
              << ", erro relativo " << std::setprecision(2) << std::scientific << (volume - kVolume) / kVolume
              << std::defaultfloat << "\n";
    std::cout << "  geometria contra o mapa exato do cilindro: max |x - x_exato| = " << std::setprecision(2)
              << std::scientific << erroX << ", max erro relativo de |det J| = " << erroW << std::defaultfloat << "\n";
    std::cout << "  pontos de Gauss por direcao: " << ng << " (pares sem no em comum), " << ngv
              << " (pares com no em comum), " << ngp << " (elemento com ele mesmo)" << std::endl;

    // ---- validacao dos pares: algoritmo anterior, simetria e caso 2D ----
    if (validar) {
        const auto tv = Relogio::now();
        // pares (e1, e2) com e1 entre os 3 primeiros elementos e e2 entre o proprio, um vizinho e o ultimo elemento
        double difSim = 0., difAnterior = 0.;
        for (int e1 = 0; e1 < std::min(nel, 3); e1++) {
            std::vector<int> outros = {e1, nel - 1};
            if (!el[e1].toca.empty()) outros.push_back(el[e1].toca.front());
            for (int e2 : outros) {
                TPZElementMatrixT<REAL> ce, ce21;
                el[e1].intel->CalcStiffC(el[e2].intel, ce);
                el[e2].intel->CalcStiffC(el[e1].intel, ce21);
                difSim = std::max(difSim, DifRelTransposta(ce.fMat, ce21.fMat));
                TPZFNMatrix<729, STATE> ant;
                CalcStiffCAnterior(el[e1].intel, el[e2].intel, mat, ant);
                difAnterior = std::max(difAnterior, DifRel(ce.fMat, ant));
            }
        }
        std::cout << "  validacao dos pares (max |dif| / max |C_par|):\n"
                  << "    CalcStiffC corrigido x algoritmo anterior (elementos iguais): " << difAnterior << "\n"
                  << "    simetria C(e1,e2) = C(e2,e1)^T: " << difSim << "\n";
        ValidaCaso2D();
        std::cout << "  (" << std::setprecision(3) << Segundos(tv) << " s)" << std::endl;
    }

    // ---- C: todos os pares (e1 <= e2), em threads; cada thread soma na sua matriz ----
    const auto t3 = Relogio::now();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(neq, neq);
    {
        nthreads = std::max(1, nthreads);
        std::vector<Eigen::MatrixXd> Cextra(nthreads - 1);
        for (auto &c : Cextra) c = Eigen::MatrixXd::Zero(neq, neq);
        std::atomic<int> proxima{0};
        auto trabalho = [&](Eigen::MatrixXd &Cl) {
            TPZFNMatrix<729, STATE> ce;
            std::vector<char> toca(nel, 0);
            for (int e1; (e1 = proxima++) < nel;) {
                const ElementoKL &E1 = el[e1];
                for (int o : E1.toca) toca[o] = 1;
                for (int e2 = e1; e2 < nel; e2++) {
                    const ElementoKL &E2 = el[e2];
                    if (e2 == e1) mat->ContributeC(E1.Proprio(), E1.Proprio(), ce);
                    else if (toca[e2]) mat->ContributeC(E1.Vizinho(), E2.Vizinho(), ce);
                    else mat->ContributeC(E1.pontos, E2.pontos, ce);
                    Espalha(ce, E1.eq, E2.eq, e2 != e1, Cl);
                }
                for (int o : E1.toca) toca[o] = 0;
            }
        };
        std::vector<std::thread> th;
        for (int t = 1; t < nthreads; t++) th.emplace_back(trabalho, std::ref(Cextra[t - 1]));
        trabalho(C);
        for (auto &t : th) t.join();
        for (auto &c : Cextra) C += c;
    }
    res.tC = Segundos(t3);

    // ---- validacao global: C montada com CalcStiffC de cada par ordenado (e1, e2), sem usar a simetria ----
    if (validar && nel <= 64) {
        const auto tv = Relogio::now();
        Eigen::MatrixXd Cl = Eigen::MatrixXd::Zero(neq, neq);
        for (int e1 = 0; e1 < nel; e1++)
            for (int e2 = 0; e2 < nel; e2++) {
                const bool vizinho = std::binary_search(el[e1].toca.begin(), el[e1].toca.end(), e2);
                const TPZIntCube3D &r = e1 == e2 ? regraP : (vizinho ? regraV : regra);
                TPZElementMatrixT<REAL> ce;
                if (ng == nElem && !vizinho && e1 != e2) el[e1].intel->CalcStiffC(el[e2].intel, ce); // regras dos elementos
                else el[e1].intel->CalcStiffC(el[e2].intel, r, r, ce);
                Espalha(ce.fMat, el[e1].eq, el[e2].eq, false, Cl);
            }
        std::cout << "  validacao global: max |C - C_CalcStiffC| / max |C| = " << std::setprecision(3) << std::scientific
                  << (C - Cl).cwiseAbs().maxCoeff() / C.cwiseAbs().maxCoeff() << std::defaultfloat << " ("
                  << std::setprecision(3) << Segundos(tv) << " s)" << std::endl;
    }

    // ---- autovalores de C v = lambda B v: B = L L^T, A = L^-1 C L^-T (no lugar de B e C) ----
    const auto t4 = Relogio::now();
    {
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(B);
        if (llt.info() != Eigen::Success) {
            std::cerr << "B nao e positiva definida\n";
            std::exit(1);
        }
        llt.matrixL().solveInPlace(C);
        llt.matrixU().solveInPlace<Eigen::OnTheRight>(C);
    }
    res.traco = C.trace();
    Eigen::VectorXd lam;
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C, Eigen::EigenvaluesOnly);
        if (es.info() != Eigen::Success) {
            std::cerr << "o calculo dos autovalores nao convergiu\n";
            std::exit(1);
        }
        lam = es.eigenvalues(); // ordem crescente
    }
    res.tAuto = Segundos(t4);
    res.soma = lam.sum();
    res.menor = lam(0);
    for (int k = 0; k < kNModos && k < lam.size(); k++) res.lambda.push_back(lam(lam.size() - 1 - k));
    res.tTotal = Segundos(t0);
    res.memoria = PicoMemoriaMB();

    std::cout << "  modo    NeoPZ            Galerkin 2       erro rel.     Galerkin 1       (G1-G2)/G2\n";
    for (int k = 0; k < (int)res.lambda.size(); k++)
        std::cout << "  " << std::setw(4) << k + 1 << std::fixed << std::setprecision(7) << std::setw(17) << res.lambda[k]
                  << std::setw(17) << kGalerkin2[k] << std::scientific << std::setprecision(2) << std::setw(13)
                  << (res.lambda[k] - kGalerkin2[k]) / kGalerkin2[k] << std::fixed << std::setprecision(7)
                  << std::setw(17) << kGalerkin1[k] << std::scientific << std::setprecision(2) << std::setw(13)
                  << (kGalerkin1[k] - kGalerkin2[k]) / kGalerkin2[k] << std::defaultfloat << "\n";
    std::cout << std::setprecision(10) << "  soma de todos os autovalores = " << res.soma << " (traco de L^-1 C L^-T = "
              << res.traco << "); sigma^2 |D| = 270 pi = " << kVolume << "; razao " << std::setprecision(6)
              << res.soma / kVolume << "\n";
    std::cout << "  menor autovalor = " << std::setprecision(3) << std::scientific << res.menor << std::defaultfloat
              << "\n";
    std::cout << std::fixed << std::setprecision(2) << "  tempos (s): malhas " << res.tMalha << ", pontos " << res.tPontos
              << ", B " << res.tB << ", C " << res.tC << " (" << nthreads << " threads), autovalores " << res.tAuto
              << ", total " << res.tTotal << "; pico de memoria " << std::setprecision(0) << res.memoria << " MB"
              << std::defaultfloat << std::endl;

    // a malha computacional usa o material e a geometrica
    delete cmesh;
    delete gmesh;
    return res;
}

} // namespace

int main() {
    const int p = envi("P", 2);
    if (p < 1 || p > 3) {
        std::cerr << "P deve ser 1, 2 ou 3\n";
        return 2;
    }
    const std::string padrao = p == 1 ? "16x2x8,24x3x12,32x4x16" : "8x1x4,12x2x6,16x2x8";
    const std::vector<Malha> malhas = LeMalhas(envs("MALHAS", padrao));
    const int ng = envi("NG", 0), ngv = envi("NG_VIZINHO", 0), ngp = envi("NG_PROPRIO", 0);
    const int nthreads = envi("THREADS", 2);
    const bool validar = envi("VALIDAR", 0) != 0;
    const std::string csv = envs("CSV", "");

    std::cout << "Karhunen-Loeve de Galerkin no meio cilindro vazado de Mika et al. (2021), Secao 6.2, Exemplo 1-1\n"
              << "  8 <= rho <= 10, 0 <= theta <= pi, 0 <= z <= 15; G = exp(-|x - x'| / 5), sigma^2 = 1\n"
              << "  referencia: Galerkin 2 da Tabela 4 (NURBS p = 2, 3800 elementos, 6642 graus de liberdade)\n";

    std::vector<Resultado> res;
    for (const auto &m : malhas) {
        res.push_back(Resolve(m, p, ng, ngv, ngp, nthreads, validar));
        if (!csv.empty()) {
            const bool novo = !std::ifstream(csv).good();
            std::ofstream f(csv, std::ios::app);
            if (novo) {
                f << "p,nt,nr,nz,elementos,equacoes,h,ng,ng_vizinho,ng_proprio";
                for (int k = 1; k <= kNModos; k++) f << ",lambda" << k;
                f << ",soma,menor,t_pontos,t_B,t_C,t_autovalores,t_total,threads,memoria_MB\n";
            }
            const Resultado &r = res.back();
            const int ngE = ng > 0 ? ng : p + 1;
            f << p << ',' << m.nt << ',' << m.nr << ',' << m.nz << ',' << r.nel << ',' << r.neq << ','
              << std::setprecision(6) << r.h << ',' << ngE << ',' << (ngv > 0 ? ngv : ngE) << ','
              << (ngp > 0 ? ngp : ngE) << std::setprecision(12);
            for (double l : r.lambda) f << ',' << l;
            f << ',' << r.soma << ',' << r.menor << std::setprecision(6) << ',' << r.tPontos << ',' << r.tB << ','
              << r.tC << ',' << r.tAuto << ',' << r.tTotal << ',' << nthreads << ',' << r.memoria << '\n';
        }
    }

    // ---- resumo da convergencia ----
    std::cout << "\nResumo (erros relativos contra Galerkin 2; soma dos autovalores / 270 pi)\n"
              << "  malha        elementos  equacoes  lambda_1       lambda_10      lambda_20      max |erro| 1-20  "
                 "soma/270pi  tempo (s)\n";
    for (const auto &r : res) {
        double emax = 0.;
        for (int k = 0; k < (int)r.lambda.size(); k++)
            emax = std::max(emax, std::fabs(r.lambda[k] - kGalerkin2[k]) / kGalerkin2[k]);
        std::ostringstream nome;
        nome << r.m.nt << 'x' << r.m.nr << 'x' << r.m.nz;
        auto erro = [&](int k) { return (r.lambda[k] - kGalerkin2[k]) / kGalerkin2[k]; };
        std::cout << "  " << std::left << std::setw(12) << nome.str() << std::right << std::setw(10) << r.nel
                  << std::setw(10) << r.neq << std::scientific << std::setprecision(3) << std::setw(15) << erro(0)
                  << std::setw(15) << erro(9) << std::setw(15) << erro(19) << std::setw(17) << emax << std::fixed
                  << std::setprecision(5) << std::setw(12) << r.soma / kVolume << std::setprecision(1) << std::setw(11)
                  << r.tTotal << std::defaultfloat << "\n";
    }
    return 0;
}
