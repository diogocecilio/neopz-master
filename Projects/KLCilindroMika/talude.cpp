// KLTaludeCeron - expansao de Karhunen-Loeve de Galerkin do NeoPZ com o nucleo exponencial separavel de Vargas Ceron,
// Cecilio, Linn e Maghous (2025, eq. 64), G(x, x') = exp(-|x - x'| / Lx - |y - y'| / Ly), Lx = 20 m, Ly = 2 m
// (TPZKarhunenLoeveMat do tipo 3, CalcStiffB e ContributeC nos pares de elementos, o mesmo esquema de main.cpp), comparada
// com a expansao em forma fechada no retangulo que o SlopeRandomField usa (../SlopeRandomField/KLSeparavel.h,
// klsep::KL2D: os autopares 2D sao produtos dos autopares 1D do nucleo exponencial).
//   (a) retangulo [0,25] x [0,10]: autovalores numericos contra os da forma fechada (produtos dos autovalores 1D), em
//       malhas estruturadas de triangulos (cada celula nx x ny dividida em dois triangulos, diagonais alternadas).
//   (b) talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10), area 187,5 m2, malhas do Gmsh de
//       LearningMohrCoulomb2/cluster/malhas: a KL numerica no proprio talude contra a KL do retangulo [0,25] x [0,10]
//       (a caixa dos nos da malha) avaliada no talude, com o mesmo numero de termos M: erro medio quadratico de
//       truncamento 1 - media de sigma_M^2 no talude (sigma_M^2(x) = sum_{k<=M} lambda_k phi_k(x)^2), minimo e maximo
//       de sigma_M^2 em pontos de cada elemento, e a covariancia truncada sum_{k<=M} lambda_k phi_k(x) phi_k(x') contra
//       o nucleo exato em pares de pontos aleatorios no talude.
// Parametros (variaveis de ambiente):
//   PARTE=ab              a, b ou ab
//   LX=20 LY=2            comprimentos de correlacao (m)
//   RETANGULOS=25x10:1,50x20:1,100x40:1,25x10:2,50x20:2   malhas do retangulo, nx x ny celulas : ordem
//   TALUDES=ceron_h5_a45.msh:1,ceron_h5_a45_h025.msh:1,ceron_h5_a45_fina.msh:1,ceron_h5_a45.msh:2,
//           ceron_h5_a45_h025.msh:2                       malhas do talude (Gmsh 4.1) : ordem
//   MALHAS_DIR=<fonte>/../LearningMohrCoulomb2/cluster/malhas   diretorio dos arquivos de TALUDES sem caminho
//   M=50,200,1000,2000    numeros de termos comparados no talude (os maiores que o numero de equacoes ficam de fora)
//   ORDEM=0               ordem da regra de integracao nos pares de elementos sem no em comum (0: a do elemento, 2P)
//   ORDEM_VIZINHO=0       idem, nos pares de elementos distintos com no em comum (0: a de ORDEM)
//   ORDEM_PROPRIO=0       idem, no par de um elemento com ele mesmo (0: a de ORDEM_VIZINHO)
//   PARES=500             pares de pontos aleatorios no talude (uniformes e proximos)
//   SEMENTE=1             semente dos pontos aleatorios
//   REDE=4                pontos por aresta - 1 na rede de pontos de cada elemento onde sigma_M^2 e avaliada
//   RET_MAX=20000         termos da forma fechada usados para achar o M equivalente do retangulo
//   THREADS=2             threads na montagem de C
//   CSV_RET=, CSV_TALUDE= arquivos CSV onde acrescentar os resultados de (a) e de (b)
//   AUTOVALORES=          prefixo: grava todos os autovalores numericos de cada malha, decrescentes, um por linha, em
//                         <prefixo>_ret_<nx>x<ny>_p<P>_o<ordem>.txt (a) e <prefixo>_<malha>_p<P>_o<ordem>.txt (b)
//   PONTOS=               prefixo: grava, para cada malha do talude, <prefixo>_<malha>_p<P>_o<ordem>.csv com os pontos da
//                         rede de cada elemento (x, y), sigma_M^2 da KL numerica e da forma fechada para cada M e os
//                         PONTOS_MODOS (9) primeiros modos numericos phi_k (B-ortonormais, sinal arbitrario)
// Ver README.md.
#include <Eigen/Dense>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeotriangle.h"
#include "pzinterpolationspace.h"
#include "TPZElementMatrixT.h"
#include "TPZMaterialDataT.h"
#include "TPZKarhunenLoeveMat.h"
#include "TPZGmshReader.h"
#include "pzquad.h"
#include "tpzintrulet.h"

#include "../SlopeRandomField/KLSeparavel.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifndef MALHAS_DIR
#define MALHAS_DIR "."
#endif

namespace {

constexpr int kMatSolo = 1;
// retangulo do SlopeRandomField (a caixa dos nos da malha do talude) e o talude de Vargas Ceron et al. (2025), Secao 6:
// base em y = 0, pe em (15, 5), crista em (10, 10)
constexpr double kX0 = 0., kY0 = 0., kDx = 25., kDy = 10.;
constexpr double kAreaRetangulo = kDx * kDy, kAreaTalude = 187.5;

// ponto no talude: 0 <= x <= 25, 0 <= y <= 10 e, acima do terreno do pe (y > 5), a esquerda da face x + y = 20
bool NoTalude(double x, double y) {
    if (x < kX0 || x > kX0 + kDx || y < kY0 || y > kY0 + kDy) return false;
    return y <= 5. || x + y <= 20.;
}

int envi(const char *k, int d) { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
double envd(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }

using Relogio = std::chrono::steady_clock;
double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }

// pico da memoria residente do processo (VmHWM), em MB, desde o ultimo ZeraPicoMemoria
double PicoMemoriaMB() {
    std::ifstream f("/proc/self/status");
    std::string l;
    while (std::getline(f, l))
        if (l.rfind("VmHWM:", 0) == 0) return std::atof(l.c_str() + 6) / 1024.;
    return -1.;
}

// recomeca o pico de memoria (Linux: "5" em /proc/self/clear_refs), para medir cada malha separadamente
void ZeraPicoMemoria() { std::ofstream("/proc/self/clear_refs") << "5"; }

std::vector<std::string> Separa(const std::string &s, char c) {
    std::vector<std::string> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, c))
        if (!t.empty()) v.push_back(t);
    return v;
}

struct MalhaRet { int nx = 0, ny = 0, p = 1; };
struct MalhaTalude { std::string arquivo; int p = 1; };

// "25x10:1,50x20:2"
std::vector<MalhaRet> LeRetangulos(const std::string &s) {
    std::vector<MalhaRet> v;
    for (const auto &t : Separa(s, ',')) {
        MalhaRet m;
        if (std::sscanf(t.c_str(), "%dx%d:%d", &m.nx, &m.ny, &m.p) != 3 || m.nx < 1 || m.ny < 1 || m.p < 1 || m.p > 3) {
            std::cerr << "malha invalida em RETANGULOS: " << t << " (use nx x ny : ordem, por exemplo 50x20:1)\n";
            std::exit(2);
        }
        v.push_back(m);
    }
    return v;
}

// "ceron_h5_a45.msh:1,ceron_h5_a45.msh:2"; arquivo sem '/' fica relativo a MALHAS_DIR
std::vector<MalhaTalude> LeTaludes(const std::string &s, const std::string &dir) {
    std::vector<MalhaTalude> v;
    for (const auto &t : Separa(s, ',')) {
        const auto pos = t.rfind(':');
        MalhaTalude m;
        m.arquivo = pos == std::string::npos ? t : t.substr(0, pos);
        m.p = pos == std::string::npos ? 1 : std::atoi(t.c_str() + pos + 1);
        if (m.arquivo.find('/') == std::string::npos) m.arquivo = dir + "/" + m.arquivo;
        if (m.p < 1 || m.p > 3 || !std::ifstream(m.arquivo).good()) {
            std::cerr << "malha invalida em TALUDES: " << t << " (arquivo " << m.arquivo << ", ordem 1 a 3)\n";
            std::exit(2);
        }
        v.push_back(m);
    }
    return v;
}

std::string Nome(const std::string &caminho) { return caminho.substr(caminho.find_last_of('/') + 1); }

// n pontos de Gauss-Legendre em [-1, 1] (Newton nas raizes de P_n)
void GaussLegendre(int n, std::vector<double> &x, std::vector<double> &w) {
    x.assign(n, 0.);
    w.assign(n, 0.);
    for (int i = 0; i < (n + 1) / 2; i++) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 0.;
        for (int it = 0; it < 100; it++) {
            double p0 = 1., p1 = 0.;
            for (int k = 1; k <= n; k++) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2. * k - 1.) * z * p1 - (k - 1.) * p2) / k;
            }
            dp = n * (z * p0 - p1) / (z * z - 1.);
            const double dz = p0 / dp;
            z -= dz;
            if (std::fabs(dz) < 1.e-15) break;
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = w[n - 1 - i] = 2. / ((1. - z * z) * dp * dp);
    }
}

// ---- malhas ----

// retangulo [0,25] x [0,10] com nx x ny celulas, cada uma dividida em dois triangulos (diagonais alternadas)
TPZGeoMesh *MalhaRetangulo(int nx, int ny) {
    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    std::vector<int64_t> no((size_t)(nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) {
            const int64_t id = gmesh->NodeVec().AllocateNewElement();
            TPZManVector<REAL, 3> co = {kX0 + kDx * i / nx, kY0 + kDy * j / ny, 0.};
            gmesh->NodeVec()[id].Initialize(co, *gmesh);
            no[(size_t)i + (nx + 1) * (size_t)j] = id;
        }
    using Tri = TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>;
    auto n = [&](int i, int j) { return no[(size_t)i + (nx + 1) * (size_t)j]; };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const int64_t a = n(i, j), b = n(i + 1, j), c = n(i + 1, j + 1), d = n(i, j + 1);
            const bool ac = (i + j) % 2 == 0; // diagonal a-c ou b-d
            TPZManVector<int64_t, 3> t1 = {a, b, ac ? c : d}, t2 = {ac ? a : b, c, d};
            new Tri(t1, kMatSolo, *gmesh);
            new Tri(t2, kMatSolo, *gmesh);
        }
    gmesh->BuildConnectivity();
    return gmesh;
}

// malha do talude do Gmsh (.msh 4.1), com os grupos fisicos de LearningMohrCoulomb2/cluster/malhas; so o grupo "Solo"
// (dimensao 2) recebe elementos computacionais
TPZGeoMesh *MalhaTaludeGmsh(const std::string &arquivo) {
    TPZManVector<std::map<std::string, int>, 4> nomes(4);
    nomes[2]["Solo"] = kMatSolo;
    nomes[1]["Base"] = 2;
    nomes[1]["Direita"] = 3;
    nomes[1]["Superficie"] = 4;
    nomes[1]["Esquerda"] = 5;
    TPZGmshReader leitor;
    leitor.SetDimNamePhysical(nomes);
    TPZGeoMesh *gmesh = leitor.GeometricGmshMesh(arquivo);
    gmesh->BuildConnectivity();
    return gmesh;
}

// ---- KL de Galerkin ----

struct ElementoKL {
    TPZInterpolationSpace *intel = nullptr;
    std::vector<int64_t> eq;  // equacao global de cada funcao de forma (na ordem das linhas de phi)
    std::vector<int> toca;    // elementos distintos com algum no em comum
    double caixa[4] = {0., 0., 0., 0.}; // xmin, xmax, ymin, ymax
    // dados nos pontos das regras dos pares sem no em comum, com no em comum e do elemento com ele mesmo (pontosVizinho
    // e pontosProprio ficam vazios quando a regra e a mesma de pontos)
    TPZKLPointData pontos, pontosVizinho, pontosProprio;
    int fonteVizinho = 0, fonteProprio = 0; // 0: pontos, 1: pontosVizinho, 2: pontosProprio
    const TPZKLPointData &Fonte(int f) const { return f == 0 ? pontos : (f == 1 ? pontosVizinho : pontosProprio); }
    const TPZKLPointData &Vizinho() const { return Fonte(fonteVizinho); }
    const TPZKLPointData &Proprio() const { return Fonte(fonteProprio); }
};

struct KLNumerica {
    TPZGeoMesh *gmesh = nullptr;
    TPZCompMesh *cmesh = nullptr;
    TPZKarhunenLoeveMat *mat = nullptr;
    std::vector<ElementoKL> el;
    int p = 1;
    int64_t neq = 0;
    int ordem = 0, ordemV = 0, ordemP = 0;       // ordens das regras usadas
    int npontos = 0, npontosV = 0, npontosP = 0; // pontos por elemento em cada regra (do primeiro elemento)
    double area = 0., hmax = 0.;
    Eigen::VectorXd lambda; // todos os autovalores, em ordem decrescente
    // os primeiros nv autovetores generalizados (B-ortonormais: int phi_k^2 = 1), transpostos: Vt(k, i) e o coeficiente
    // da equacao i no modo k (os nv coeficientes de uma equacao ficam contiguos)
    Eigen::MatrixXd Vt;
    int64_t NModos() const { return Vt.rows(); }
    double tMalha = 0., tPontos = 0., tB = 0., tC = 0., tAuto = 0.;
    KLNumerica() = default;
    KLNumerica(const KLNumerica &) = delete;
    KLNumerica &operator=(const KLNumerica &) = delete;
    ~KLNumerica() {
        delete cmesh; // a malha computacional apaga o material
        delete gmesh;
    }
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

// M(eq1, eq2) += ce e, se transposta, M(eq2, eq1) += ce^T
void Espalha(const TPZFMatrix<STATE> &ce, const std::vector<int64_t> &eq1, const std::vector<int64_t> &eq2,
             bool transposta, Eigen::MatrixXd &M) {
    const int64_t n1 = (int64_t)eq1.size(), n2 = (int64_t)eq2.size();
    if (n1 == 0 || n2 == 0) return;
    if (ce.Rows() != n1 || ce.Cols() != n2) {
        std::cerr << "matriz do par com dimensoes diferentes das equacoes dos elementos\n";
        std::exit(1);
    }
    const STATE *v = &ce.g(0, 0); // por colunas
    for (int64_t j = 0; j < n2; j++) {
        double *col = M.data() + eq2[j] * M.rows();
        for (int64_t i = 0; i < n1; i++) col[eq1[i]] += v[i + j * n1];
    }
    if (!transposta) return;
    for (int64_t i = 0; i < n1; i++) {
        double *col = M.data() + eq1[i] * M.rows();
        for (int64_t j = 0; j < n2; j++) col[eq2[j]] += v[i + j * n1];
    }
}

// monta B e C, resolve C v = lambda B v (Cholesky de B e autovalores de L^-1 C L^-T) e guarda os nvetores primeiros
// autovetores (nvetores = 0: so os autovalores)
void ResolveKL(KLNumerica &kl, TPZGeoMesh *gmesh, int p, double Lx, double Ly, int ordem, int ordemV, int ordemP,
               int nthreads, int64_t nvetores) {
    const auto t0 = Relogio::now();
    kl.gmesh = gmesh;
    kl.p = p;
    auto *cmesh = new TPZCompMesh(gmesh);
    kl.cmesh = cmesh;
    cmesh->SetDefaultOrder(p);
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsContinuous();
    // tipo 3: G = exp(-|x1 - x2| / Lx - |y1 - y2| / Ly) (Vargas Ceron et al. 2025, eq. 64)
    kl.mat = new TPZKarhunenLoeveMat(kMatSolo, Lx, Ly, 1., 2, 3);
    cmesh->InsertMaterialObject(kl.mat);
    std::set<int> ids = {kMatSolo};
    cmesh->AutoBuild(ids);
    kl.neq = cmesh->NEquations();
    for (auto *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 2) continue;
        ElementoKL e;
        e.intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        e.eq = Equacoes(cmesh, cel);
        if (!e.intel || (int64_t)e.eq.size() != e.intel->NShapeF()) {
            std::cerr << "elemento inesperado na malha computacional\n";
            std::exit(1);
        }
        TPZGeoEl *gel = cel->Reference();
        e.caixa[0] = e.caixa[2] = 1.e300;
        e.caixa[1] = e.caixa[3] = -1.e300;
        for (int a = 0; a < gel->NCornerNodes(); a++) {
            TPZManVector<REAL, 3> xa(3);
            gel->Node(a).GetCoordinates(xa);
            e.caixa[0] = std::min(e.caixa[0], (double)xa[0]);
            e.caixa[1] = std::max(e.caixa[1], (double)xa[0]);
            e.caixa[2] = std::min(e.caixa[2], (double)xa[1]);
            e.caixa[3] = std::max(e.caixa[3], (double)xa[1]);
            for (int b = a + 1; b < gel->NCornerNodes(); b++) {
                TPZManVector<REAL, 3> xb(3);
                gel->Node(b).GetCoordinates(xb);
                kl.hmax = std::max(kl.hmax, std::hypot(xa[0] - xb[0], xa[1] - xb[1]));
            }
        }
        kl.el.push_back(std::move(e));
    }
    const int nel = (int)kl.el.size();
    {
        std::vector<std::vector<int>> doNo(gmesh->NNodes());
        for (int e = 0; e < nel; e++) {
            TPZGeoEl *gel = kl.el[e].intel->Reference();
            for (int in = 0; in < gel->NCornerNodes(); in++) doNo[gel->NodeIndex(in)].push_back(e);
        }
        for (int e = 0; e < nel; e++) {
            TPZGeoEl *gel = kl.el[e].intel->Reference();
            auto &t = kl.el[e].toca;
            for (int in = 0; in < gel->NCornerNodes(); in++)
                for (int o : doNo[gel->NodeIndex(in)])
                    if (o != e) t.push_back(o);
            std::sort(t.begin(), t.end());
            t.erase(std::unique(t.begin(), t.end()), t.end());
        }
    }
    kl.tMalha = Segundos(t0);

    // ---- dados nos pontos de integracao (coordenadas, pesos x |det J| e funcoes de forma), uma vez por elemento ----
    const auto t1 = Relogio::now();
    if (ordem <= 0) ordem = 2 * p; // a regra do elemento
    if (ordemV <= 0) ordemV = ordem;
    if (ordemP <= 0) ordemP = ordemV;
    for (int o : {ordem, ordemV, ordemP})
        if (o > TPZIntRuleT::NRULESTRIANGLE_ORDER) {
            std::cerr << "ordem de integracao " << o << " acima da maxima do triangulo ("
                      << TPZIntRuleT::NRULESTRIANGLE_ORDER << ")\n";
            std::exit(2);
        }
    kl.ordem = ordem;
    kl.ordemV = ordemV;
    kl.ordemP = ordemP;
    std::map<std::pair<int, int>, std::unique_ptr<TPZIntPoints>> regras; // (tipo do elemento, ordem) -> regra
    auto regra = [&](TPZGeoEl *gel, int o) -> const TPZIntPoints & {
        const auto chave = std::make_pair((int)gel->Type(), o);
        auto it = regras.find(chave);
        if (it == regras.end())
            it = regras.emplace(chave, std::unique_ptr<TPZIntPoints>(gel->CreateSideIntegrationRule(gel->NSides() - 1, o)))
                     .first;
        return *it->second;
    };
    kl.area = 0.;
    for (auto &e : kl.el) {
        TPZGeoEl *gel = e.intel->Reference();
        e.intel->KLPointData(regra(gel, ordem), e.pontos);
        e.fonteVizinho = ordemV == ordem ? 0 : 1;
        e.fonteProprio = ordemP == ordem ? 0 : (ordemP == ordemV ? 1 : 2);
        if (e.fonteVizinho == 1) e.intel->KLPointData(regra(gel, ordemV), e.pontosVizinho);
        if (e.fonteProprio == 2) e.intel->KLPointData(regra(gel, ordemP), e.pontosProprio);
        for (int ip = 0; ip < e.pontos.fW.size(); ip++) kl.area += e.pontos.fW[ip];
    }
    if (nel > 0) {
        kl.npontos = (int)kl.el[0].pontos.fW.size();
        kl.npontosV = (int)kl.el[0].Vizinho().fW.size();
        kl.npontosP = (int)kl.el[0].Proprio().fW.size();
    }
    kl.tPontos = Segundos(t1);

    // ---- B: CalcStiffB de cada elemento ----
    const auto t2 = Relogio::now();
    const int64_t neq = kl.neq;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(neq, neq);
    for (auto &e : kl.el) {
        TPZElementMatrixT<REAL> be;
        e.intel->CalcStiffB(be);
        Espalha(be.fMat, e.eq, e.eq, false, B);
    }
    kl.tB = Segundos(t2);

    // ---- C: todos os pares (e1 <= e2), em threads; cada thread soma na sua matriz ----
    const auto t3 = Relogio::now();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(neq, neq);
    {
        nthreads = std::max(1, nthreads);
        std::vector<Eigen::MatrixXd> Cextra(nthreads - 1);
        for (auto &c : Cextra) c = Eigen::MatrixXd::Zero(neq, neq);
        std::atomic<int> proxima{0};
        const TPZKarhunenLoeveMat *mat = kl.mat;
        const auto &el = kl.el;
        auto trabalho = [&](Eigen::MatrixXd &Cl) {
            TPZFNMatrix<400, STATE> ce;
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
    kl.tC = Segundos(t3);

    // ---- autopares de C v = lambda B v: B = L L^T, A = L^-1 C L^-T (no lugar de B e C), v = L^-T y ----
    const auto t4 = Relogio::now();
    Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(B);
    if (llt.info() != Eigen::Success) {
        std::cerr << "B nao e positiva definida\n";
        std::exit(1);
    }
    llt.matrixL().solveInPlace(C);
    llt.matrixU().solveInPlace<Eigen::OnTheRight>(C);
    const int64_t nv = std::min(nvetores, neq);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(C, nv > 0 ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
    C.resize(0, 0);
    if (es.info() != Eigen::Success) {
        std::cerr << "o calculo dos autovalores nao convergiu\n";
        std::exit(1);
    }
    kl.lambda = es.eigenvalues().reverse(); // decrescente
    if (nv > 0) {
        Eigen::MatrixXd V = es.eigenvectors().rightCols(nv).rowwise().reverse();
        llt.matrixU().solveInPlace(V);
        kl.Vt = V.transpose();
    }
    kl.tAuto = Segundos(t4);
}

// dados do elemento no ponto qsi (funcoes de forma em data.phi e coordenadas em data.x)
void DadosNoPonto(TPZInterpolationSpace *intel, TPZVec<REAL> &qsi, TPZMaterialDataT<REAL> &data) {
    intel->ComputeRequiredData(data, qsi);
}

void PreparaDados(TPZInterpolationSpace *intel, TPZMaterialDataT<REAL> &data) {
    intel->InitMaterialData(data);
    data.fNeedsSol = false;
    data.p = intel->MaxOrder();
}

// phi_k(x) = sum_s v_k(eq_s) N_s(x), k < nv, a partir das funcoes de forma do elemento no ponto
void ModosNoPonto(const KLNumerica &kl, const ElementoKL &e, const TPZMaterialDataT<REAL> &data, std::vector<double> &phi) {
    const int64_t nv = kl.NModos(), ns = (int64_t)e.eq.size();
    phi.assign(nv, 0.);
    for (int64_t s = 0; s < ns; s++) {
        const double N = data.phi(s, 0);
        const double *v = kl.Vt.data() + e.eq[s] * nv;
        for (int64_t k = 0; k < nv; k++) phi[k] += N * v[k];
    }
}

// elemento que contem (x, y) e o ponto qsi correspondente
bool Localiza(const KLNumerica &kl, double x, double y, int &elemento, TPZManVector<REAL, 3> &qsi) {
    const double tol = 1.e-9;
    TPZManVector<REAL, 3> X = {x, y, 0.};
    for (int e = 0; e < (int)kl.el.size(); e++) {
        const auto &c = kl.el[e].caixa;
        if (x < c[0] - tol || x > c[1] + tol || y < c[2] - tol || y > c[3] + tol) continue;
        TPZGeoEl *gel = kl.el[e].intel->Reference();
        qsi.Resize(2);
        if (gel->Type() == ETriangle) qsi[0] = qsi[1] = 1. / 3.;
        else qsi[0] = qsi[1] = 0.;
        gel->ComputeXInverse(X, qsi, 1.e-12);
        if (gel->IsInParametricDomain(qsi, 1.e-10)) {
            elemento = e;
            return true;
        }
    }
    return false;
}

// pontos de referencia da rede de um elemento: triangulo (a/n, b/n), a + b <= n; quadrilatero (-1 + 2a/n, -1 + 2b/n)
std::vector<std::array<double, 2>> RedeReferencia(MElementType tipo, int n) {
    std::vector<std::array<double, 2>> r;
    for (int b = 0; b <= n; b++)
        for (int a = 0; a <= n; a++) {
            if (tipo == ETriangle) {
                if (a + b <= n) r.push_back({(double)a / n, (double)b / n});
            } else {
                r.push_back({-1. + 2. * a / n, -1. + 2. * b / n});
            }
        }
    return r;
}

// ---- forma fechada (klsep::KL2D) no talude ----

// integral de f_i(s)^2 para s de s0 a s1 (s medido a partir do centro do intervalo), f = n cos(w s) ou n sin(w s)
double IntegralQuadrado(const klsep::Modos1D &m, int i, double s0, double s1) {
    const double w = m.omega[i], n2 = m.norma[i] * m.norma[i], sinal = m.par[i] ? 1. : -1.;
    return 0.5 * n2 * ((s1 - s0) + sinal * (std::sin(2. * w * s1) - std::sin(2. * w * s0)) / (2. * w));
}

// integral no talude de phi_k^2 = f_i(x)^2 g_j(y)^2 (modo (i, j) do retangulo): na faixa 0 <= y <= 5 o produto das
// integrais 1D; em 5 <= y <= 10, x de 0 a 20 - y: a integral em x em forma fechada e Gauss-Legendre composto em y, com um
// subintervalo de 10 pontos por periodo de cos(2 (w_x + w_y) y)
double IntegralNoTalude(const klsep::KL2D &kl, int i, int j, const std::vector<double> &gx,
                        const std::vector<double> &gw) {
    const double cx = kl.x0 + 0.5 * kl.Dx, cy = kl.y0 + 0.5 * kl.Dy;
    double I = IntegralQuadrado(kl.mx, i, kl.x0 - cx, kl.x0 + kl.Dx - cx) * IntegralQuadrado(kl.my, j, 0. - cy, 5. - cy);
    const double wmax = kl.mx.omega[i] + kl.my.omega[j];
    const int nsub = 2 + (int)std::ceil(5. * wmax / M_PI);
    const double h = 5. / nsub;
    double I2 = 0.;
    for (int s = 0; s < nsub; s++) {
        const double ya = 5. + s * h;
        for (size_t q = 0; q < gx.size(); q++) {
            const double y = ya + 0.5 * h * (gx[q] + 1.), g = kl.my.Avalia(j, y - cy);
            I2 += 0.5 * h * gw[q] * g * g * IntegralQuadrado(kl.mx, i, kl.x0 - cx, 20. - y - cx);
        }
    }
    return I + I2;
}

double Quadrado(double v) { return v * v; }

// nome do arquivo sem o caminho e sem a extensao .msh
std::string Base(const std::string &caminho) {
    std::string b = Nome(caminho);
    if (b.size() > 4 && b.compare(b.size() - 4, 4, ".msh") == 0) b.resize(b.size() - 4);
    return b;
}

// AUTOVALORES: todos os autovalores numericos, decrescentes, um por linha
void GravaAutovalores(const std::string &arquivo, const Eigen::VectorXd &lambda) {
    std::ofstream f(arquivo);
    f << std::setprecision(17);
    for (int64_t k = 0; k < lambda.size(); k++) f << lambda(k) << '\n';
}

// ---- (a) retangulo ----

struct ResultadoRet {
    MalhaRet m;
    int64_t nel = 0, neq = 0;
    double hmax = 0., area = 0.;
    std::vector<double> lambda; // ate 1000 autovalores numericos, decrescentes
    double traco = 0.;
    double tTotal = 0., tC = 0., tAuto = 0., memoria = 0.;
    int ordem = 0, ordemV = 0, ordemP = 0;
};

ResultadoRet ParteA(const MalhaRet &m, double Lx, double Ly, const std::vector<double> &lamRef, int ordem, int ordemV,
                    int ordemP, int nthreads) {
    ZeraPicoMemoria();
    const auto t0 = Relogio::now();
    ResultadoRet r;
    r.m = m;
    KLNumerica kl;
    ResolveKL(kl, MalhaRetangulo(m.nx, m.ny), m.p, Lx, Ly, ordem, ordemV, ordemP, nthreads, 0);
    r.nel = (int64_t)kl.el.size();
    r.neq = kl.neq;
    r.hmax = kl.hmax;
    r.area = kl.area;
    r.ordem = kl.ordem;
    r.ordemV = kl.ordemV;
    r.ordemP = kl.ordemP;
    r.traco = kl.lambda.sum();
    for (int k = 0; k < kl.lambda.size() && k < (int)lamRef.size(); k++) r.lambda.push_back(kl.lambda(k));
    r.tC = kl.tC;
    r.tAuto = kl.tAuto;
    r.tTotal = Segundos(t0);
    r.memoria = PicoMemoriaMB();
    const std::string prefixoAuto = envs("AUTOVALORES", "");
    if (!prefixoAuto.empty())
        GravaAutovalores(prefixoAuto + "_ret_" + std::to_string(m.nx) + "x" + std::to_string(m.ny) + "_p" +
                             std::to_string(m.p) + "_o" + std::to_string(kl.ordem) + ".txt",
                         kl.lambda);

    std::cout << "\nRetangulo " << m.nx << " x " << m.ny << " celulas (" << r.nel << " triangulos), ordem " << m.p << ": "
              << r.neq << " equacoes, maior aresta " << std::setprecision(4) << r.hmax << " m; area (soma dos pesos) "
              << std::setprecision(10) << r.area << "\n";
    std::cout << "  regras (ordem/pontos): pares sem no em comum " << kl.ordem << "/" << kl.npontos
              << ", com no em comum " << kl.ordemV << "/" << kl.npontosV << ", elemento com ele mesmo " << kl.ordemP << "/"
              << kl.npontosP << "\n";
    std::cout << "  modo   numerico         forma fechada    erro rel.\n";
    for (int k = 0; k < 20 && k < (int)r.lambda.size(); k++)
        std::cout << "  " << std::setw(4) << k + 1 << std::fixed << std::setprecision(8) << std::setw(17) << r.lambda[k]
                  << std::setw(17) << lamRef[k] << std::scientific << std::setprecision(2) << std::setw(13)
                  << (r.lambda[k] - lamRef[k]) / lamRef[k] << std::defaultfloat << "\n";
    std::cout << "  max |erro rel.| dos modos 1-n:";
    for (int n : {20, 50, 100, 200, 500}) {
        if (n > (int)r.lambda.size()) break;
        double e = 0.;
        for (int k = 0; k < n; k++) e = std::max(e, std::fabs(r.lambda[k] - lamRef[k]) / lamRef[k]);
        std::cout << "  n = " << n << ": " << std::scientific << std::setprecision(2) << e << std::defaultfloat;
    }
    std::cout << "\n  fracao de variancia sum_{k<=M} lambda_k / 250 (numerica / forma fechada):";
    for (int n : {20, 50, 200, 1000}) {
        if (n > (int)r.lambda.size()) break;
        double sn = 0., sf = 0.;
        for (int k = 0; k < n; k++) { sn += r.lambda[k]; sf += lamRef[k]; }
        std::cout << "  M = " << n << ": " << std::fixed << std::setprecision(5) << sn / kAreaRetangulo << " / "
                  << sf / kAreaRetangulo << std::defaultfloat;
    }
    std::cout << "\n  soma de todos os " << r.neq << " autovalores / 250 = " << std::setprecision(6) << r.traco / kAreaRetangulo
              << "\n";
    std::cout << std::fixed << std::setprecision(2) << "  tempos (s): pontos " << kl.tPontos << ", B " << kl.tB << ", C "
              << kl.tC << " (" << nthreads << " threads), autovalores " << kl.tAuto << ", total " << r.tTotal
              << "; pico de memoria " << std::setprecision(0) << r.memoria << " MB" << std::defaultfloat << std::endl;
    return r;
}

// ---- (b) talude ----

struct LinhaTalude {
    int M = 0;
    // 0: KL numerica no talude; 1: forma fechada do retangulo avaliada no talude
    double erroMQ[2] = {0., 0.}, s2min[2] = {0., 0.}, s2max[2] = {0., 0.};
    double mediaQuad[2] = {0., 0.}; // media de sigma_M^2 por quadratura nos pontos da regra ORDEM (conferencia)
    double covUnifMax[2] = {0., 0.}, covUnifMed[2] = {0., 0.}, covProxMax[2] = {0., 0.}, covProxMed[2] = {0., 0.};
    int Mequivalente = -1; // menor M do retangulo com erro medio quadratico <= o da KL numerica com M termos
};

struct ResultadoTalude {
    MalhaTalude m;
    int64_t nel = 0, neq = 0;
    double hmax = 0., area = 0., traco = 0.;
    std::vector<double> lambda10;
    std::vector<LinhaTalude> linhas;
    double tTotal = 0., tC = 0., tAuto = 0., tAvalia = 0., memoria = 0.;
    int ordem = 0, ordemV = 0, ordemP = 0;
};

// pontos uniformes no talude (rejeicao no retangulo) e pares proximos (o segundo ponto a ate Lx/10 em x e Ly/4 em y)
struct Sorteio {
    std::mt19937_64 g;
    explicit Sorteio(uint64_t s) : g(s) {}
    double U() { return (double)(g() >> 11) * (1.0 / 9007199254740992.0); } // [0, 1)
    std::array<double, 2> NoTaludeUniforme() {
        for (;;) {
            const double x = kX0 + kDx * U(), y = kY0 + kDy * U();
            if (NoTalude(x, y)) return {x, y};
        }
    }
};

ResultadoTalude ParteB(const MalhaTalude &m, double Lx, double Ly, const std::vector<int> &Ms, const klsep::KL2D &ret,
                       const std::vector<double> &erroRet, int ordem, int ordemV, int ordemP, int nthreads, int npares,
                       uint64_t semente, int nrede) {
    ZeraPicoMemoria();
    const auto t0 = Relogio::now();
    ResultadoTalude r;
    r.m = m;
    const int Mmax = *std::max_element(Ms.begin(), Ms.end());
    KLNumerica kl;
    TPZGeoMesh *gmesh = MalhaTaludeGmsh(m.arquivo);
    // caixa dos nos: o retangulo da forma fechada (como no SlopeRandomField)
    double bx0 = 1.e300, bx1 = -1.e300, by0 = 1.e300, by1 = -1.e300;
    for (int64_t i = 0; i < gmesh->NNodes(); i++) {
        TPZManVector<REAL, 3> c(3);
        gmesh->NodeVec()[i].GetCoordinates(c);
        bx0 = std::min(bx0, (double)c[0]); bx1 = std::max(bx1, (double)c[0]);
        by0 = std::min(by0, (double)c[1]); by1 = std::max(by1, (double)c[1]);
    }
    if (std::fabs(bx0 - kX0) > 1.e-9 || std::fabs(bx1 - kX0 - kDx) > 1.e-9 || std::fabs(by0 - kY0) > 1.e-9 ||
        std::fabs(by1 - kY0 - kDy) > 1.e-9) {
        std::cerr << "a caixa dos nos de " << m.arquivo << " nao e [0,25] x [0,10]\n";
        std::exit(1);
    }
    // a integral da forma fechada no talude (IntegralNoTalude) e os pontos aleatorios (NoTalude) supoem o talude
    // (0,0) (25,0) (25,5) (15,5) (10,10) (0,10): a area dos elementos do solo deve ser 187,5 (lados retos)
    {
        double areaMalha = 0.;
        for (TPZGeoEl *gel : gmesh->ElementVec()) {
            if (!gel || gel->Dimension() != 2 || gel->MaterialId() != kMatSolo || gel->HasSubElement()) continue;
            TPZManVector<REAL, 3> a(3), b(3), c(3);
            gel->Node(0).GetCoordinates(a);
            for (int t = 1; t + 1 < gel->NCornerNodes(); t++) { // leque a partir do no 0
                gel->Node(t).GetCoordinates(b);
                gel->Node(t + 1).GetCoordinates(c);
                areaMalha += 0.5 * std::fabs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
            }
        }
        if (std::fabs(areaMalha - kAreaTalude) > 1.e-9 * kAreaTalude) {
            std::cerr << "a area dos elementos de " << m.arquivo << " e " << areaMalha
                      << ", nao a do talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10) (187,5)\n";
            std::exit(1);
        }
    }
    ResolveKL(kl, gmesh, m.p, Lx, Ly, ordem, ordemV, ordemP, nthreads, Mmax);
    const int64_t nv = kl.NModos();
    r.nel = (int64_t)kl.el.size();
    r.neq = kl.neq;
    r.hmax = kl.hmax;
    r.area = kl.area;
    r.traco = kl.lambda.sum();
    r.ordem = kl.ordem;
    r.ordemV = kl.ordemV;
    r.ordemP = kl.ordemP;
    for (int k = 0; k < 10 && k < kl.lambda.size(); k++) r.lambda10.push_back(kl.lambda(k));
    r.tC = kl.tC;
    r.tAuto = kl.tAuto;
    std::vector<int> Mv; // os M desta malha (M <= numero de equacoes)
    for (int M : Ms)
        if (M <= nv) Mv.push_back(M);
    const auto ta = Relogio::now();

    // ---- sigma_M^2 na rede de pontos de cada elemento (minimo e maximo) e nos pontos da regra ORDEM (media por
    // quadratura, que confere sum lambda_k / area na numerica e a integral exata da forma fechada) ----
    std::vector<double> lam(kl.lambda.data(), kl.lambda.data() + nv);
    const size_t nM = Mv.size();
    std::vector<double> s2minG(nM, 1.e300), s2maxG(nM, -1.e300), s2minR(nM, 1.e300), s2maxR(nM, -1.e300);
    std::vector<double> quadG(nM, 0.), quadR(nM, 0.);
    std::vector<double> xs, ys;     // pontos da rede
    std::vector<double> xq, yq, wq; // pontos e pesos x |det J| da regra ORDEM
    // PONTOS (opcional): sigma_M^2 nos pontos da rede, numerica e forma fechada ([ponto * nM + im]), e os primeiros modos
    const std::string prefixoPontos = envs("PONTOS", "");
    const int64_t nModosPontos = prefixoPontos.empty() ? 0 : std::min<int64_t>(std::max(0, envi("PONTOS_MODOS", 9)), nv);
    std::vector<double> s2RedeG, s2RedeR, modosRede;
    {
        TPZMaterialDataT<REAL> data;
        std::vector<double> phi;
        // soma acumulada de lambda_k phi_k^2 nos M de Mv: chama f(im, s) para cada M
        auto acumula = [&](auto f) {
            double s = 0.;
            size_t im = 0;
            for (int64_t k = 0; k < nv && im < nM; k++) {
                s += lam[k] * phi[k] * phi[k];
                if (k + 1 == Mv[im]) f(im++, s);
            }
        };
        for (const auto &e : kl.el) {
            TPZGeoEl *gel = e.intel->Reference();
            // a regra ORDEM do elemento: as coordenadas de referencia dos pontos de e.pontos
            std::unique_ptr<TPZIntPoints> rg(gel->CreateSideIntegrationRule(gel->NSides() - 1, kl.ordem));
            if (rg->NPoints() != (int)e.pontos.fW.size()) {
                std::cerr << "regra de integracao diferente da usada em C\n";
                std::exit(1);
            }
            PreparaDados(e.intel, data);
            for (const auto &q : RedeReferencia(gel->Type(), nrede)) {
                TPZManVector<REAL, 3> qsi = {q[0], q[1]};
                DadosNoPonto(e.intel, qsi, data);
                ModosNoPonto(kl, e, data, phi);
                acumula([&](size_t im, double v) {
                    s2minG[im] = std::min(s2minG[im], v);
                    s2maxG[im] = std::max(s2maxG[im], v);
                    if (!prefixoPontos.empty()) s2RedeG.push_back(v);
                });
                for (int64_t k = 0; k < nModosPontos; k++) modosRede.push_back(phi[k]);
                xs.push_back(data.x[0]);
                ys.push_back(data.x[1]);
            }
            TPZManVector<REAL, 3> qsi(2);
            for (int ip = 0; ip < rg->NPoints(); ip++) {
                REAL w;
                rg->Point(ip, qsi, w);
                DadosNoPonto(e.intel, qsi, data);
                ModosNoPonto(kl, e, data, phi);
                const double wp = e.pontos.fW[ip];
                acumula([&](size_t im, double v) { quadG[im] += wp * v; });
                xq.push_back(data.x[0]);
                yq.push_back(data.x[1]);
                wq.push_back(wp);
            }
        }
    }
    // forma fechada nos mesmos pontos, em blocos de pontos (as tabelas f_i(x_p) e g_j(y_p) de PreparaPontos)
    auto formaFechada = [&](const std::vector<double> &x, const std::vector<double> &y, auto f) {
        for (size_t p0 = 0; p0 < x.size(); p0 += 4096) {
            const size_t p1 = std::min(x.size(), p0 + 4096);
            const std::vector<double> bx(x.begin() + p0, x.begin() + p1), by(y.begin() + p0, y.begin() + p1);
            const auto P = ret.PreparaPontos(bx, by);
            std::vector<double> s(bx.size(), 0.);
            size_t im = 0;
            for (int k = 0; k < ret.M() && im < nM; k++) {
                const auto &fx = P.fx[P.mapx[ret.ix[k]]], &gy = P.gy[P.mapy[ret.iy[k]]];
                const double l = ret.raiz[k] * ret.raiz[k];
                for (size_t q = 0; q < bx.size(); q++) s[q] += l * fx[q] * fx[q] * gy[q] * gy[q];
                if (k + 1 == Mv[im]) f(im++, p0, s);
            }
        }
    };
    if (!prefixoPontos.empty()) s2RedeR.assign(xs.size() * nM, 0.);
    formaFechada(xs, ys, [&](size_t im, size_t p0, const std::vector<double> &v) {
        s2minR[im] = std::min(s2minR[im], *std::min_element(v.begin(), v.end()));
        s2maxR[im] = std::max(s2maxR[im], *std::max_element(v.begin(), v.end()));
        if (!s2RedeR.empty())
            for (size_t q = 0; q < v.size(); q++) s2RedeR[(p0 + q) * nM + im] = v[q];
    });
    formaFechada(xq, yq, [&](size_t im, size_t p0, const std::vector<double> &v) {
        for (size_t q = 0; q < v.size(); q++) quadR[im] += wq[p0 + q] * v[q];
    });

    // ---- covariancia truncada em pares de pontos aleatorios ----
    Sorteio sorteio(semente);
    std::vector<std::array<double, 2>> pts; // 4 npares pontos: pares uniformes (0..2n-1) e proximos (2n..4n-1)
    for (int i = 0; i < npares; i++) {
        pts.push_back(sorteio.NoTaludeUniforme());
        pts.push_back(sorteio.NoTaludeUniforme());
    }
    for (int i = 0; i < npares; i++) {
        const auto a = sorteio.NoTaludeUniforme();
        std::array<double, 2> b;
        do {
            b = {a[0] + Lx / 10. * (2. * sorteio.U() - 1.), a[1] + Ly / 4. * (2. * sorteio.U() - 1.)};
        } while (!NoTalude(b[0], b[1]));
        pts.push_back(a);
        pts.push_back(b);
    }
    const size_t np = pts.size();
    std::vector<std::vector<double>> phiG(np); // modos numericos nos pontos
    {
        TPZMaterialDataT<REAL> data;
        for (size_t i = 0; i < np; i++) {
            int e = -1;
            TPZManVector<REAL, 3> qsi(2);
            if (!Localiza(kl, pts[i][0], pts[i][1], e, qsi)) {
                std::cerr << "ponto (" << pts[i][0] << ", " << pts[i][1] << ") fora da malha " << m.arquivo << "\n";
                std::exit(1);
            }
            PreparaDados(kl.el[e].intel, data);
            DadosNoPonto(kl.el[e].intel, qsi, data);
            if (std::hypot(data.x[0] - pts[i][0], data.x[1] - pts[i][1]) > 1.e-8) {
                std::cerr << "ponto localizado com erro: (" << pts[i][0] << ", " << pts[i][1] << ")\n";
                std::exit(1);
            }
            ModosNoPonto(kl, kl.el[e], data, phiG[i]);
        }
    }
    std::vector<double> px(np), py(np);
    for (size_t i = 0; i < np; i++) { px[i] = pts[i][0]; py[i] = pts[i][1]; }
    const auto PR = ret.PreparaPontos(px, py);
    r.linhas.resize(Mv.size());
    for (size_t im = 0; im < Mv.size(); im++) {
        LinhaTalude &L = r.linhas[im];
        L.M = Mv[im];
        L.s2min[0] = s2minG[im];
        L.s2max[0] = s2maxG[im];
        L.s2min[1] = s2minR[im];
        L.s2max[1] = s2maxR[im];
        L.mediaQuad[0] = quadG[im] / kl.area;
        L.mediaQuad[1] = quadR[im] / kl.area;
    }
    for (int conjunto = 0; conjunto < 2; conjunto++) {
        for (int i = 0; i < npares; i++) {
            const size_t a = (size_t)conjunto * 2 * npares + 2 * i, b = a + 1;
            const double G = std::exp(-std::fabs(px[a] - px[b]) / Lx - std::fabs(py[a] - py[b]) / Ly);
            double cg = 0., cr = 0.;
            size_t im = 0;
            for (int k = 0; k < Mmax && im < Mv.size(); k++) {
                if (k < nv) cg += lam[k] * phiG[a][k] * phiG[b][k];
                const auto &f = PR.fx[PR.mapx[ret.ix[k]]], &g = PR.gy[PR.mapy[ret.iy[k]]];
                cr += ret.raiz[k] * ret.raiz[k] * f[a] * g[a] * f[b] * g[b];
                if (k + 1 == Mv[im]) {
                    LinhaTalude &L = r.linhas[im];
                    const double eg = std::fabs(cg - G), er = std::fabs(cr - G);
                    if (conjunto == 0) {
                        L.covUnifMax[0] = std::max(L.covUnifMax[0], eg);
                        L.covUnifMax[1] = std::max(L.covUnifMax[1], er);
                        L.covUnifMed[0] += eg / npares;
                        L.covUnifMed[1] += er / npares;
                    } else {
                        L.covProxMax[0] = std::max(L.covProxMax[0], eg);
                        L.covProxMax[1] = std::max(L.covProxMax[1], er);
                        L.covProxMed[0] += eg / npares;
                        L.covProxMed[1] += er / npares;
                    }
                    im++;
                }
            }
        }
    }
    // erro medio quadratico: numerica 1 - sum_{k<=M} lambda_k / area (int phi_k^2 = 1 no talude); forma fechada
    // 1 - sum_{k<=M} lambda_k int_talude phi_k^2 / area, calculado em main (erroRet)
    for (auto &L : r.linhas) {
        double s = 0.;
        for (int k = 0; k < L.M; k++) s += lam[k];
        L.erroMQ[0] = 1. - s / kAreaTalude;
        L.erroMQ[1] = erroRet[L.M - 1];
        for (size_t k = 0; k < erroRet.size(); k++)
            if (erroRet[k] <= L.erroMQ[0]) { L.Mequivalente = (int)k + 1; break; }
    }
    r.tAvalia = Segundos(ta);
    r.tTotal = Segundos(t0);
    r.memoria = PicoMemoriaMB();
    // saidas opcionais (depois das medidas de tempo e de memoria)
    const std::string sufixo = "_" + Base(m.arquivo) + "_p" + std::to_string(m.p) + "_o" + std::to_string(kl.ordem);
    const std::string prefixoAuto = envs("AUTOVALORES", "");
    if (!prefixoAuto.empty()) GravaAutovalores(prefixoAuto + sufixo + ".txt", kl.lambda);
    if (!prefixoPontos.empty()) {
        if (s2RedeG.size() != xs.size() * nM || modosRede.size() != xs.size() * (size_t)nModosPontos) {
            std::cerr << "PONTOS: numero de valores inesperado\n";
            std::exit(1);
        }
        std::ofstream f(prefixoPontos + sufixo + ".csv");
        f << "x,y";
        for (int M : Mv) f << ",s2_num_M" << M;
        for (int M : Mv) f << ",s2_ret_M" << M;
        for (int64_t k = 0; k < nModosPontos; k++) f << ",phi" << k + 1;
        f << '\n' << std::setprecision(10);
        for (size_t p = 0; p < xs.size(); p++) {
            f << xs[p] << ',' << ys[p];
            for (size_t im = 0; im < nM; im++) f << ',' << s2RedeG[p * nM + im];
            for (size_t im = 0; im < nM; im++) f << ',' << s2RedeR[p * nM + im];
            for (int64_t k = 0; k < nModosPontos; k++) f << ',' << modosRede[p * nModosPontos + k];
            f << '\n';
        }
    }

    std::cout << "\nTalude " << Nome(m.arquivo) << ", ordem " << m.p << ": " << r.nel << " triangulos, " << r.neq
              << " equacoes, maior aresta " << std::setprecision(4) << r.hmax << " m; area (soma dos pesos) "
              << std::setprecision(10) << r.area << " (187,5)\n";
    std::cout << "  regras (ordem/pontos): pares sem no em comum " << kl.ordem << "/" << kl.npontos
              << ", com no em comum " << kl.ordemV << "/" << kl.npontosV << ", elemento com ele mesmo " << kl.ordemP << "/"
              << kl.npontosP << "\n";
    std::cout << "  10 maiores autovalores, KL numerica no talude:";
    for (double l : r.lambda10) std::cout << " " << std::setprecision(6) << l;
    std::cout << "\n  10 maiores autovalores, forma fechada no retangulo:";
    for (int k = 0; k < 10; k++) std::cout << " " << std::setprecision(6) << ret.raiz[k] * ret.raiz[k];
    std::cout << "\n  soma de todos os " << r.neq << " autovalores numericos / 187,5 = " << std::setprecision(6)
              << r.traco / kAreaTalude << "\n";
    std::cout << "  conferencia, media de sigma_M^2 no talude por quadratura (regra ORDEM) contra o valor usado:";
    for (const auto &L : r.linhas)
        std::cout << "  M = " << L.M << ": numerica " << std::setprecision(8) << L.mediaQuad[0] << " / "
                  << 1. - L.erroMQ[0] << ", retangulo " << L.mediaQuad[1] << " / " << 1. - L.erroMQ[1];
    std::cout << "\n";
    std::cout << "  M     KL               1-media(s2)  min(s2)    max(s2)    cov. pares uniformes   cov. pares proximos\n"
              << "                                                         max|erro|  media      max|erro|  media\n";
    for (const auto &L : r.linhas)
        for (int q = 0; q < 2; q++)
            std::cout << "  " << std::left << std::setw(6) << (q == 0 ? std::to_string(L.M) : std::string(""))
                      << std::setw(17) << (q == 0 ? "numerica talude" : "retangulo") << std::right << std::fixed
                      << std::setprecision(5) << std::setw(11) << L.erroMQ[q] << std::setw(11) << L.s2min[q]
                      << std::setw(11) << L.s2max[q] << std::scientific << std::setprecision(2) << std::setw(11)
                      << L.covUnifMax[q] << std::setw(11) << L.covUnifMed[q] << std::setw(11) << L.covProxMax[q]
                      << std::setw(11) << L.covProxMed[q] << std::defaultfloat << "\n";
    std::cout << "  termos do retangulo com o mesmo erro medio quadratico da KL numerica com M termos:";
    for (const auto &L : r.linhas)
        std::cout << "  M = " << L.M << ": " << (L.Mequivalente > 0 ? std::to_string(L.Mequivalente)
                                                                     : "> " + std::to_string(erroRet.size()));
    std::cout << "\n" << std::fixed << std::setprecision(2) << "  tempos (s): pontos " << kl.tPontos << ", B " << kl.tB
              << ", C " << kl.tC << " (" << nthreads << " threads), autopares " << kl.tAuto << ", avaliacoes "
              << r.tAvalia << ", total " << r.tTotal << "; pico de memoria " << std::setprecision(0) << r.memoria << " MB"
              << std::defaultfloat << std::endl;
    return r;
}

} // namespace

int main() {
    const std::string parte = envs("PARTE", "ab");
    const double Lx = envd("LX", 20.), Ly = envd("LY", 2.);
    const int ordem = envi("ORDEM", 0), ordemV = envi("ORDEM_VIZINHO", 0), ordemP = envi("ORDEM_PROPRIO", 0);
    const int nthreads = envi("THREADS", 2);
    const std::string csvRet = envs("CSV_RET", ""), csvTalude = envs("CSV_TALUDE", "");
    std::cout << "Karhunen-Loeve de Galerkin (NeoPZ, TPZKarhunenLoeveMat tipo 3) contra a forma fechada do retangulo\n"
              << "  G = exp(-|x - x'| / Lx - |y - y'| / Ly), Lx = " << Lx << " m, Ly = " << Ly << " m, sigma^2 = 1\n";

    if (parte.find('a') != std::string::npos) {
        const auto malhas = LeRetangulos(envs("RETANGULOS", "25x10:1,50x20:1,100x40:1,25x10:2,50x20:2"));
        // forma fechada: os 1000 maiores produtos lambda_i^x lambda_j^y
        const klsep::KL2D ref(kX0, kY0, kDx, kDy, Lx, Ly, 1., 1000);
        std::vector<double> lamRef;
        for (int k = 0; k < ref.M(); k++) lamRef.push_back(ref.raiz[k] * ref.raiz[k]);
        std::cout << "\n(a) Retangulo [0,25] x [0,10]: autovalores numericos contra a forma fechada (produtos dos "
                     "autovalores 1D)\n";
        std::vector<ResultadoRet> res;
        for (const auto &m : malhas) {
            res.push_back(ParteA(m, Lx, Ly, lamRef, ordem, ordemV, ordemP, nthreads));
            const auto &r = res.back();
            if (!csvRet.empty()) {
                const bool novo = !std::ifstream(csvRet).good();
                std::ofstream f(csvRet, std::ios::app);
                if (novo) {
                    f << "nx,ny,ordem,elementos,equacoes,h,ordem_regra,ordem_vizinho,ordem_proprio";
                    for (int k = 1; k <= 20; k++) f << ",lambda" << k;
                    f << ",soma_20,soma_200,soma_1000,soma_todos,t_C,t_autovalores,t_total,threads,memoria_MB\n";
                }
                f << m.nx << ',' << m.ny << ',' << m.p << ',' << r.nel << ',' << r.neq << ',' << std::setprecision(6)
                  << r.hmax << ',' << r.ordem << ',' << r.ordemV << ',' << r.ordemP << std::setprecision(12);
                for (int k = 0; k < 20; k++) f << ',' << (k < (int)r.lambda.size() ? r.lambda[k] : NAN);
                for (int n : {20, 200, 1000}) {
                    double s = 0.;
                    for (int k = 0; k < n && k < (int)r.lambda.size(); k++) s += r.lambda[k];
                    f << ',' << (n <= (int)r.lambda.size() ? s : NAN);
                }
                f << ',' << r.traco << std::setprecision(6) << ',' << r.tC << ',' << r.tAuto << ',' << r.tTotal << ','
                  << nthreads << ',' << r.memoria << '\n';
            }
        }
        std::cout << "\nResumo de (a) (erros relativos contra a forma fechada; fracoes de variancia sum lambda / 250)\n"
                  << "  malha       ordem  equacoes  lambda_1    lambda_10   lambda_20   max 1-20    max 1-100   "
                     "M=200 num/ff      todos/250  tempo (s)\n";
        double f200 = 0.;
        for (int k = 0; k < 200; k++) f200 += lamRef[k];
        for (const auto &r : res) {
            auto erro = [&](int k) { return k < (int)r.lambda.size() ? (r.lambda[k] - lamRef[k]) / lamRef[k] : NAN; };
            auto emax = [&](int n) {
                if (n > (int)r.lambda.size()) return (double)NAN;
                double e = 0.;
                for (int k = 0; k < n; k++) e = std::max(e, std::fabs(erro(k)));
                return e;
            };
            double s200 = 0.;
            for (int k = 0; k < 200 && k < (int)r.lambda.size(); k++) s200 += r.lambda[k];
            std::ostringstream nome;
            nome << r.m.nx << 'x' << r.m.ny;
            std::cout << "  " << std::left << std::setw(12) << nome.str() << std::right << std::setw(5) << r.m.p
                      << std::setw(10) << r.neq << std::scientific << std::setprecision(2) << std::setw(12) << erro(0)
                      << std::setw(12) << erro(9) << std::setw(12) << erro(19) << std::setw(12) << emax(20)
                      << std::setw(12) << emax(100) << std::fixed << std::setprecision(5) << std::setw(9)
                      << (r.lambda.size() >= 200 ? s200 / kAreaRetangulo : NAN) << "/" << std::setw(7)
                      << f200 / kAreaRetangulo << std::setw(11) << r.traco / kAreaRetangulo << std::setprecision(1)
                      << std::setw(11) << r.tTotal << std::defaultfloat << "\n";
        }
    }

    if (parte.find('b') != std::string::npos) {
        const auto malhas = LeTaludes(envs("TALUDES", "ceron_h5_a45.msh:1,ceron_h5_a45_h025.msh:1,"
                                                      "ceron_h5_a45_fina.msh:1,ceron_h5_a45.msh:2,ceron_h5_a45_h025.msh:2"),
                                      envs("MALHAS_DIR", MALHAS_DIR));
        std::vector<int> Ms;
        for (const auto &t : Separa(envs("M", "50,200,1000,2000"), ',')) Ms.push_back(std::atoi(t.c_str()));
        std::sort(Ms.begin(), Ms.end());
        Ms.erase(std::unique(Ms.begin(), Ms.end()), Ms.end());
        if (Ms.empty() || Ms.front() < 1) {
            std::cerr << "M invalido\n";
            return 2;
        }
        const int npares = envi("PARES", 500), nrede = envi("REDE", 4);
        const uint64_t semente = (uint64_t)envi("SEMENTE", 1);
        const int retMax = std::max(envi("RET_MAX", 20000), Ms.back());
        const auto tr = Relogio::now();
        // forma fechada no retangulo [0,25] x [0,10] com os Ms.back() maiores termos (a do SlopeRandomField com
        // KL_M = Ms.back()) e, para o M equivalente, com retMax termos; os primeiros M termos sao os mesmos nas duas
        const klsep::KL2D ret(kX0, kY0, kDx, kDy, Lx, Ly, 1., Ms.back());
        const klsep::KL2D retLongo(kX0, kY0, kDx, kDy, Lx, Ly, 1., retMax);
        for (int k = 0; k < ret.M(); k++)
            if (ret.ix[k] != retLongo.ix[k] || ret.iy[k] != retLongo.iy[k]) {
                std::cerr << "os termos da forma fechada mudaram com o numero de termos\n";
                return 1;
            }
        // erro medio quadratico da forma fechada no talude: 1 - sum_{k<=M} lambda_k int_talude phi_k^2 / 187,5
        std::vector<double> erroRet(retLongo.M());
        {
            std::vector<double> gx, gw;
            GaussLegendre(10, gx, gw);
            double s = 0.;
            for (int k = 0; k < retLongo.M(); k++) {
                s += retLongo.raiz[k] * retLongo.raiz[k] * IntegralNoTalude(retLongo, retLongo.ix[k], retLongo.iy[k], gx, gw);
                erroRet[k] = 1. - s / kAreaTalude;
            }
        }
        std::cout << "\n(b) Talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10), area 187,5 m2: KL numerica no talude contra a "
                     "forma fechada do retangulo [0,25] x [0,10] avaliada no talude\n"
                  << "  forma fechada: fracao de variancia retida no retangulo (sum lambda / 250) com M termos:";
        for (int M : Ms) {
            double s = 0.;
            for (int k = 0; k < M; k++) s += retLongo.raiz[k] * retLongo.raiz[k];
            std::cout << "  M = " << M << ": " << std::fixed << std::setprecision(5) << s / kAreaRetangulo
                      << std::defaultfloat;
        }
        std::cout << "\n  forma fechada: 1 - media de sigma_M^2 no talude:";
        for (int M : Ms) std::cout << "  M = " << M << ": " << std::fixed << std::setprecision(5) << erroRet[M - 1];
        std::cout << std::defaultfloat << " (" << std::setprecision(3) << Segundos(tr) << " s, " << retMax
                  << " termos para o M equivalente)" << std::endl;

        std::vector<ResultadoTalude> res;
        for (const auto &m : malhas) {
            res.push_back(ParteB(m, Lx, Ly, Ms, ret, erroRet, ordem, ordemV, ordemP, nthreads, npares, semente, nrede));
            const auto &r = res.back();
            if (!csvTalude.empty()) {
                const bool novo = !std::ifstream(csvTalude).good();
                std::ofstream f(csvTalude, std::ios::app);
                if (novo)
                    f << "malha,ordem,elementos,equacoes,h,ordem_regra,ordem_vizinho,ordem_proprio,M,kl,erro_mq,s2_min,"
                         "s2_max,cov_unif_max,cov_unif_media,cov_prox_max,cov_prox_media,M_ret_equivalente,traco_area,"
                         "t_C,t_autopares,t_total,threads,memoria_MB\n";
                for (const auto &L : r.linhas)
                    for (int q = 0; q < 2; q++)
                        f << Nome(m.arquivo) << ',' << m.p << ',' << r.nel << ',' << r.neq << ',' << std::setprecision(6)
                          << r.hmax << ',' << r.ordem << ',' << r.ordemV << ',' << r.ordemP << ',' << L.M << ','
                          << (q == 0 ? "numerica_talude" : "retangulo") << std::setprecision(10) << ',' << L.erroMQ[q]
                          << ',' << L.s2min[q] << ',' << L.s2max[q] << ',' << L.covUnifMax[q] << ',' << L.covUnifMed[q]
                          << ',' << L.covProxMax[q] << ',' << L.covProxMed[q] << ','
                          << (q == 0 ? L.Mequivalente : L.M) << ',' << r.traco / kAreaTalude << std::setprecision(6)
                          << ',' << r.tC << ',' << r.tAuto << ',' << r.tTotal << ',' << nthreads << ',' << r.memoria
                          << '\n';
            }
        }
        std::cout << "\nResumo de (b): erro medio quadratico de truncamento no talude, 1 - media de sigma_M^2 (numerica "
                     "no talude / forma fechada do retangulo)\n  malha                     ordem  equacoes  todos/187,5";
        for (int M : Ms) std::cout << "   M = " << std::left << std::setw(16) << M << std::right;
        std::cout << "\n";
        for (const auto &r : res) {
            std::cout << "  " << std::left << std::setw(26) << Nome(r.m.arquivo) << std::right << std::setw(5) << r.m.p
                      << std::setw(10) << r.neq << std::fixed << std::setprecision(5) << std::setw(13)
                      << r.traco / kAreaTalude;
            for (int M : Ms) {
                const LinhaTalude *L = nullptr;
                for (const auto &l : r.linhas)
                    if (l.M == M) L = &l;
                if (L) std::cout << "   " << std::setw(8) << L->erroMQ[0] << " / " << std::setw(8) << L->erroMQ[1];
                else std::cout << "   " << std::setw(19) << "(M > equacoes)";
            }
            std::cout << std::defaultfloat << "\n";
        }
    }
    return 0;
}
