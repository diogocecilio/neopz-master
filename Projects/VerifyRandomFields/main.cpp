// VerifyRandomFields - verificacao do campo aleatorio de coesao do talude de Vargas Ceron, Cecilio, Linn e Maghous (2025,
// Secao 6) gerado de tres formas com os MESMOS numeros aleatorios e comparado em todos os pontos de integracao de uma malha
// uniforme fina do talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10), area 187,5 m2:
//   G  Karhunen-Loeve de Galerkin do NeoPZ no proprio talude: TPZKarhunenLoeveMat do tipo 3, B por CalcStiffB, C por
//      ContributeC em todos os pares de elementos (o esquema de ../KLCilindroMika/talude.cpp), elementos H1 de ordem
//      ORDEM_P; C v = lambda B v pelo Cholesky de B e os maiores autopares de L^-1 C L^-T. Os modos sao B-ortonormais.
//   A  forma fechada (analitica) do retangulo [0,25] x [0,10] que contem o talude (a caixa dos nos da malha), a mesma do
//      SlopeRandomField (../SlopeRandomField/KLSeparavel.h, klsep::KL2D): os autopares 2D sao produtos dos autopares 1D do
//      nucleo exponencial; com M termos e com os 2000 termos de producao (coluna A2000).
//   N  Nystrom, independente da malha de Galerkin: o talude e dividido no retangulo de baixo [0,25] x [0,5] e no trapezio
//      de cima {5 <= y <= 10, 0 <= x <= 20 - y} (x = s (20 - y), s em [0,1], jacobiano 20 - y), cada um em celulas de
//      cerca de NYS_CELULA com NYS_Q x NYS_Q pontos de Gauss-Legendre; autopares de W^1/2 K W^1/2 (simetrica),
//      f_k(x_j) = v_jk / sqrt(w_j), e os modos fora dos nos pela interpolacao de Nystrom
//      f_k(x) = (1/lambda_k) sum_j w_j C(x, x_j) f_k(x_j).
//      O nucleo tem dobra em y = y' (e em x = x'); com Ly = 2 m, a regra de Gauss nao integra a dobra em y e os
//      autovalores ficam O(h^2) acima dos exatos (cerca de 1e-2 relativo do modo 10 em diante com celulas de 0,5 m).
//      Com NYS_PRODUTO=1 (padrao), o fator em y dos pares de nos da mesma linha de celulas (as que a dobra y = y' corta)
//      e integrado exatamente (integracao-produto: os nos de Gauss de uma linha de celulas tem as mesmas ordenadas, e a
//      integral dupla de l_a(y) exp(-|y - y'| / Ly) l_b(y') na celula, l_a os polinomios de Lagrange nos pontos de Gauss,
//      e calculada uma vez), o que mantem a matriz simetrica; na interpolacao, o mesmo fator usa
//      int exp(-|y - y'| / Ly) l_b(y') dy' no ponto. NYS_PRODUTO=0: Gauss-Legendre puro.
// Nucleo (Vargas Ceron et al. 2025, eq. 64): C(x, x') = exp(-|x - x'| / Lx - |y - y'| / Ly), variancia 1.
// Campo gaussiano padrao truncado: Z(x) = sum_{k<M} sqrt(lambda_k) xi_k f_k(x); coesao lognormal (klsep::Lognormal):
//   c = exp(mu_ln + sigma_ln Z), sigma_ln^2 = ln(1 + CoV^2), mu_ln = ln(media) - sigma_ln^2 / 2.
// Mesmos numeros aleatorios: xi = klsep::Normal(SEMENTE * 0x9E3779B97F4A7C15 + 3 REAL + 0).Vetor(...), a sequencia da coesao
// do SlopeRandomField na realizacao REAL; os M primeiros valores alimentam os tres metodos.
// Sinais: cada modo de Galerkin tem o sinal que torna int f_k p > 0, p(x, y) = exp((x - 12,5) / 10 + (y - 5) / 4) (um peso
// fixo, positivo e sem simetria); cada modo de Nystrom tem o sinal do produto interno <f_k^G, f_k^N> nos pontos de
// integracao. A forma fechada tem os sinais fixos dos cossenos e senos. A regra e determinista, mas os modos com
// |int f_k p| pequeno (contados no resumo) podem trocar de sinal entre malhas ou valores de ORDEM: so G e N da mesma
// corrida sao comparados ponto a ponto.
// Reprodutibilidade: com o mesmo THREADS, duas corridas dao os mesmos arquivos bit a bit (cada thread monta as linhas
// e1 = t, t + THREADS, ... de C, e as parcelas sao somadas na ordem das threads); com outro THREADS, so o arredondamento
// muda (cerca de 1e-15 relativo em C).
// Atencao: os mesmos xi_k em bases diferentes nao dao a mesma realizacao. G e N resolvem o mesmo problema de autovalores
// no talude e concordam ponto a ponto nos modos bem resolvidos; A usa as autofuncoes do RETANGULO e da outra realizacao de um
// campo com a mesma covariancia no talude: c_G - c_A nao e erro. As comparacoes de estatisticas (sigma_M^2, covariancia
// com pontos de referencia, media e desvio de Monte Carlo) valem para os tres.
//
// Parametros (variaveis de ambiente):
//   MALHA=<fonte>/malhas/talude_uniforme_h025.msh   malha do Gmsh 4.1 (grupos "Solo", "Base", "Direita", "Superficie",
//                         "Esquerda"); sem '/', procurada tambem em <fonte>/malhas. Deve ser o talude acima.
//   LX=20 LY=2            comprimentos de correlacao (m)
//   MEDIA=10 COV=0.3      media (kPa) e coeficiente de variacao da coesao
//   ORDEM_P=1             ordem dos elementos de Galerkin (1 a 3)
//   ORDEM=6               ordem da regra do triangulo nas integrais duplas de C (todos os pares); 0: a do elemento (2P)
//   M=200                 termos da comparacao principal (os tres metodos)
//   M_VERIF=20            termos dos campos extras c_G - c_N (lista: 20,30): truncamento em que os modos de G e N
//                         concordam (modos 1-20 com acordo >= 0,99 nas duas malhas); campos dGN_M<m> nos CSV e VTK
//   MS=10,20,50,100,200,500,1000   termos da varredura (limitados pelos modos disponiveis em G e N)
//   SEMENTE=20250416 REAL=0        semente e realizacao dos mapas
//   R=1000                realizacoes de Monte Carlo (REAL = 0..R-1, os mesmos xi nos tres metodos); 0: sem Monte Carlo
//   REFS=5:2.5,12:6,20:4  pontos de referencia da covariancia (x:y,x:y,...)
//   THREADS=4             threads (montagem de C, produtos de matrizes, interpolacao de Nystrom)
//   NYS_CELULA=0.5 NYS_Q=3 NYS_PRODUTO=1   celulas (m) e pontos de Gauss por direcao do Nystrom; integracao-produto em y
//   SAIDA=saida_<malha>   diretorio das saidas (criado se preciso)
//   VTK=1                 0: sem os arquivos VTK
//   SO_AUTOVALORES=0      1: so os autovalores de G e N (lambda_1..20 e as somas) e termina; para escolher ORDEM e NYS_*
//   VALIDAR=0             1: confere o solver dos maiores autopares contra o SelfAdjointEigenSolver completo do Eigen
// Saidas em SAIDA: pontos.csv (todos os pontos de integracao e campos; x, y, peso e os campos Z com 17 algarismos, o resto
// com 9), modos.csv (autovalores e concordancia dos modos), varredura.csv (diferencas contra M), resumo.txt, e os VTK
// campos_pontos.vtk (nuvem de pontos), campos_elementos.vtk (media de cada campo por elemento) e
// campos_elementos_<campo>.vtk (TPZVTKGeoMesh::PrintGMeshVTK). Ver README.md.
#include <Eigen/Dense>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "TPZElementMatrixT.h"
#include "TPZMaterialDataT.h"
#include "TPZKarhunenLoeveMat.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzquad.h"
#include "tpzintrulet.h"

#include "../SlopeRandomField/KLSeparavel.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#ifndef VRF_DIR
#define VRF_DIR "."
#endif

namespace {

constexpr int kMatSolo = 1;
constexpr double kAreaTalude = 187.5;

// ---------------------------------------------------------------------------------------------------------------------
// utilidades

int envi(const char *k, int d) { const char *v = std::getenv(k); return v ? std::atoi(v) : d; }
double envd(const char *k, double d) { const char *v = std::getenv(k); return v ? std::atof(v) : d; }
std::string envs(const char *k, const std::string &d) { const char *v = std::getenv(k); return v ? std::string(v) : d; }
uint64_t envu(const char *k, uint64_t d) { const char *v = std::getenv(k); return v ? std::strtoull(v, nullptr, 10) : d; }

using Relogio = std::chrono::steady_clock;
double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }

// erro fatal: a mensagem e o fim imediato do processo. std::_Exit nao roda os destrutores estaticos, que poderiam correr
// junto com a thread do Nystrom (std::exit fora da thread principal, ou com ela ainda rodando, e indefinido).
[[noreturn]] void Falha(const std::string &msg, int codigo = 1) {
    std::cout << std::flush;
    std::cerr << msg << std::endl;
    std::_Exit(codigo);
}

// pico da memoria residente do processo (VmHWM), em MB
double PicoMemoriaMB() {
    std::ifstream f("/proc/self/status");
    std::string l;
    while (std::getline(f, l))
        if (l.rfind("VmHWM:", 0) == 0) return std::atof(l.c_str() + 6) / 1024.;
    return -1.;
}

std::vector<std::string> Separa(const std::string &s, char c) {
    std::vector<std::string> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, c))
        if (!t.empty()) v.push_back(t);
    return v;
}

// n pontos de Gauss-Legendre em [0, 1] (Newton nas raizes de P_n)
void GaussLegendre01(int n, std::vector<double> &x, std::vector<double> &w) {
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
        const double wi = 2. / ((1. - z * z) * dp * dp);
        x[i] = 0.5 * (1. - z);
        x[n - 1 - i] = 0.5 * (1. + z);
        w[i] = w[n - 1 - i] = 0.5 * wi;
    }
}

// executa f(i0, i1) em nthreads faixas de [0, n)
void EmFaixas(int64_t n, int nthreads, const std::function<void(int64_t, int64_t)> &f) {
    nthreads = (int)std::max<int64_t>(1, std::min<int64_t>(nthreads, n));
    if (nthreads == 1) { f(0, n); return; }
    std::vector<std::thread> th;
    for (int t = 0; t < nthreads; t++) {
        const int64_t i0 = n * t / nthreads, i1 = n * (t + 1) / nthreads;
        th.emplace_back(f, i0, i1);
    }
    for (auto &t : th) t.join();
}

// Y = X Z, com as linhas de X divididas entre as threads (o produto do Eigen em cada uma)
template <class TX, class TZ> void Produto(const TX &X, const TZ &Z, Eigen::MatrixXd &Y, int nthreads) {
    Y.resize(X.rows(), Z.cols());
    EmFaixas(X.rows(), nthreads, [&](int64_t i0, int64_t i1) {
        if (i1 > i0) Y.middleRows(i0, i1 - i0).noalias() = X.middleRows(i0, i1 - i0) * Z;
    });
}

// ---------------------------------------------------------------------------------------------------------------------
// maiores autopares de uma matriz simetrica densa: tridiagonalizacao de Householder (Eigen, A = Q T Q^T), todos os
// autovalores de T (QR implicito sem vetores), os nv autovetores de T dos maiores autovalores por iteracao inversa (LU de
// T - lambda I com pivoteamento parcial, como dlagtf/dlagts do LAPACK; reortogonalizacao contra os vetores ja calculados com
// autovalor a menos de 1e-7 ||T||) e V = Q Z (reflexoes de Householder em blocos). Custo: 4/3 n^3 + 2 n^2 nv, contra cerca
// de 9 n^3 do SelfAdjointEigenSolver com todos os autovetores.

struct Autopares {
    Eigen::VectorXd lambda; // todos os autovalores, decrescentes
    Eigen::MatrixXd V;      // n x nv, ortonormais
    double residuo = -1.;   // max_k ||A v_k - lambda_k v_k|| / lambda_1
    double ortogonalidade = -1.; // max |V^T V - I|
    double tTri = 0., tVetores = 0., tConfere = 0.;
    bool ok = true;         // false: o QR da tridiagonal nao convergiu (quem chama decide; nada de exit numa thread)
};

Autopares MaioresAutopares(const Eigen::MatrixXd &A, int64_t nv, int nthreads) {
    const auto t0 = Relogio::now();
    const int64_t n = A.rows();
    nv = std::max<int64_t>(0, std::min(nv, n));
    Autopares r;
    Eigen::Tridiagonalization<Eigen::MatrixXd> tri(A);
    const Eigen::VectorXd d = tri.diagonal();
    const Eigen::VectorXd e = tri.subDiagonal();
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.computeFromTridiagonal(d, e, Eigen::EigenvaluesOnly);
        if (es.info() != Eigen::Success) {
            r.ok = false;
            return r;
        }
        r.lambda = es.eigenvalues().reverse();
    }
    r.tTri = Segundos(t0);
    if (nv == 0) return r;
    const auto t1 = Relogio::now();
    double tnorm = 0.;
    for (int64_t i = 0; i < n; i++)
        tnorm = std::max(tnorm, std::fabs(d(i)) + (i > 0 ? std::fabs(e(i - 1)) : 0.) + (i + 1 < n ? std::fabs(e(i)) : 0.));
    const double eps = std::numeric_limits<double>::epsilon();
    const double pivoMin = eps * tnorm, ortol = 1.e-7 * tnorm, pertol = 10. * eps * tnorm;
    Eigen::MatrixXd Z(n, nv);
    std::vector<double> a(n), b(n), c(n), dd(n);
    std::vector<char> troca(n);
    Eigen::VectorXd x(n);
    int64_t grupo = 0;
    double lamAnt = 0.;
    for (int64_t k = 0; k < nv; k++) {
        double lam = r.lambda(k);
        if (k > 0 && lamAnt - lam > ortol) grupo = k;
        if (k > 0 && lamAnt - lam < pertol) lam = lamAnt - pertol; // autovalores (quase) iguais: perturba (dstein)
        lamAnt = lam;
        // T - lam I = P L U: U com diagonal a, superdiagonais b e dd; multiplicadores em c
        for (int64_t i = 0; i < n; i++) a[i] = d(i) - lam;
        for (int64_t i = 0; i + 1 < n; i++) b[i] = c[i] = e(i);
        for (int64_t i = 0; i + 1 < n; i++) {
            if (std::fabs(a[i]) >= std::fabs(c[i])) {
                troca[i] = 0;
                if (a[i] == 0.) a[i] = pivoMin;
                const double mult = c[i] / a[i];
                c[i] = mult;
                a[i + 1] -= mult * b[i];
                dd[i] = 0.;
            } else {
                troca[i] = 1;
                const double mult = a[i] / c[i];
                a[i] = c[i];
                const double temp = a[i + 1];
                a[i + 1] = b[i] - mult * temp;
                if (i + 2 < n) {
                    dd[i] = b[i + 1];
                    b[i + 1] = -mult * dd[i];
                } else {
                    dd[i] = 0.;
                }
                b[i] = temp;
                c[i] = mult;
            }
        }
        auto resolve = [&](Eigen::VectorXd &v) {
            for (int64_t i = 0; i + 1 < n; i++) {
                if (!troca[i]) {
                    v[i + 1] -= c[i] * v[i];
                } else {
                    const double t = v[i];
                    v[i] = v[i + 1];
                    v[i + 1] = t - c[i] * v[i];
                }
            }
            for (int64_t i = n - 1; i >= 0; i--) {
                double s = v[i];
                if (i + 1 < n) s -= b[i] * v[i + 1];
                if (i + 2 < n) s -= dd[i] * v[i + 2];
                double piv = a[i];
                if (std::fabs(piv) < pivoMin) piv = piv < 0. ? -pivoMin : pivoMin;
                v[i] = s / piv;
            }
        };
        std::mt19937_64 g(0x5eed0000ULL + (uint64_t)k);
        std::normal_distribution<double> nd;
        for (int64_t i = 0; i < n; i++) x[i] = nd(g);
        x.normalize();
        for (int it = 0; it < 3; it++) {
            resolve(x);
            x /= x.norm();
            if (k > grupo) {
                auto Zg = Z.middleCols(grupo, k - grupo);
                for (int pass = 0; pass < 2; pass++) x.noalias() -= Zg * (Zg.transpose() * x);
                x.normalize();
            }
        }
        Z.col(k) = x;
    }
    r.V.noalias() = tri.matrixQ() * Z;
    r.tVetores = Segundos(t1);
    // conferencia: residuos e ortogonalidade
    const auto t2 = Relogio::now();
    {
        Eigen::MatrixXd AV;
        Produto(A, r.V, AV, nthreads);
        double res = 0.;
        for (int64_t k = 0; k < nv; k++) res = std::max(res, (AV.col(k) - r.lambda(k) * r.V.col(k)).norm());
        r.residuo = res / std::fabs(r.lambda(0));
        Eigen::MatrixXd VV;
        Eigen::MatrixXd Vt = r.V.transpose();
        Produto(Vt, r.V, VV, nthreads);
        r.ortogonalidade = (VV - Eigen::MatrixXd::Identity(nv, nv)).cwiseAbs().maxCoeff();
    }
    r.tConfere = Segundos(t2);
    return r;
}

// VALIDAR=1: os maiores autopares contra o SelfAdjointEigenSolver completo (autovalores e |<v, v_ref>|)
void ConfereAutopares(const Eigen::MatrixXd &A, const Autopares &r, const std::string &nome, std::ostream &out) {
    const auto t0 = Relogio::now();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Eigen::ComputeEigenvectors);
    const int64_t n = A.rows(), nv = r.V.cols();
    double el = 0., ev = 0.;
    for (int64_t k = 0; k < n; k++) el = std::max(el, std::fabs(r.lambda(k) - es.eigenvalues()(n - 1 - k)));
    for (int64_t k = 0; k < nv; k++) ev = std::max(ev, 1. - std::fabs(r.V.col(k).dot(es.eigenvectors().col(n - 1 - k))));
    out << "  VALIDAR " << nome << ": contra o SelfAdjointEigenSolver completo (" << std::fixed << std::setprecision(1)
        << Segundos(t0) << " s): max |lambda - lambda_ref| / lambda_1 = " << std::scientific << std::setprecision(2)
        << el / std::fabs(r.lambda(0)) << ", max (1 - |<v, v_ref>|) nos " << nv << " vetores = " << ev
        << std::defaultfloat << "\n";
}

// ---------------------------------------------------------------------------------------------------------------------
// malha e KL de Galerkin (adaptado de ../KLCilindroMika/talude.cpp)

TPZGeoMesh *LeMalha(const std::string &arquivo) {
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

struct ElementoKL {
    TPZInterpolationSpace *intel = nullptr;
    std::vector<int64_t> eq;            // equacao global de cada funcao de forma
    double caixa[4] = {0., 0., 0., 0.}; // xmin, xmax, ymin, ymax
    TPZKLPointData pontos;              // dados nos pontos da regra ORDEM (integrais de C)
};

struct KLGalerkin {
    TPZGeoMesh *gmesh = nullptr;
    TPZCompMesh *cmesh = nullptr;
    TPZKarhunenLoeveMat *mat = nullptr;
    std::vector<ElementoKL> el;
    int p = 1, ordem = 0, npontosC = 0;
    int64_t neq = 0;
    double hmax = 0.;
    Eigen::VectorXd lambda; // todos os autovalores, decrescentes
    Eigen::MatrixXd Vt;     // nv x neq: coeficientes dos nv primeiros modos (B-ortonormais), transpostos
    int64_t NModos() const { return Vt.rows(); }
    double tPontos = 0., tB = 0., tC = 0., tReducao = 0., residuo = -1., ortogonalidade = -1.;
    Autopares autopares;    // tempos do solver
    KLGalerkin() = default;
    KLGalerkin(const KLGalerkin &) = delete;
    KLGalerkin &operator=(const KLGalerkin &) = delete;
    ~KLGalerkin() {
        delete cmesh; // a malha computacional apaga o material
    }
};

std::vector<int64_t> Equacoes(TPZCompMesh *cmesh, TPZCompEl *cel) {
    std::vector<int64_t> eq;
    for (int ic = 0; ic < cel->NConnects(); ic++) {
        TPZConnect &c = cel->Connect(ic);
        if (c.HasDependency() || c.IsCondensed()) {
            Falha("conector com dependencia ou condensado: a malha devia ser conforme");
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
        Falha("matriz do par com dimensoes diferentes das equacoes dos elementos");
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

// malha computacional H1 de ordem p no grupo "Solo" e os dados dos elementos (equacoes, caixas e pontos da regra ORDEM)
void PreparaGalerkin(KLGalerkin &kl, TPZGeoMesh *gmesh, int p, double Lx, double Ly, int ordem) {
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
    if (ordem <= 0) ordem = 2 * p;
    if (ordem > TPZIntRuleT::NRULESTRIANGLE_ORDER) {
        Falha("ORDEM " + std::to_string(ordem) + " acima da maxima do triangulo (" +
                  std::to_string(TPZIntRuleT::NRULESTRIANGLE_ORDER) + ")", 2);
    }
    kl.ordem = ordem;
    std::map<int, std::unique_ptr<TPZIntPoints>> regras; // tipo do elemento -> regra de ordem ORDEM
    for (auto *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 2) continue;
        ElementoKL e;
        e.intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        e.eq = Equacoes(cmesh, cel);
        if (!e.intel || (int64_t)e.eq.size() != e.intel->NShapeF()) {
            Falha("elemento inesperado na malha computacional");
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
        auto it = regras.find((int)gel->Type());
        if (it == regras.end())
            it = regras.emplace((int)gel->Type(),
                                std::unique_ptr<TPZIntPoints>(gel->CreateSideIntegrationRule(gel->NSides() - 1, ordem)))
                     .first;
        e.intel->KLPointData(*it->second, e.pontos);
        kl.el.push_back(std::move(e));
    }
    if (!kl.el.empty()) kl.npontosC = (int)kl.el[0].pontos.fW.size();
    kl.tPontos = Segundos(t0);
}

// monta B (CalcStiffB) e C (ContributeC em todos os pares e1 <= e2, em threads), resolve C v = lambda B v pelo Cholesky
// B = L L^T e os maiores autopares de A = L^-1 C L^-T; v = L^-T y
void ResolveGalerkin(KLGalerkin &kl, int nthreads, int64_t nvetores, bool validar, std::ostream &out) {
    const int nel = (int)kl.el.size();
    const int64_t neq = kl.neq;
    const auto t2 = Relogio::now();
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(neq, neq);
    for (auto &e : kl.el) {
        TPZElementMatrixT<REAL> be;
        e.intel->CalcStiffB(be);
        Espalha(be.fMat, e.eq, e.eq, false, B);
    }
    kl.tB = Segundos(t2);
    const auto t3 = Relogio::now();
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(neq, neq);
    {
        nthreads = std::max(1, nthreads);
        std::vector<Eigen::MatrixXd> Cextra(nthreads - 1);
        for (auto &c : Cextra) c = Eigen::MatrixXd::Zero(neq, neq);
        const TPZKarhunenLoeveMat *mat = kl.mat;
        const auto &el = kl.el;
        // a thread t monta as linhas e1 = t, t + THREADS, t + 2 THREADS, ... (distribuicao ciclica fixa, e nao a proxima
        // linha livre): as parcelas de cada thread e a ordem da soma final nao dependem do tempo de cada thread, e C sai
        // igual bit a bit em duas corridas com o mesmo THREADS. Os pares por thread diferem menos de 0,1%.
        auto trabalho = [&](int t, Eigen::MatrixXd &Cl) {
            TPZFNMatrix<400, STATE> ce;
            for (int e1 = t; e1 < nel; e1 += nthreads) {
                const ElementoKL &E1 = el[e1];
                for (int e2 = e1; e2 < nel; e2++) {
                    mat->ContributeC(E1.pontos, el[e2].pontos, ce);
                    Espalha(ce, E1.eq, el[e2].eq, e2 != e1, Cl);
                }
            }
        };
        std::vector<std::thread> th;
        for (int t = 1; t < nthreads; t++) th.emplace_back(trabalho, t, std::ref(Cextra[t - 1]));
        trabalho(0, C);
        for (auto &t : th) t.join();
        for (auto &c : Cextra) C += c; // na ordem das threads
    }
    kl.tC = Segundos(t3);
    const auto t4 = Relogio::now();
    Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(B);
    if (llt.info() != Eigen::Success) {
        Falha("B nao e positiva definida");
    }
    llt.matrixL().solveInPlace(C);
    llt.matrixU().solveInPlace<Eigen::OnTheRight>(C);
    // simetriza (arredondamento das solucoes triangulares)
    for (int64_t j = 0; j < neq; j++)
        for (int64_t i = j + 1; i < neq; i++) C(i, j) = C(j, i) = 0.5 * (C(i, j) + C(j, i));
    kl.tReducao = Segundos(t4);
    kl.autopares = MaioresAutopares(C, nvetores, nthreads);
    if (!kl.autopares.ok) Falha("Galerkin: os autovalores da tridiagonal nao convergiram");
    if (validar) ConfereAutopares(C, kl.autopares, "Galerkin", out);
    C.resize(0, 0);
    kl.lambda = kl.autopares.lambda;
    kl.residuo = kl.autopares.residuo;
    kl.ortogonalidade = kl.autopares.ortogonalidade;
    Eigen::MatrixXd V = std::move(kl.autopares.V);
    kl.autopares.V.resize(0, 0);
    llt.matrixU().solveInPlace(V);
    kl.Vt = V.transpose();
}

void PreparaDados(TPZInterpolationSpace *intel, TPZMaterialDataT<REAL> &data) {
    intel->InitMaterialData(data);
    data.fNeedsSol = false;
    data.p = intel->MaxOrder();
}

// phi_k(x) = sum_s v_k(eq_s) N_s(x), k < nv, a partir das funcoes de forma do elemento no ponto
void ModosNoPonto(const KLGalerkin &kl, const ElementoKL &e, const TPZMaterialDataT<REAL> &data, double *phi) {
    const int64_t nv = kl.NModos(), ns = (int64_t)e.eq.size();
    for (int64_t k = 0; k < nv; k++) phi[k] = 0.;
    for (int64_t s = 0; s < ns; s++) {
        const double N = data.phi(s, 0);
        const double *v = kl.Vt.data() + e.eq[s] * nv;
        for (int64_t k = 0; k < nv; k++) phi[k] += N * v[k];
    }
}

// elemento que contem (x, y) e o ponto qsi correspondente
bool Localiza(const KLGalerkin &kl, double x, double y, int &elemento, TPZManVector<REAL, 3> &qsi) {
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

// os pontos de integracao da malha de Galerkin: em cada elemento computacional, a SUA regra (GetIntegrationRule, a do
// material para a ordem do elemento); peso = w |det J|; os modos de Galerkin pelas funcoes de forma
struct PontosGalerkin {
    std::vector<int64_t> geoel; // indice do elemento geometrico
    std::vector<int> ip;        // indice do ponto na regra do elemento
    std::vector<double> x, y, w;
    Eigen::MatrixXd F;          // np x nv
    int ordemRegra = 0, pontosPorElemento = 0;
    int64_t n() const { return (int64_t)x.size(); }
};

void PontosDeIntegracao(const KLGalerkin &kl, PontosGalerkin &P) {
    int64_t np = 0;
    for (const auto &e : kl.el) np += e.intel->GetIntegrationRule().NPoints();
    const int64_t nv = kl.NModos();
    P.F.resize(np, nv);
    P.x.reserve(np); P.y.reserve(np); P.w.reserve(np); P.geoel.reserve(np); P.ip.reserve(np);
    TPZMaterialDataT<REAL> data;
    std::vector<double> phi(nv);
    int64_t p = 0;
    for (const auto &e : kl.el) {
        const TPZIntPoints &regra = e.intel->GetIntegrationRule();
        if (P.pontosPorElemento == 0) {
            P.pontosPorElemento = regra.NPoints();
            TPZManVector<int, 3> ord(2, 0);
            regra.GetOrder(ord);
            P.ordemRegra = ord[0];
        }
        PreparaDados(e.intel, data);
        TPZManVector<REAL, 3> qsi(2, 0.);
        for (int i = 0; i < regra.NPoints(); i++) {
            REAL wq;
            regra.Point(i, qsi, wq);
            e.intel->ComputeRequiredData(data, qsi);
            ModosNoPonto(kl, e, data, phi.data());
            P.geoel.push_back(e.intel->Reference()->Index());
            P.ip.push_back(i);
            P.x.push_back(data.x[0]);
            P.y.push_back(data.x[1]);
            P.w.push_back(wq * std::fabs(data.detjac));
            for (int64_t k = 0; k < nv; k++) P.F(p, k) = phi[k];
            p++;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Nystrom no talude

struct Nystrom {
    double Lx = 20., Ly = 2.;
    double celula = 0.5;
    int Q = 3;
    bool produto = true;
    int nx = 0, ny = 0, ns = 0; // celulas: retangulo de baixo nx x ny; trapezio de cima ns x ny (s e y)
    double hy = 0.;             // altura das linhas de celulas (5 / ny, a mesma nas duas partes)
    std::vector<double> u, wu;  // Gauss-Legendre em [0, 1]
    std::vector<double> x, y, w, wx;  // nos e pesos; w = wx * wy[b]
    std::vector<double> wy;           // parte em y do peso, por ponto de Gauss em y (hy * wu[b])
    std::vector<int> linha, b;        // linha de celulas (0 .. 2 ny - 1) e ponto de Gauss em y de cada no
    Eigen::MatrixXd G;                // Q x Q: int int l_a(y) exp(-|y - y'| / Ly) l_b(y') dy dy' na celula de altura hy
    std::vector<double> g20, w20;     // Gauss-Legendre de 20 pontos em [0, 1] (integrais em y)
    Eigen::VectorXd lambda;           // todos os autovalores, decrescentes
    Eigen::MatrixXd F;                // n x nv: f_k(x_j)
    Autopares autopares;              // tempos e conferencias
    double tMatriz = 0., traco = 0.;
    int64_t n() const { return (int64_t)x.size(); }

    double Lagrange(int bb, double t) const { // polinomio de Lagrange do ponto bb, t em [0, 1]
        double v = 1.;
        for (int c = 0; c < Q; c++)
            if (c != bb) v *= (t - u[c]) / (u[bb] - u[c]);
        return v;
    }
    // int_0^hy exp(-|t - t'| / Ly) l_b(t' / hy) dt', t em [0, hy] (integrando suave em [0, t] e em [t, hy])
    double Omega(int bb, double t) const {
        double s = 0.;
        for (int parte = 0; parte < 2; parte++) {
            const double a0 = parte == 0 ? 0. : t, a1 = parte == 0 ? t : hy;
            if (a1 <= a0) continue;
            for (size_t q = 0; q < g20.size(); q++) {
                const double tp = a0 + (a1 - a0) * g20[q];
                s += (a1 - a0) * w20[q] * std::exp(-std::fabs(t - tp) / Ly) * Lagrange(bb, tp / hy);
            }
        }
        return s;
    }
    // linha de celulas que contem a ordenada yp e a coordenada local t em [0, hy]
    int LinhaDe(double yp, double &t) const {
        int r;
        double y0;
        if (yp < 5.) {
            r = std::min(ny - 1, std::max(0, (int)std::floor(yp / hy)));
            y0 = r * hy;
        } else {
            const int rr = std::min(ny - 1, std::max(0, (int)std::floor((yp - 5.) / hy)));
            r = ny + rr;
            y0 = 5. + rr * hy;
        }
        t = std::min(hy, std::max(0., yp - y0));
        return r;
    }

    void MontaGrade() {
        GaussLegendre01(Q, u, wu);
        GaussLegendre01(20, g20, w20);
        nx = std::max(1, (int)std::lround(25. / celula));
        ny = std::max(1, (int)std::lround(5. / celula));
        ns = std::max(1, (int)std::ceil(15. / celula - 1.e-9)); // largura 15 m em y = 5, 10 m em y = 10
        hy = 5. / ny;
        wy.resize(Q);
        for (int bb = 0; bb < Q; bb++) wy[bb] = hy * wu[bb];
        x.clear(); y.clear(); w.clear(); wx.clear(); linha.clear(); b.clear();
        // retangulo de baixo [0,25] x [0,5]
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const double x0 = 25. * i / nx, hx = 25. / nx;
                for (int a = 0; a < Q; a++)
                    for (int bb = 0; bb < Q; bb++) {
                        x.push_back(x0 + hx * u[a]);
                        y.push_back(j * hy + hy * u[bb]);
                        wx.push_back(hx * wu[a]);
                        w.push_back(hx * wu[a] * wy[bb]);
                        linha.push_back(j);
                        b.push_back(bb);
                    }
            }
        // trapezio de cima: y em [5, 10], x = s (20 - y), s em [0, 1]
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < ns; i++)
                for (int a = 0; a < Q; a++)
                    for (int bb = 0; bb < Q; bb++) {
                        const double yy = 5. + j * hy + hy * u[bb], s = (i + u[a]) / ns, jac = 20. - yy;
                        x.push_back(s * jac);
                        y.push_back(yy);
                        wx.push_back(wu[a] / ns * jac);
                        w.push_back(wu[a] / ns * jac * wy[bb]);
                        linha.push_back(ny + j);
                        b.push_back(bb);
                    }
        // G(a, b) = int_0^hy l_a(t / hy) Omega_b(t) dt (Omega_b e suave em [0, hy]), simetrizada
        G.resize(Q, Q);
        for (int bb = 0; bb < Q; bb++) {
            std::vector<double> om(g20.size());
            for (size_t q = 0; q < g20.size(); q++) om[q] = Omega(bb, hy * g20[q]);
            for (int a = 0; a < Q; a++) {
                double s = 0.;
                for (size_t q = 0; q < g20.size(); q++) s += hy * w20[q] * Lagrange(a, g20[q]) * om[q];
                G(a, bb) = s;
            }
        }
        G = (0.5 * (G + G.transpose())).eval();
    }

    // elemento (i, j) da matriz simetrica W^1/2 K W^1/2 (com a integracao-produto em y nos pares da mesma linha)
    double Elemento(int64_t i, int64_t j) const {
        const double ex = std::exp(-std::fabs(x[i] - x[j]) / Lx);
        if (produto && linha[i] == linha[j]) return std::sqrt(wx[i] * wx[j]) * ex * G(b[i], b[j]) / std::sqrt(wy[b[i]] * wy[b[j]]);
        return std::sqrt(w[i] * w[j]) * ex * std::exp(-std::fabs(y[i] - y[j]) / Ly);
    }

    // false se os autovalores nao convergiram (roda numa thread: quem chama junta a thread e termina)
    bool Resolve(int64_t nv, int nthreads, bool validar, std::ostream &out) {
        const auto t0 = Relogio::now();
        const int64_t nn = n();
        Eigen::MatrixXd A(nn, nn);
        EmFaixas(nn, nthreads, [&](int64_t j0, int64_t j1) {
            for (int64_t j = j0; j < j1; j++)
                for (int64_t i = 0; i < nn; i++) A(i, j) = Elemento(i, j);
        });
        traco = A.trace();
        tMatriz = Segundos(t0);
        autopares = MaioresAutopares(A, nv, nthreads);
        if (!autopares.ok) return false;
        if (validar) ConfereAutopares(A, autopares, "Nystrom", out);
        lambda = autopares.lambda;
        F = std::move(autopares.V);
        autopares.V.resize(0, 0);
        if (F.rows() == nn) // sem vetores (SO_AUTOVALORES), F fica vazia
            for (int64_t j = 0; j < nn; j++) F.row(j) /= std::sqrt(w[j]);
        return true;
    }

    // linhas [p0, p1) da matriz de interpolacao: P(p, j) = w_j C(x_p, x_j), ou wx_j exp(-|x_p - x_j| / Lx) Omega_b(t_p) se
    // j esta na linha de celulas de y_p (integracao-produto); f_k(x_p) = (1 / lambda_k) sum_j P(p, j) f_k(x_j)
    void Interpolacao(const std::vector<double> &xp, const std::vector<double> &yp, int64_t p0, int64_t p1,
                      Eigen::MatrixXd &P) const {
        const int64_t nn = n();
        P.resize(p1 - p0, nn);
        std::vector<double> om(Q);
        for (int64_t p = p0; p < p1; p++) {
            double t = 0.;
            const int r = LinhaDe(yp[p], t);
            if (produto)
                for (int bb = 0; bb < Q; bb++) om[bb] = Omega(bb, t);
            for (int64_t j = 0; j < nn; j++) {
                const double ex = std::exp(-std::fabs(xp[p] - x[j]) / Lx);
                P(p - p0, j) = produto && linha[j] == r ? wx[j] * ex * om[b[j]]
                                                        : w[j] * ex * std::exp(-std::fabs(yp[p] - y[j]) / Ly);
            }
        }
    }

    // modos nos pontos (xp, yp): Fp (np x nv), em blocos de pontos distribuidos entre as threads
    void ModosNosPontos(const std::vector<double> &xp, const std::vector<double> &yp, Eigen::MatrixXd &Fp,
                        int nthreads, int64_t nvUsar = -1) const {
        const int64_t np = (int64_t)xp.size(), nv = nvUsar < 0 ? F.cols() : std::min<int64_t>(nvUsar, F.cols());
        Fp.resize(np, nv);
        const int64_t bloco = 1024;
        const int64_t nblocos = (np + bloco - 1) / bloco;
        std::atomic<int64_t> proximo{0};
        auto trabalho = [&]() {
            Eigen::MatrixXd P;
            for (int64_t ib; (ib = proximo++) < nblocos;) {
                const int64_t p0 = ib * bloco, p1 = std::min(np, p0 + bloco);
                Interpolacao(xp, yp, p0, p1, P);
                Fp.middleRows(p0, p1 - p0).noalias() = P * F.leftCols(nv);
            }
        };
        std::vector<std::thread> th;
        for (int t = 1; t < std::max(1, nthreads); t++) th.emplace_back(trabalho);
        trabalho();
        for (auto &t : th) t.join();
        for (int64_t k = 0; k < nv; k++) Fp.col(k) /= lambda(k);
    }
};

// ---------------------------------------------------------------------------------------------------------------------
// estatisticas ponderadas pela area

struct Estat {
    double media = 0., dp = 0., min = 0., max = 0., rms = 0., maxabs = 0.;
};

Estat Estatistica(const std::vector<double> &w, const Eigen::Ref<const Eigen::VectorXd> &v) {
    Estat e;
    double sa = 0., s1 = 0., s2 = 0.;
    e.min = 1.e300;
    e.max = -1.e300;
    for (int64_t p = 0; p < v.size(); p++) {
        sa += w[p];
        s1 += w[p] * v(p);
        s2 += w[p] * v(p) * v(p);
        e.min = std::min(e.min, v(p));
        e.max = std::max(e.max, v(p));
    }
    e.media = s1 / sa;
    e.rms = std::sqrt(s2 / sa);
    e.dp = std::sqrt(std::max(0., s2 / sa - e.media * e.media));
    e.maxabs = std::max(std::fabs(e.min), std::fabs(e.max));
    return e;
}

double Correlacao(const std::vector<double> &w, const Eigen::VectorXd &a, const Eigen::VectorXd &b) {
    const Estat ea = Estatistica(w, a), eb = Estatistica(w, b);
    double sa = 0., s = 0.;
    for (int64_t p = 0; p < a.size(); p++) {
        sa += w[p];
        s += w[p] * (a(p) - ea.media) * (b(p) - eb.media);
    }
    return s / sa / (ea.dp * eb.dp);
}

std::string Fmt(double v, int prec = 4, bool cientifico = false) {
    std::ostringstream o;
    if (cientifico) o << std::scientific;
    else o << std::fixed;
    o << std::setprecision(prec) << v;
    return o.str();
}

std::string Base(const std::string &caminho) {
    std::string b = caminho.substr(caminho.find_last_of('/') + 1);
    if (b.size() > 4 && b.compare(b.size() - 4, 4, ".msh") == 0) b.resize(b.size() - 4);
    return b;
}

// ---------------------------------------------------------------------------------------------------------------------
// VTK

// nuvem dos pontos de integracao: celulas VTK_VERTEX e um escalar por ponto para cada campo
void GravaVTKPontos(const std::string &arquivo, const PontosGalerkin &P,
                    const std::vector<std::pair<std::string, const Eigen::VectorXd *>> &campos) {
    std::ofstream v(arquivo);
    const int64_t np = P.n();
    v << "# vtk DataFile Version 3.0\nVerifyRandomFields: pontos de integracao da malha de Galerkin\nASCII\n"
      << "DATASET UNSTRUCTURED_GRID\nPOINTS " << np << " double\n" << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int64_t p = 0; p < np; p++) v << P.x[p] << ' ' << P.y[p] << " 0\n";
    v << "CELLS " << np << ' ' << 2 * np << '\n';
    for (int64_t p = 0; p < np; p++) v << "1 " << p << '\n';
    v << "CELL_TYPES " << np << '\n';
    for (int64_t p = 0; p < np; p++) v << "1\n";
    v << "POINT_DATA " << np << '\n' << std::setprecision(8);
    for (const auto &c : campos) {
        v << "SCALARS " << c.first << " double 1\nLOOKUP_TABLE default\n";
        for (int64_t p = 0; p < np; p++) v << (*c.second)(p) << '\n';
    }
}

// media de cada campo por elemento (sum w f / sum w nos pontos do elemento): indice do elemento geometrico -> valor
std::vector<double> MediaPorElemento(const TPZGeoMesh *gmesh, const PontosGalerkin &P, const Eigen::VectorXd &f) {
    std::vector<double> s(gmesh->NElements(), 0.), sw(gmesh->NElements(), 0.);
    for (int64_t p = 0; p < P.n(); p++) {
        s[P.geoel[p]] += P.w[p] * f(p);
        sw[P.geoel[p]] += P.w[p];
    }
    for (size_t e = 0; e < s.size(); e++) s[e] = sw[e] > 0. ? s[e] / sw[e] : -999.;
    return s;
}

// todos os campos por elemento num arquivo so (os triangulos do solo, nos compartilhados)
void GravaVTKElementos(const std::string &arquivo, TPZGeoMesh *gmesh, const PontosGalerkin &P,
                       const std::vector<std::pair<std::string, const Eigen::VectorXd *>> &campos) {
    std::vector<int64_t> els;
    {
        std::set<int64_t> s(P.geoel.begin(), P.geoel.end());
        els.assign(s.begin(), s.end());
    }
    std::ofstream v(arquivo);
    v << "# vtk DataFile Version 3.0\nVerifyRandomFields: media dos campos em cada elemento\nASCII\n"
      << "DATASET UNSTRUCTURED_GRID\nPOINTS " << gmesh->NNodes() << " double\n"
      << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int64_t i = 0; i < gmesh->NNodes(); i++) {
        TPZManVector<REAL, 3> c(3, 0.);
        gmesh->NodeVec()[i].GetCoordinates(c);
        v << c[0] << ' ' << c[1] << " 0\n";
    }
    int64_t tamanho = 0;
    for (int64_t e : els) tamanho += 1 + gmesh->Element(e)->NCornerNodes();
    v << "CELLS " << els.size() << ' ' << tamanho << '\n';
    for (int64_t e : els) {
        TPZGeoEl *g = gmesh->Element(e);
        v << g->NCornerNodes();
        for (int a = 0; a < g->NCornerNodes(); a++) v << ' ' << g->NodeIndex(a);
        v << '\n';
    }
    v << "CELL_TYPES " << els.size() << '\n';
    for (int64_t e : els) v << (gmesh->Element(e)->NCornerNodes() == 3 ? 5 : 9) << '\n';
    v << "CELL_DATA " << els.size() << '\n' << std::setprecision(8);
    v << "SCALARS elemento int 1\nLOOKUP_TABLE default\n";
    for (int64_t e : els) v << e << '\n';
    for (const auto &c : campos) {
        const auto m = MediaPorElemento(gmesh, P, *c.second);
        v << "SCALARS " << c.first << " double 1\nLOOKUP_TABLE default\n";
        for (int64_t e : els) v << m[e] << '\n';
    }
}

} // namespace

// =====================================================================================================================

int main() {
    Eigen::initParallel();
    const auto tInicio = Relogio::now();
    const double Lx = envd("LX", 20.), Ly = envd("LY", 2.);
    const double media = envd("MEDIA", 10.), cov = envd("COV", 0.3);
    const int pOrdem = envi("ORDEM_P", 1), ordemC = envi("ORDEM", 6);
    const int M = envi("M", 200);
    const uint64_t semente = envu("SEMENTE", 20250416ULL);
    const int real = envi("REAL", 0), R = envi("R", 1000);
    const int nthreads = std::max(1, envi("THREADS", 4));
    const double nysCelula = envd("NYS_CELULA", 0.5);
    const int nysQ = envi("NYS_Q", 3), nysProduto = envi("NYS_PRODUTO", 1);
    const int soAutovalores = envi("SO_AUTOVALORES", 0), validar = envi("VALIDAR", 0), gravaVTK = envi("VTK", 1);
    std::vector<int> MS, MV;
    for (const auto &t : Separa(envs("MS", "10,20,50,100,200,500,1000"), ',')) MS.push_back(std::atoi(t.c_str()));
    for (const auto &t : Separa(envs("M_VERIF", "20"), ',')) MV.push_back(std::atoi(t.c_str()));
    std::vector<std::array<double, 2>> refs;
    for (const auto &t : Separa(envs("REFS", "5:2.5,12:6,20:4"), ',')) {
        std::array<double, 2> r{};
        if (std::sscanf(t.c_str(), "%lf:%lf", &r[0], &r[1]) != 2) {
            std::cerr << "REFS invalido: " << t << " (use x:y,x:y,...)\n";
            return 2;
        }
        refs.push_back(r);
    }
    if (pOrdem < 1 || pOrdem > 3 || M < 1 || nysQ < 1 || nysQ > 12 || nysCelula <= 0. || R < 0 || MS.empty() ||
        cov <= 0. || media <= 0. || std::any_of(MV.begin(), MV.end(), [](int m) { return m < 1; })) {
        std::cerr << "parametros invalidos (ORDEM_P 1 a 3, M >= 1, M_VERIF >= 1, NYS_Q 1 a 12, NYS_CELULA > 0, R >= 0, "
                     "COV > 0, MEDIA > 0)\n";
        return 2;
    }
    std::string malha = envs("MALHA", std::string(VRF_DIR) + "/malhas/talude_uniforme_h025.msh");
    if (!std::filesystem::exists(malha) && std::filesystem::exists(std::string(VRF_DIR) + "/malhas/" + malha))
        malha = std::string(VRF_DIR) + "/malhas/" + malha;
    if (!std::filesystem::exists(malha)) {
        std::cerr << "malha nao encontrada: " << malha << "\n";
        return 2;
    }
    const std::string saida = envs("SAIDA", "saida_" + Base(malha));
    std::filesystem::create_directories(saida);
    // numero de modos numericos: o maior M usado
    int Kmax = M;
    for (int m : MS) Kmax = std::max(Kmax, m);
    for (int m : MV) Kmax = std::max(Kmax, m);

    std::ostringstream res; // resumo (tela e resumo.txt)
    std::string resumo;     // o que ja foi mostrado
    res << "VerifyRandomFields: campo de coesao do talude por Galerkin (G), forma fechada do retangulo (A) e Nystrom (N)\n"
        << "  nucleo C = exp(-|dx| / Lx - |dy| / Ly), Lx = " << Lx << " m, Ly = " << Ly << " m; coesao lognormal, media "
        << media << " kPa, CoV " << cov << "\n"
        << "  malha " << malha << "\n";
    std::cout << res.str() << std::flush;
    resumo = res.str();

    // ---- malha ----
    TPZGeoMesh *gmesh = LeMalha(malha);
    std::unique_ptr<TPZGeoMesh> donoMalha(gmesh); // apagada depois da malha computacional (kl e declarado depois)
    double bx0 = 1.e300, bx1 = -1.e300, by0 = 1.e300, by1 = -1.e300;
    for (int64_t i = 0; i < gmesh->NNodes(); i++) {
        TPZManVector<REAL, 3> c(3);
        gmesh->NodeVec()[i].GetCoordinates(c);
        bx0 = std::min(bx0, (double)c[0]); bx1 = std::max(bx1, (double)c[0]);
        by0 = std::min(by0, (double)c[1]); by1 = std::max(by1, (double)c[1]);
    }
    double areaMalha = 0.;
    int64_t ntri = 0;
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->MaterialId() != kMatSolo || gel->HasSubElement()) continue;
        ntri++;
        TPZManVector<REAL, 3> a(3), b(3), c(3);
        gel->Node(0).GetCoordinates(a);
        for (int t = 1; t + 1 < gel->NCornerNodes(); t++) {
            gel->Node(t).GetCoordinates(b);
            gel->Node(t + 1).GetCoordinates(c);
            areaMalha += 0.5 * std::fabs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
        }
    }
    // o Nystrom e a varredura supoem o talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10)
    if (std::fabs(bx0) > 1.e-9 || std::fabs(bx1 - 25.) > 1.e-9 || std::fabs(by0) > 1.e-9 || std::fabs(by1 - 10.) > 1.e-9 ||
        std::fabs(areaMalha - kAreaTalude) > 1.e-9 * kAreaTalude) {
        std::cerr << "a malha nao e o talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10): caixa [" << bx0 << "," << bx1 << "] x ["
                  << by0 << "," << by1 << "], area " << areaMalha << "\n";
        return 1;
    }

    // ---- Nystrom (numa thread, enquanto a KL de Galerkin e montada) ----
    Nystrom nys;
    nys.Lx = Lx;
    nys.Ly = Ly;
    nys.celula = nysCelula;
    nys.Q = nysQ;
    nys.produto = nysProduto != 0;
    nys.MontaGrade();
    const int64_t nvN = soAutovalores ? 0 : std::min<int64_t>(Kmax, nys.n());
    std::ostringstream resNys;
    // sem exit nem excecao escapando da thread: o resultado e conferido depois do join
    bool nysOk = false;
    std::string nysErro = "os autovalores da tridiagonal nao convergiram";
    std::thread threadNys([&]() {
        try {
            nysOk = nys.Resolve(nvN, 1, validar != 0, resNys);
        } catch (const std::exception &ex) {
            nysErro = ex.what();
        } catch (...) {
            nysErro = "excecao desconhecida";
        }
    });

    // ---- Galerkin ----
    KLGalerkin kl;
    PreparaGalerkin(kl, gmesh, pOrdem, Lx, Ly, ordemC);
    const int64_t nvG = soAutovalores ? 0 : std::min<int64_t>(Kmax, kl.neq);
    std::ostringstream resGal;
    ResolveGalerkin(kl, nthreads, nvG, validar != 0, resGal);
    const double tGalerkin = Segundos(tInicio);
    threadNys.join();
    if (!nysOk) {
        std::cerr << "Nystrom: " << nysErro << "\n";
        return 1;
    }
    const double tAmbos = Segundos(tInicio);

    res << "  Galerkin: " << ntri << " triangulos, ordem " << pOrdem << ", " << kl.neq << " equacoes, maior aresta "
        << Fmt(kl.hmax, 3) << " m; regra de C: ordem " << kl.ordem << " (" << kl.npontosC
        << " pontos por triangulo) em todos os pares\n"
        << "    tempos (s): B " << Fmt(kl.tB, 2) << ", C " << Fmt(kl.tC, 2) << " (" << nthreads << " threads), reducao "
        << Fmt(kl.tReducao, 2) << ", tridiagonalizacao " << Fmt(kl.autopares.tTri, 2) << ", " << nvG << " vetores "
        << Fmt(kl.autopares.tVetores, 2) << ", conferencia " << Fmt(kl.autopares.tConfere, 2) << "\n";
    res << "  Nystrom: celulas de cerca de " << nysCelula << " m, " << nysQ << " x " << nysQ
        << " pontos de Gauss-Legendre por celula, " << (nys.produto ? "integracao-produto em y na linha de celulas" : "Gauss-Legendre puro")
        << "; retangulo " << nys.nx << " x " << nys.ny << " celulas, trapezio " << nys.ns << " x " << nys.ny << " (s x y), "
        << nys.n() << " nos, soma dos pesos " << Fmt(std::accumulate(nys.w.begin(), nys.w.end(), 0.), 10) << "\n"
        << "    tempos (s): matriz " << Fmt(nys.tMatriz, 2) << ", tridiagonalizacao " << Fmt(nys.autopares.tTri, 2) << ", "
        << nvN << " vetores " << Fmt(nys.autopares.tVetores, 2) << ", conferencia " << Fmt(nys.autopares.tConfere, 2)
        << " (em paralelo com Galerkin: " << Fmt(tGalerkin, 1) << " s e " << Fmt(tAmbos, 1) << " s)\n";
    res << resGal.str() << resNys.str();
    // autovalores
    const klsep::KL2D ffRef(bx0, by0, bx1 - bx0, by1 - by0, Lx, Ly, 1., std::max(2000, Kmax));
    {
        res << "  autovalores (G: Galerkin no talude; N: Nystrom no talude; A: forma fechada do retangulo [0,25] x [0,10])\n"
            << "     k       lambda_G       lambda_N      (N-G)/G      lambda_A\n";
        for (int k = 0; k < 20 && k < kl.lambda.size() && k < nys.lambda.size(); k++)
            res << "  " << std::setw(4) << k + 1 << std::setw(15) << Fmt(kl.lambda(k), 6) << std::setw(15)
                << Fmt(nys.lambda(k), 6) << std::setw(13) << Fmt((nys.lambda(k) - kl.lambda(k)) / kl.lambda(k), 2, true)
                << std::setw(14) << Fmt(ffRef.raiz[k] * ffRef.raiz[k], 6) << "\n";
        res << "  soma de todos os autovalores / 187,5: G " << Fmt(kl.lambda.sum() / kAreaTalude, 6) << " (" << kl.lambda.size()
            << "), N " << Fmt(nys.lambda.sum() / kAreaTalude, 6) << " (" << nys.lambda.size() << "; traco da matriz / 187,5 = "
            << Fmt(nys.traco / kAreaTalude, 6) << ")\n";
        for (int m : {20, 200, 1000})
            if (m <= kl.lambda.size() && m <= nys.lambda.size())
                res << "  fracao de variancia sum_{k<=" << m << "} lambda_k / 187,5: G " << Fmt(kl.lambda.head(m).sum() / kAreaTalude, 6)
                    << ", N " << Fmt(nys.lambda.head(m).sum() / kAreaTalude, 6) << "\n";
    }
    if (soAutovalores) {
        std::cout << res.str().substr(resumo.size());
        std::ofstream(saida + "/autovalores_resumo.txt") << res.str();
        std::ofstream fa(saida + "/autovalores.csv");
        fa << "k,lambda_G,lambda_N\n" << std::setprecision(12);
        for (int64_t k = 0; k < std::max(kl.lambda.size(), nys.lambda.size()); k++)
            fa << k + 1 << ',' << (k < kl.lambda.size() ? kl.lambda(k) : NAN) << ','
               << (k < nys.lambda.size() ? nys.lambda(k) : NAN) << '\n';
        return 0;
    }
    res << "  conferencia dos autopares (max ||A v - lambda v|| / lambda_1; max |V^T V - I|): G "
        << Fmt(kl.residuo, 2, true) << ", " << Fmt(kl.ortogonalidade, 2, true) << "; N " << Fmt(nys.autopares.residuo, 2, true)
        << ", " << Fmt(nys.autopares.ortogonalidade, 2, true) << "\n";
    std::cout << res.str().substr(resumo.size()) << std::flush;
    resumo = res.str();
    res.str("");

    // ---- pontos de integracao e modos de Galerkin ----
    const auto tp0 = Relogio::now();
    PontosGalerkin P;
    PontosDeIntegracao(kl, P);
    const int64_t np = P.n();
    double somaPesos = 0.;
    for (double v : P.w) somaPesos += v;
    // sinal dos modos de Galerkin: int f_k p > 0
    const auto pesoSinal = [](double x, double y) { return std::exp((x - 12.5) / 10. + (y - 5.) / 4.); };
    Eigen::VectorXd integralSinal(nvG);
    double normaPeso = 0.;
    for (int64_t p = 0; p < np; p++) normaPeso += P.w[p] * pesoSinal(P.x[p], P.y[p]) * pesoSinal(P.x[p], P.y[p]);
    normaPeso = std::sqrt(normaPeso);
    for (int64_t k = 0; k < nvG; k++) {
        double s = 0.;
        for (int64_t p = 0; p < np; p++) s += P.w[p] * P.F(p, k) * pesoSinal(P.x[p], P.y[p]);
        if (s < 0.) {
            P.F.col(k) *= -1.;
            kl.Vt.row(k) *= -1.;
            s = -s;
        }
        integralSinal(k) = s / normaPeso;
    }
    // modos de Galerkin nos pontos de referencia
    const int nref = (int)refs.size();
    Eigen::MatrixXd FGref(nref, nvG);
    {
        TPZMaterialDataT<REAL> data;
        std::vector<double> phi(nvG);
        for (int r = 0; r < nref; r++) {
            int e = -1;
            TPZManVector<REAL, 3> qsi(2);
            if (!Localiza(kl, refs[r][0], refs[r][1], e, qsi)) {
                std::cerr << "ponto de referencia (" << refs[r][0] << ", " << refs[r][1] << ") fora da malha\n";
                return 2;
            }
            PreparaDados(kl.el[e].intel, data);
            kl.el[e].intel->ComputeRequiredData(data, qsi);
            ModosNoPonto(kl, kl.el[e], data, phi.data());
            for (int64_t k = 0; k < nvG; k++) FGref(r, k) = phi[k];
        }
    }
    const double tPontosG = Segundos(tp0);

    // ---- Nystrom nos pontos (e nos de referencia), alinhado a Galerkin ----
    const auto tn0 = Relogio::now();
    Eigen::MatrixXd FN, FNref;
    nys.ModosNosPontos(P.x, P.y, FN, nthreads);
    {
        std::vector<double> rx, ry;
        for (const auto &r : refs) { rx.push_back(r[0]); ry.push_back(r[1]); }
        nys.ModosNosPontos(rx, ry, FNref, 1);
    }
    const int64_t K = std::min(nvG, nvN); // modos comparaveis
    Eigen::VectorXd normaG(K), normaN(K), acordo(K);
    for (int64_t k = 0; k < K; k++) {
        double sGN = 0., sGG = 0., sNN = 0.;
        for (int64_t p = 0; p < np; p++) {
            sGN += P.w[p] * P.F(p, k) * FN(p, k);
            sGG += P.w[p] * P.F(p, k) * P.F(p, k);
            sNN += P.w[p] * FN(p, k) * FN(p, k);
        }
        if (sGN < 0.) {
            FN.col(k) *= -1.;
            FNref.col(k) *= -1.;
            nys.F.col(k) *= -1.;
            sGN = -sGN;
        }
        normaG(k) = std::sqrt(sGG);
        normaN(k) = std::sqrt(sNN);
        acordo(k) = sGN / (normaG(k) * normaN(k));
    }
    // matriz de cossenos entre todos os pares de modos (G_k, N_j): troca de indice e mistura de modos quase degenerados
    Eigen::MatrixXd Cos;
    {
        Eigen::MatrixXd GWt = (P.F.leftCols(K).array().colwise() * Eigen::Map<const Eigen::ArrayXd>(P.w.data(), np)).matrix().transpose();
        Produto(GWt, FN.leftCols(K), Cos, nthreads);
        for (int64_t j = 0; j < K; j++) Cos.col(j) /= normaN(j);
        for (int64_t k = 0; k < K; k++) Cos.row(k) /= normaG(k);
    }
    // conferencia da interpolacao de Nystrom nos proprios nos (modos 1 a 20 e 1 a M)
    double interpNos20 = 0., interpNosM = 0.;
    {
        const int64_t kk = std::min<int64_t>(std::min<int64_t>(M, K), nys.F.cols());
        Eigen::MatrixXd Fn;
        nys.ModosNosPontos(nys.x, nys.y, Fn, nthreads, kk);
        for (int64_t k = 0; k < kk; k++) {
            const double escala = nys.F.col(k).cwiseAbs().maxCoeff();
            const double d = (Fn.col(k) - nys.F.col(k)).cwiseAbs().maxCoeff() / escala;
            if (k < 20) interpNos20 = std::max(interpNos20, d);
            interpNosM = std::max(interpNosM, d);
        }
    }
    const double tNysPontos = Segundos(tn0);

    // ---- forma fechada nos pontos ----
    const auto ta0 = Relogio::now();
    const int Kff = std::max(2000, Kmax);
    const klsep::KL2D &ff = ffRef;
    {
        // os primeiros M termos com M ou com Kff termos sao os mesmos (a expansao do SlopeRandomField com KL_M = M)
        const klsep::KL2D ffM(bx0, by0, bx1 - bx0, by1 - by0, Lx, Ly, 1., M);
        for (int k = 0; k < ffM.M(); k++)
            if (ffM.ix[k] != ff.ix[k] || ffM.iy[k] != ff.iy[k]) {
                std::cerr << "os termos da forma fechada mudaram com o numero de termos\n";
                return 1;
            }
    }
    Eigen::MatrixXd FA(np, Kff), FAref(nref, Kff);
    Eigen::VectorXd lamA(Kff);
    {
        const auto PA = ff.PreparaPontos(P.x, P.y);
        for (int k = 0; k < Kff; k++) {
            const auto &fx = PA.fx[PA.mapx[ff.ix[k]]], &gy = PA.gy[PA.mapy[ff.iy[k]]];
            for (int64_t p = 0; p < np; p++) FA(p, k) = fx[p] * gy[p];
            lamA(k) = ff.raiz[k] * ff.raiz[k];
        }
        std::vector<double> rx, ry;
        for (const auto &r : refs) { rx.push_back(r[0]); ry.push_back(r[1]); }
        const auto PR = ff.PreparaPontos(rx, ry);
        for (int k = 0; k < Kff; k++)
            for (int r = 0; r < nref; r++) FAref(r, k) = PR.fx[PR.mapx[ff.ix[k]]][r] * PR.gy[PR.mapy[ff.iy[k]]][r];
    }
    const double tFF = Segundos(ta0);

    // M da comparacao principal e da varredura, limitados pelos modos de G e N
    const int Mp = (int)std::min<int64_t>(M, K);
    std::vector<int> Ms;
    for (int m : MS) Ms.push_back((int)std::min<int64_t>(std::max(1, m), K));
    std::sort(Ms.begin(), Ms.end());
    Ms.erase(std::unique(Ms.begin(), Ms.end()), Ms.end());
    std::vector<int> Mv; // M_VERIF
    for (int m : MV) Mv.push_back((int)std::min<int64_t>(m, K));
    std::sort(Mv.begin(), Mv.end());
    Mv.erase(std::unique(Mv.begin(), Mv.end()), Mv.end());

    // ---- campos da realizacao REAL ----
    const klsep::Lognormal ln(media, cov);
    const int Kxi = std::max(Kff, Kmax);
    const auto semente_r = [&](int r) { return semente * 0x9E3779B97F4A7C15ULL + 3ULL * (uint64_t)r + 0ULL; };
    const std::vector<double> xi = klsep::Normal(semente_r(real)).Vetor(Kxi);
    Eigen::VectorXd sqG = kl.lambda.head(K).cwiseMax(0.).cwiseSqrt(), sqN = nys.lambda.head(K).cwiseMax(0.).cwiseSqrt();
    Eigen::VectorXd sqA = lamA.cwiseSqrt();
    const Eigen::Map<const Eigen::VectorXd> xiv(xi.data(), Kxi);
    auto campoZ = [&](const Eigen::MatrixXd &F, const Eigen::VectorXd &sq, int m) -> Eigen::VectorXd {
        return F.leftCols(m) * (sq.head(m).cwiseProduct(xiv.head(m)));
    };
    auto coesao = [&](const Eigen::VectorXd &z) -> Eigen::VectorXd { return z.unaryExpr([&](double v) { return ln(v); }); };
    const Eigen::VectorXd zG = campoZ(P.F, sqG, Mp), zN = campoZ(FN, sqN, Mp), zA = campoZ(FA, sqA, Mp),
                          zA2000 = campoZ(FA, sqA, 2000);
    const Eigen::VectorXd cG = coesao(zG), cN = coesao(zN), cA = coesao(zA), cA2000 = coesao(zA2000);
    const Eigen::VectorXd dGN = cG - cN, dGA = cG - cA, dNA = cN - cA, dGA2000 = cG - cA2000;
    // c_G - c_N com M_VERIF termos: o truncamento em que os modos das duas discretizacoes concordam (a verificacao ponto a
    // ponto de G contra N; com M = 200 a diferenca e dominada pelos modos que nenhuma das duas resolve)
    const size_t nMv = Mv.size();
    std::vector<Eigen::VectorXd> zGv(nMv), zNv(nMv), dGNv(nMv);
    std::vector<std::string> nomeMv(nMv);
    for (size_t i = 0; i < nMv; i++) {
        zGv[i] = campoZ(P.F, sqG, Mv[i]);
        zNv[i] = campoZ(FN, sqN, Mv[i]);
        dGNv[i] = coesao(zGv[i]) - coesao(zNv[i]);
        nomeMv[i] = "_M" + std::to_string(Mv[i]);
    }
    // conferencia: o campo A2000 como o SlopeRandomField o calcula (KL2D::Campo)
    double difCampo = 0.;
    {
        klsep::KL2D ff2000(bx0, by0, bx1 - bx0, by1 - by0, Lx, Ly, 1., 2000);
        const auto PA2 = ff2000.PreparaPontos(P.x, P.y);
        const std::vector<double> xi2000(xi.begin(), xi.begin() + 2000);
        const auto z = ff2000.Campo(PA2, xi2000);
        for (int64_t p = 0; p < np; p++) difCampo = std::max(difCampo, std::fabs(z[p] - zA2000(p)));
    }
    // variancias truncadas sigma_M^2 e diferencas esperadas (sobre as realizacoes) dos campos Z
    auto variancia = [&](const Eigen::MatrixXd &F, const Eigen::VectorXd &sq, int m) -> Eigen::VectorXd {
        return (F.leftCols(m) * sq.head(m).asDiagonal()).rowwise().squaredNorm();
    };
    const Eigen::VectorXd s2G = variancia(P.F, sqG, Mp), s2N = variancia(FN, sqN, Mp), s2A = variancia(FA, sqA, Mp),
                          s2A2000 = variancia(FA, sqA, 2000);
    auto difEsperada = [&](const Eigen::MatrixXd &F1, const Eigen::VectorXd &q1, const Eigen::MatrixXd &F2,
                           const Eigen::VectorXd &q2, int m) -> Eigen::VectorXd {
        return (F1.leftCols(m) * q1.head(m).asDiagonal() - F2.leftCols(m) * q2.head(m).asDiagonal()).rowwise().norm();
    };
    const Eigen::VectorXd ezGN = difEsperada(P.F, sqG, FN, sqN, Mp), ezGA = difEsperada(P.F, sqG, FA, sqA, Mp),
                          ezNA = difEsperada(FN, sqN, FA, sqA, Mp);
    // covariancia com os pontos de referencia: sum_{k<M} lambda_k f_k(x) f_k(x_r) contra o nucleo
    std::vector<Eigen::VectorXd> covG(nref), covN(nref), covA(nref), covAlvo(nref);
    for (int r = 0; r < nref; r++) {
        covG[r] = P.F.leftCols(Mp) * (kl.lambda.head(Mp).cwiseProduct(FGref.row(r).head(Mp).transpose()));
        covN[r] = FN.leftCols(Mp) * (nys.lambda.head(Mp).cwiseProduct(FNref.row(r).head(Mp).transpose()));
        covA[r] = FA.leftCols(Mp) * (lamA.head(Mp).cwiseProduct(FAref.row(r).head(Mp).transpose()));
        covAlvo[r].resize(np);
        for (int64_t p = 0; p < np; p++)
            covAlvo[r](p) = std::exp(-std::fabs(P.x[p] - refs[r][0]) / Lx - std::fabs(P.y[p] - refs[r][1]) / Ly);
    }
    // alvos lognormais com a variancia truncada: E[c] = exp(mu + s^2 sigma_M^2 / 2), DP[c] = E[c] sqrt(exp(s^2 sigma_M^2) - 1)
    auto alvoMedia = [&](const Eigen::VectorXd &s2) -> Eigen::VectorXd {
        return s2.unaryExpr([&](double v) { return std::exp(ln.mu_ln + 0.5 * ln.sigma_ln * ln.sigma_ln * v); });
    };
    auto alvoDP = [&](const Eigen::VectorXd &s2) -> Eigen::VectorXd {
        return s2.unaryExpr([&](double v) {
            const double s = ln.sigma_ln * ln.sigma_ln * v;
            return std::exp(ln.mu_ln + 0.5 * s) * std::sqrt(std::max(0., std::exp(s) - 1.));
        });
    };
    const Eigen::VectorXd amG = alvoMedia(s2G), amN = alvoMedia(s2N), amA = alvoMedia(s2A), amA2000 = alvoMedia(s2A2000);
    const Eigen::VectorXd adG = alvoDP(s2G), adN = alvoDP(s2N), adA = alvoDP(s2A), adA2000 = alvoDP(s2A2000);
    // ortonormalidade dos modos nos pontos de integracao, max |F^T W F - I| (W: os pesos w |det J|). Os de Galerkin sao
    // B-ortonormais e a regra do elemento integra P1 x P1 exatamente; os de Nystrom interpolados so sao ortonormais a menos
    // do erro da interpolacao. As normas entram na variancia truncada: a media na area de sigma_M^2 e
    // sum_k lambda_k ||f_k||^2 / 187,5, e nao sum_k lambda_k / 187,5.
    auto maxGram = [&](const Eigen::MatrixXd &F, int64_t m) {
        const Eigen::MatrixXd FWt =
            (F.leftCols(m).array().colwise() * Eigen::Map<const Eigen::ArrayXd>(P.w.data(), np)).matrix().transpose();
        Eigen::MatrixXd Gm;
        Produto(FWt, F.leftCols(m), Gm, nthreads);
        return (Gm - Eigen::MatrixXd::Identity(m, m)).cwiseAbs().maxCoeff();
    };
    const int m20 = std::min(20, Mp);
    const double gramG20 = maxGram(P.F, m20), gramGM = maxGram(P.F, Mp), gramN20 = maxGram(FN, m20), gramNM = maxGram(FN, Mp);
    const double areaRet = (bx1 - bx0) * (by1 - by0);
    const double defLamG = 1. - kl.lambda.head(Mp).sum() / kAreaTalude, defLamN = 1. - nys.lambda.head(Mp).sum() / kAreaTalude;
    const double defLamA = 1. - lamA.head(Mp).sum() / areaRet, defLamA2000 = 1. - lamA.head(2000).sum() / areaRet;
    double normaContribG = 0., normaContribN = 0.; // sum_{k<M} lambda_k (||f_k||^2 - 1) / 187,5
    for (int k = 0; k < Mp; k++) {
        normaContribG += kl.lambda(k) * (normaG(k) * normaG(k) - 1.) / kAreaTalude;
        normaContribN += nys.lambda(k) * (normaN(k) * normaN(k) - 1.) / kAreaTalude;
    }

    // ---- Monte Carlo: R realizacoes, os mesmos xi nos tres metodos ----
    const auto tmc0 = Relogio::now();
    const int nmc = 4; // G, N, A, A2000
    std::vector<Eigen::VectorXd> mcS1(nmc, Eigen::VectorXd::Zero(np)), mcS2(nmc, Eigen::VectorXd::Zero(np));
    std::vector<Eigen::VectorXd> mcD2(3, Eigen::VectorXd::Zero(np)); // (cG - cN)^2, (cG - cA)^2, (cN - cA)^2
    if (R > 0) {
        const Eigen::MatrixXd FGs = P.F.leftCols(Mp) * sqG.head(Mp).asDiagonal(), FNs = FN.leftCols(Mp) * sqN.head(Mp).asDiagonal();
        const Eigen::MatrixXd FAs = FA * sqA.asDiagonal();
        const int bloco = 50;
        for (int r0 = 0; r0 < R; r0 += bloco) {
            const int nb = std::min(bloco, R - r0);
            Eigen::MatrixXd Xi(Kff, nb);
            for (int j = 0; j < nb; j++) {
                const auto v = klsep::Normal(semente_r(r0 + j)).Vetor(Kff);
                for (int k = 0; k < Kff; k++) Xi(k, j) = v[k];
            }
            Eigen::MatrixXd Z[4];
            const Eigen::MatrixXd XiM = Xi.topRows(Mp), Xi2000 = Xi.topRows(2000);
            Produto(FGs, XiM, Z[0], nthreads);
            Produto(FNs, XiM, Z[1], nthreads);
            Produto(FAs.leftCols(Mp), XiM, Z[2], nthreads);
            Produto(FAs.leftCols(2000), Xi2000, Z[3], nthreads);
            EmFaixas(np, nthreads, [&](int64_t p0, int64_t p1) {
                for (int j = 0; j < nb; j++)
                    for (int64_t p = p0; p < p1; p++) {
                        double c[4];
                        for (int m = 0; m < nmc; m++) {
                            c[m] = ln(Z[m](p, j));
                            mcS1[m](p) += c[m];
                            mcS2[m](p) += c[m] * c[m];
                        }
                        mcD2[0](p) += (c[0] - c[1]) * (c[0] - c[1]);
                        mcD2[1](p) += (c[0] - c[2]) * (c[0] - c[2]);
                        mcD2[2](p) += (c[1] - c[2]) * (c[1] - c[2]);
                    }
            });
        }
    }
    std::vector<Eigen::VectorXd> mcMedia(nmc), mcDP(nmc), mcRmsD(3);
    for (int m = 0; m < nmc; m++) {
        mcMedia[m] = mcS1[m] / std::max(1, R);
        mcDP[m] = ((mcS2[m] - mcS1[m].cwiseProduct(mcMedia[m])) / std::max(1, R - 1)).cwiseMax(0.).cwiseSqrt();
    }
    for (int m = 0; m < 3; m++) mcRmsD[m] = (mcD2[m] / std::max(1, R)).cwiseSqrt();
    const double tMC = Segundos(tmc0);

    // ---- varredura em M (realizacao REAL) ----
    struct LinhaVarredura {
        int M = 0;
        Estat dGN, dGA, dNA, cG;
        double ezGN = 0., ezGA = 0., ezNA = 0., erroVarG = 0., erroVarN = 0., erroVarA = 0., corrGN = 0., corrGA = 0.;
        double erroVarLamG = 0., erroVarLamN = 0.; // 1 - sum_{k<M} lambda_k / 187,5
    };
    std::vector<LinhaVarredura> varredura;
    {
        // somas acumuladas em k (uma passada pelas colunas): campos Z da realizacao REAL, quadrados das diferencas
        // esperadas e variancias truncadas, registradas em cada M da varredura
        Eigen::VectorXd zg = Eigen::VectorXd::Zero(np), zn = zg, za = zg, sGN = zg, sGA = zg, sNA = zg, vG = zg, vN = zg, vA = zg;
        double somaLamG = 0., somaLamN = 0.;
        size_t im = 0;
        for (int k = 0; k < Ms.back(); k++) {
            somaLamG += kl.lambda(k);
            somaLamN += nys.lambda(k);
            const double qg = sqG(k), qn = sqN(k), qa = sqA(k);
            for (int64_t p = 0; p < np; p++) {
                const double ag = qg * P.F(p, k), an = qn * FN(p, k), aa = qa * FA(p, k);
                zg(p) += ag * xi[k];
                zn(p) += an * xi[k];
                za(p) += aa * xi[k];
                sGN(p) += (ag - an) * (ag - an);
                sGA(p) += (ag - aa) * (ag - aa);
                sNA(p) += (an - aa) * (an - aa);
                vG(p) += ag * ag;
                vN(p) += an * an;
                vA(p) += aa * aa;
            }
            if (k + 1 != Ms[im]) continue;
            LinhaVarredura L;
            L.M = k + 1;
            const Eigen::VectorXd g = coesao(zg), nn = coesao(zn), a = coesao(za);
            L.dGN = Estatistica(P.w, g - nn);
            L.dGA = Estatistica(P.w, g - a);
            L.dNA = Estatistica(P.w, nn - a);
            L.cG = Estatistica(P.w, g);
            L.corrGN = Correlacao(P.w, g, nn);
            L.corrGA = Correlacao(P.w, g, a);
            L.ezGN = std::sqrt(Estatistica(P.w, sGN).media);
            L.ezGA = std::sqrt(Estatistica(P.w, sGA).media);
            L.ezNA = std::sqrt(Estatistica(P.w, sNA).media);
            L.erroVarG = 1. - Estatistica(P.w, vG).media;
            L.erroVarN = 1. - Estatistica(P.w, vN).media;
            L.erroVarA = 1. - Estatistica(P.w, vA).media;
            L.erroVarLamG = 1. - somaLamG / kAreaTalude;
            L.erroVarLamN = 1. - somaLamN / kAreaTalude;
            varredura.push_back(L);
            if (++im == Ms.size()) break;
        }
    }

    // ---- resumo ----
    res << "\nPontos de integracao da malha de Galerkin (regra de cada elemento: ordem " << P.ordemRegra << ", "
        << P.pontosPorElemento << " pontos por triangulo): " << np << " pontos; soma dos pesos w |det J| = "
        << Fmt(somaPesos, 12) << " (187,5; diferenca " << Fmt(somaPesos - kAreaTalude, 2, true) << ")\n";
    {
        int fracos3 = 0, fracos3M = 0, fracos4 = 0, fracos4M = 0; // |int f_k p| / ||p|| < 1e-3 e < 1e-4 (nos Mp primeiros)
        for (int64_t k = 0; k < nvG; k++) {
            fracos3 += integralSinal(k) < 1.e-3;
            fracos4 += integralSinal(k) < 1.e-4;
            if (k < Mp) {
                fracos3M += integralSinal(k) < 1.e-3;
                fracos4M += integralSinal(k) < 1.e-4;
            }
        }
        res << "  sinal dos modos de Galerkin: int f_k p > 0, p = exp((x - 12,5) / 10 + (y - 5) / 4); menor |int f_k p| / ||p|| nos "
            << nvG << " modos: " << Fmt(integralSinal.cwiseAbs().minCoeff(), 2, true) << "; modos com |int f_k p| / ||p|| < 1e-3: "
            << fracos3 << " (" << fracos3M << " nos " << Mp << " primeiros), < 1e-4: " << fracos4 << " (" << fracos4M
            << ") (sinais que podem mudar com a malha ou ORDEM; G e N da mesma corrida tem sinais alinhados)\n";
    }
    res << "  modos de Nystrom interpolados nos proprios nos contra os valores nodais, max |dif| / max |f_k|: modos 1-20 "
        << Fmt(interpNos20, 2, true) << ", modos 1-" << std::min<int64_t>(M, K) << " " << Fmt(interpNosM, 2, true)
        << (nys.produto ? " (com a integracao-produto, a interpolacao aplica a integral exata em y a representacao polinomial da "
                          "linha de celulas e difere dos valores nodais pelo erro dessa representacao; com NYS_PRODUTO=0 a "
                          "diferenca e so de arredondamento)\n"
                        : " (Gauss-Legendre puro: so arredondamento)\n");
    res << "  forma fechada: o campo A2000 pela matriz de modos contra klsep::KL2D::Campo (o calculo do SlopeRandomField): max |dZ| = "
        << Fmt(difCampo, 2, true) << "\n";
    // concordancia dos modos
    {
        res << "\nConcordancia dos modos de Galerkin e de Nystrom nos pontos de integracao, |<f_k^G, f_k^N>| / (||f_k^G|| ||f_k^N||)"
               " (modos de Nystrom com o sinal alinhado):\n     k  lambda_G     (N-G)/G     acordo       melhor j  |cos(G_k, N_j)|  ||f^N||\n";
        for (int64_t k = 0; k < std::min<int64_t>(K, 30); k++) {
            Eigen::Index jm = 0;
            Cos.row(k).cwiseAbs().maxCoeff(&jm);
            res << "  " << std::setw(4) << k + 1 << std::setw(11) << Fmt(kl.lambda(k), 5) << std::setw(12)
                << Fmt((nys.lambda(k) - kl.lambda(k)) / kl.lambda(k), 2, true) << std::setw(14) << Fmt(acordo(k), 8)
                << std::setw(8) << jm + 1 << std::setw(14) << Fmt(std::fabs(Cos(k, jm)), 6) << std::setw(12)
                << Fmt(normaN(k), 6) << "\n";
        }
        int primeiro999 = -1, primeiro99 = -1, primeiro9 = -1;
        for (int64_t k = 0; k < K; k++) {
            if (primeiro999 < 0 && acordo(k) < 0.999) primeiro999 = (int)k + 1;
            if (primeiro99 < 0 && acordo(k) < 0.99) primeiro99 = (int)k + 1;
            if (primeiro9 < 0 && acordo(k) < 0.9) primeiro9 = (int)k + 1;
        }
        auto conta = [&](int64_t m, double lim) { int c = 0; for (int64_t k = 0; k < m; k++) c += acordo(k) >= lim; return c; };
        res << "  primeiro modo com acordo < 0,999: " << primeiro999 << "; < 0,99: " << primeiro99 << "; < 0,9: " << primeiro9 << "\n";
        for (int m : {20, 50, 100, 200, 500, 1000})
            if (m <= K)
                res << "  modos 1-" << m << ": acordo >= 0,99 em " << conta(m, 0.99) << ", >= 0,9 em " << conta(m, 0.9)
                    << "; menor acordo " << Fmt(acordo.head(m).minCoeff(), 4) << "\n";
    }
    // comparacao principal
    {
        const Estat eG = Estatistica(P.w, cG), eN = Estatistica(P.w, cN), eA = Estatistica(P.w, cA), eA2 = Estatistica(P.w, cA2000);
        res << "\nRealizacao REAL = " << real << " (SEMENTE = " << semente << "), M = " << Mp << " termos (A2000: 2000 termos), "
            << "coesao em kPa, estatisticas ponderadas pela area\n"
            << "  campo   media      desvio     minimo     maximo\n";
        const std::pair<const char *, Estat> lin[] = {{"c_G", eG}, {"c_N", eN}, {"c_A", eA}, {"c_A2000", eA2}};
        for (const auto &l : lin)
            res << "  " << std::left << std::setw(8) << l.first << std::right << std::setw(9) << Fmt(l.second.media, 4)
                << std::setw(11) << Fmt(l.second.dp, 4) << std::setw(11) << Fmt(l.second.min, 4) << std::setw(11)
                << Fmt(l.second.max, 4) << "\n";
        res << "  diferenca        RMS        max |.|    correlacao dos campos\n";
        const std::tuple<const char *, Eigen::VectorXd, double> dif[] = {
            {"c_G - c_N", dGN, Correlacao(P.w, cG, cN)}, {"c_G - c_A", dGA, Correlacao(P.w, cG, cA)},
            {"c_N - c_A", dNA, Correlacao(P.w, cN, cA)}, {"c_G - c_A2000", dGA2000, Correlacao(P.w, cG, cA2000)}};
        for (const auto &d : dif) {
            const Estat e = Estatistica(P.w, std::get<1>(d));
            res << "  " << std::left << std::setw(15) << std::get<0>(d) << std::right << std::setw(10) << Fmt(e.rms, 4, false)
                << std::setw(11) << Fmt(e.maxabs, 4) << std::setw(12) << Fmt(std::get<2>(d), 4) << "\n";
        }
        res << "  diferenca esperada sobre as realizacoes (exata), RMS na area de sqrt(E[(Z_1 - Z_2)^2]) (Z normal padrao): G-N "
            << Fmt(Estatistica(P.w, ezGN).rms, 4) << ", G-A " << Fmt(Estatistica(P.w, ezGA).rms, 4) << ", N-A "
            << Fmt(Estatistica(P.w, ezNA).rms, 4) << "\n";
        if (R > 0)
            res << "  o mesmo para a coesao, por Monte Carlo (R = " << R << "), RMS na area de sqrt(E[(c_1 - c_2)^2]) (kPa): G-N "
                << Fmt(Estatistica(P.w, mcRmsD[0]).rms, 4) << ", G-A " << Fmt(Estatistica(P.w, mcRmsD[1]).rms, 4) << ", N-A "
                << Fmt(Estatistica(P.w, mcRmsD[2]).rms, 4) << "\n";
        for (size_t i = 0; i < nMv; i++) {
            const Estat e = Estatistica(P.w, dGNv[i]);
            const Eigen::VectorXd ez = difEsperada(P.F, sqG, FN, sqN, Mv[i]);
            res << "  c_G - c_N com M = " << Mv[i] << " (M_VERIF; menor acordo dos modos 1-" << Mv[i] << ": "
                << Fmt(acordo.head(Mv[i]).minCoeff(), 6) << "): RMS " << Fmt(e.rms, 5) << ", max |.| " << Fmt(e.maxabs, 4)
                << ", correlacao " << Fmt(Correlacao(P.w, coesao(zGv[i]), coesao(zNv[i])), 8)
                << "; RMS na area de sqrt(E[(Z_G - Z_N)^2]) " << Fmt(Estatistica(P.w, ez).rms, 5) << " (campo dGN"
                << nomeMv[i] << ")\n";
        }
    }
    // variancias
    {
        res << "\nVariancia truncada sigma_M^2(x) = sum_{k<M} lambda_k f_k(x)^2 (alvo 1) nos pontos, M = " << Mp << "\n"
            << "  campo   1 - media  minimo     maximo     RMS(sigma^2 - 1)   1 - sum_{k<M} lambda_k / area\n";
        const std::tuple<const char *, const Eigen::VectorXd *, std::string> lin[] = {
            {"G", &s2G, Fmt(defLamG, 5)},
            {"N", &s2N, Fmt(defLamN, 5)},
            {"A", &s2A, Fmt(defLamA, 5) + " (area do retangulo, " + Fmt(areaRet, 1) + ")"},
            {"A2000", &s2A2000, Fmt(defLamA2000, 5) + " (idem)"}};
        for (const auto &l : lin) {
            const Estat e = Estatistica(P.w, *std::get<1>(l));
            const Estat e1 = Estatistica(P.w, (std::get<1>(l)->array() - 1.).matrix());
            res << "  " << std::left << std::setw(8) << std::get<0>(l) << std::right << std::setw(9) << Fmt(1. - e.media, 5)
                << std::setw(11) << Fmt(e.min, 5) << std::setw(11) << Fmt(e.max, 5) << std::setw(12) << Fmt(e1.rms, 5)
                << "       " << std::get<2>(l) << "\n";
        }
        res << "  1 - media de sigma_M^2 = (1 - sum_{k<M} lambda_k / 187,5) - sum_{k<M} lambda_k (||f_k||^2 - 1) / 187,5, com as "
               "normas nos pontos de integracao: o segundo termo e "
            << Fmt(normaContribG, 2, true) << " em G e " << Fmt(normaContribN, 5) << " em N (modos interpolados de Nystrom sem "
            << "norma 1)\n";
        res << "  ortonormalidade dos modos nos pontos de integracao, max |F^T W F - I|: G " << Fmt(gramG20, 2, true) << " (M = "
            << m20 << "), " << Fmt(gramGM, 2, true) << " (M = " << Mp << "); N " << Fmt(gramN20, 2, true) << " (M = " << m20
            << "), " << Fmt(gramNM, 2, true) << " (M = " << Mp << ")\n";
        const int64_t k100 = std::min<int64_t>(100, Mp);
        res << "  ||f_k^N|| nos pontos: " << Fmt(normaN.head(k100).minCoeff(), 4) << " a " << Fmt(normaN.head(k100).maxCoeff(), 4)
            << " (k = 1.." << k100 << ")";
        if (Mp > 100)
            res << ", " << Fmt(normaN.segment(100, Mp - 100).minCoeff(), 4) << " a "
                << Fmt(normaN.segment(100, Mp - 100).maxCoeff(), 4) << " (k = 101.." << Mp << ")";
        res << "; ||f_k^G||: " << Fmt(normaG.head(Mp).minCoeff(), 12) << " a " << Fmt(normaG.head(Mp).maxCoeff(), 12) << "\n";
    }
    // covariancias
    {
        res << "\nCovariancia truncada com os pontos de referencia, sum_{k<M} lambda_k f_k(x) f_k(x_r), contra o nucleo C(x, x_r), M = "
            << Mp << ": RMS na area / max |erro|\n";
        for (int r = 0; r < nref; r++) {
            res << "  (" << refs[r][0] << "; " << refs[r][1] << "):";
            const std::pair<const char *, const Eigen::VectorXd *> lin[] = {{"G", &covG[r]}, {"N", &covN[r]}, {"A", &covA[r]}};
            for (const auto &l : lin) {
                const Estat e = Estatistica(P.w, *l.second - covAlvo[r]);
                res << "  " << l.first << " " << Fmt(e.rms, 4, true) << " / " << Fmt(e.maxabs, 4, true);
            }
            res << "\n";
        }
    }
    // Monte Carlo
    if (R > 0) {
        const double sem = media * cov / std::sqrt((double)R);
        res << "\nMonte Carlo, R = " << R << " realizacoes (REAL = 0.." << R - 1 << ", os mesmos xi nos tres metodos), M = " << Mp
            << " (A2000: 2000): RMS na area e max |.| (kPa) dos erros da media e do desvio de c em cada ponto\n"
            << "  (ruido amostral esperado: erro padrao da media " << Fmt(sem, 3) << " kPa; do desvio cerca de "
            << Fmt(media * cov * std::sqrt((std::exp(4. * ln.sigma_ln * ln.sigma_ln) + 2. * std::exp(3. * ln.sigma_ln * ln.sigma_ln) +
                                           3. * std::exp(2. * ln.sigma_ln * ln.sigma_ln) - 4.) / (4. * R)), 3)
            << " kPa)\n"
            << "  campo   media - " << media << "     media - E_M[c]    desvio - " << media * cov << "     desvio - DP_M[c]\n";
        const char *nomes[] = {"G", "N", "A", "A2000"};
        const Eigen::VectorXd *am[] = {&amG, &amN, &amA, &amA2000}, *ad[] = {&adG, &adN, &adA, &adA2000};
        for (int m = 0; m < nmc; m++) {
            const Estat e1 = Estatistica(P.w, (mcMedia[m].array() - media).matrix()), e2 = Estatistica(P.w, mcMedia[m] - *am[m]);
            const Estat e3 = Estatistica(P.w, (mcDP[m].array() - media * cov).matrix()), e4 = Estatistica(P.w, mcDP[m] - *ad[m]);
            res << "  " << std::left << std::setw(7) << nomes[m] << std::right << std::setw(8) << Fmt(e1.rms, 4) << " / "
                << Fmt(e1.maxabs, 3) << std::setw(9) << Fmt(e2.rms, 4) << " / " << Fmt(e2.maxabs, 3) << std::setw(10)
                << Fmt(e3.rms, 4) << " / " << Fmt(e3.maxabs, 3) << std::setw(9) << Fmt(e4.rms, 4) << " / " << Fmt(e4.maxabs, 3)
                << "   (media na area: media " << Fmt(Estatistica(P.w, mcMedia[m]).media, 4) << ", desvio "
                << Fmt(Estatistica(P.w, mcDP[m]).media, 4) << ")\n";
        }
    }
    // varredura
    {
        res << "\nVarredura em M (realizacao REAL = " << real << "; RMS na area e max |.| das diferencas de coesao, kPa; E[dZ]: RMS na "
               "area da diferenca esperada dos campos Z; 1 - media de sigma_M^2)\n"
            << "      M  c_G-c_N RMS/max     c_G-c_A RMS/max     c_N-c_A RMS/max     E[dZ] G-N  G-A    N-A     "
               "1-s2 G  N       A       desvio c_G   1-sum lambda/187,5 G  N\n";
        for (const auto &L : varredura)
            res << "  " << std::setw(5) << L.M << std::setw(10) << Fmt(L.dGN.rms, 4) << "/" << std::left << std::setw(8)
                << Fmt(L.dGN.maxabs, 3) << std::right << std::setw(10) << Fmt(L.dGA.rms, 4) << "/" << std::left << std::setw(8)
                << Fmt(L.dGA.maxabs, 3) << std::right << std::setw(10) << Fmt(L.dNA.rms, 4) << "/" << std::left << std::setw(8)
                << Fmt(L.dNA.maxabs, 3) << std::right << std::setw(10) << Fmt(L.ezGN, 4) << std::setw(7) << Fmt(L.ezGA, 3)
                << std::setw(7) << Fmt(L.ezNA, 3) << std::setw(9) << Fmt(L.erroVarG, 4) << std::setw(8) << Fmt(L.erroVarN, 4)
                << std::setw(8) << Fmt(L.erroVarA, 4) << std::setw(10) << Fmt(L.cG.dp, 4) << std::setw(14)
                << Fmt(L.erroVarLamG, 4) << std::setw(8) << Fmt(L.erroVarLamN, 4) << "\n";
    }
    res << "\nTempos (s): pontos e modos de Galerkin " << Fmt(tPontosG, 2) << ", Nystrom nos pontos " << Fmt(tNysPontos, 2)
        << ", forma fechada " << Fmt(tFF, 2) << ", Monte Carlo " << Fmt(tMC, 2) << "; total ate aqui " << Fmt(Segundos(tInicio), 1)
        << "; pico de memoria " << Fmt(PicoMemoriaMB(), 0) << " MB\n";

    // ---- arquivos ----
    {
        std::ofstream f(saida + "/pontos.csv");
        f << "elemento,ponto,x,y,peso,zG,zN,zA,zA2000,cG,cN,cA,cA2000,dGN,dGA,dNA,dGA2000,s2G,s2N,s2A,s2A2000,"
             "ezGN,ezGA,ezNA";
        for (int r = 0; r < nref; r++) f << ",cov" << r + 1 << "_G,cov" << r + 1 << "_N,cov" << r + 1 << "_A,cov" << r + 1 << "_alvo";
        f << ",mcmG,mcmN,mcmA,mcmA2000,mcdG,mcdN,mcdA,mcdA2000,alvo_mG,alvo_mN,alvo_mA,alvo_mA2000,alvo_dG,alvo_dN,alvo_dA,"
             "alvo_dA2000,mc_rms_dGN,mc_rms_dGA,mc_rms_dNA";
        for (size_t i = 0; i < nMv; i++) f << ",zG" << nomeMv[i] << ",zN" << nomeMv[i] << ",dGN" << nomeMv[i];
        f << '\n';
        // x, y, peso e os campos Z com 17 algarismos (o double exato: conferencias de fora sem o limite da impressao);
        // o resto com 9
        const int pc = std::numeric_limits<double>::max_digits10, pn = 9;
        f << std::setprecision(pn);
        for (int64_t p = 0; p < np; p++) {
            f << P.geoel[p] << ',' << P.ip[p] << ',' << std::setprecision(pc) << P.x[p] << ',' << P.y[p] << ',' << P.w[p] << ','
              << zG(p) << ',' << zN(p) << ',' << zA(p) << ',' << zA2000(p) << std::setprecision(pn) << ',' << cG(p) << ','
              << cN(p) << ',' << cA(p) << ',' << cA2000(p) << ',' << dGN(p) << ',' << dGA(p) << ',' << dNA(p) << ','
              << dGA2000(p) << ',' << s2G(p) << ',' << s2N(p) << ',' << s2A(p) << ',' << s2A2000(p) << ',' << ezGN(p) << ','
              << ezGA(p) << ',' << ezNA(p);
            for (int r = 0; r < nref; r++)
                f << ',' << covG[r](p) << ',' << covN[r](p) << ',' << covA[r](p) << ',' << covAlvo[r](p);
            for (int m = 0; m < nmc; m++) f << ',' << mcMedia[m](p);
            for (int m = 0; m < nmc; m++) f << ',' << mcDP[m](p);
            f << ',' << amG(p) << ',' << amN(p) << ',' << amA(p) << ',' << amA2000(p) << ',' << adG(p) << ',' << adN(p) << ','
              << adA(p) << ',' << adA2000(p) << ',' << mcRmsD[0](p) << ',' << mcRmsD[1](p) << ',' << mcRmsD[2](p);
            for (size_t i = 0; i < nMv; i++)
                f << ',' << std::setprecision(pc) << zGv[i](p) << ',' << zNv[i](p) << std::setprecision(pn) << ',' << dGNv[i](p);
            f << '\n';
        }
    }
    {
        std::ofstream f(saida + "/modos.csv");
        f << "k,lambda_G,lambda_N,lambda_A,dif_rel_NG,acordo,melhor_j,cos_melhor,subespaco3,norma_N,sinal_G\n"
          << std::setprecision(12);
        for (int64_t k = 0; k < K; k++) {
            Eigen::Index jm = 0;
            Cos.row(k).cwiseAbs().maxCoeff(&jm);
            double sub = 0.;
            for (int64_t j = std::max<int64_t>(0, k - 3); j <= std::min<int64_t>(K - 1, k + 3); j++) sub += Cos(k, j) * Cos(k, j);
            f << k + 1 << ',' << kl.lambda(k) << ',' << nys.lambda(k) << ',' << lamA(k) << ','
              << (nys.lambda(k) - kl.lambda(k)) / kl.lambda(k) << ',' << acordo(k) << ',' << jm + 1 << ','
              << std::fabs(Cos(k, jm)) << ',' << std::sqrt(sub) << ',' << normaN(k) << ',' << integralSinal(k) << '\n';
        }
    }
    {
        std::ofstream f(saida + "/varredura.csv");
        f << "M,rms_dGN,max_dGN,rms_dGA,max_dGA,rms_dNA,max_dNA,corr_GN,corr_GA,ez_GN,ez_GA,ez_NA,erro_var_G,erro_var_N,"
             "erro_var_A,desvio_cG,erro_var_lambda_G,erro_var_lambda_N\n"
          << std::setprecision(10);
        for (const auto &L : varredura)
            f << L.M << ',' << L.dGN.rms << ',' << L.dGN.maxabs << ',' << L.dGA.rms << ',' << L.dGA.maxabs << ',' << L.dNA.rms
              << ',' << L.dNA.maxabs << ',' << L.corrGN << ',' << L.corrGA << ',' << L.ezGN << ',' << L.ezGA << ',' << L.ezNA
              << ',' << L.erroVarG << ',' << L.erroVarN << ',' << L.erroVarA << ',' << L.cG.dp << ',' << L.erroVarLamG << ','
              << L.erroVarLamN << '\n';
    }
    // VTK
    if (gravaVTK) {
        const auto tv0 = Relogio::now();
        const Eigen::VectorXd evG = (s2G.array() - 1.).matrix(), evN = (s2N.array() - 1.).matrix(), evA = (s2A.array() - 1.).matrix();
        const Eigen::VectorXd emG = (mcMedia[0].array() - media).matrix(), emN = (mcMedia[1].array() - media).matrix(),
                              emA = (mcMedia[2].array() - media).matrix();
        const Eigen::VectorXd edG = (mcDP[0].array() - media * cov).matrix(), edN = (mcDP[1].array() - media * cov).matrix(),
                              edA = (mcDP[2].array() - media * cov).matrix();
        std::vector<Eigen::VectorXd> ecG(nref), ecN(nref), ecA(nref);
        for (int r = 0; r < nref; r++) {
            ecG[r] = covG[r] - covAlvo[r];
            ecN[r] = covN[r] - covAlvo[r];
            ecA[r] = covA[r] - covAlvo[r];
        }
        std::vector<std::pair<std::string, const Eigen::VectorXd *>> campos = {
            {"cG", &cG}, {"cN", &cN}, {"cA", &cA}, {"cA2000", &cA2000},
            {"dGN", &dGN}, {"dGA", &dGA}, {"dNA", &dNA}, {"dGA2000", &dGA2000},
            {"erro_var_G", &evG}, {"erro_var_N", &evN}, {"erro_var_A", &evA},
            {"ezGN", &ezGN}, {"ezGA", &ezGA}, {"ezNA", &ezNA}};
        // c_G - c_N com M_VERIF termos, logo depois das diferencas da comparacao principal
        for (size_t i = 0; i < nMv; i++) campos.insert(campos.begin() + 8 + (std::ptrdiff_t)i, {"dGN" + nomeMv[i], &dGNv[i]});
        for (int r = 0; r < nref; r++) {
            campos.push_back({"erro_cov" + std::to_string(r + 1) + "_G", &ecG[r]});
            campos.push_back({"erro_cov" + std::to_string(r + 1) + "_N", &ecN[r]});
            campos.push_back({"erro_cov" + std::to_string(r + 1) + "_A", &ecA[r]});
        }
        if (R > 0) {
            campos.insert(campos.end(), {{"erro_mc_media_G", &emG}, {"erro_mc_media_N", &emN}, {"erro_mc_media_A", &emA},
                                         {"erro_mc_desvio_G", &edG}, {"erro_mc_desvio_N", &edN}, {"erro_mc_desvio_A", &edA},
                                         {"mc_rms_dGN", &mcRmsD[0]}, {"mc_rms_dGA", &mcRmsD[1]}, {"mc_rms_dNA", &mcRmsD[2]}});
        }
        GravaVTKPontos(saida + "/campos_pontos.vtk", P, campos);
        GravaVTKElementos(saida + "/campos_elementos.vtk", gmesh, P, campos);
        // um arquivo por campo principal, com o TPZVTKGeoMesh do NeoPZ (o valor de cada triangulo em ElData; os elementos
        // de contorno ficam de fora com -999)
        std::set<std::string> principais = {"cG", "cN", "cA", "cA2000", "dGN", "dGA", "dNA", "dGA2000", "erro_var_G",
                                            "erro_var_N", "erro_var_A", "erro_mc_media_G", "erro_mc_media_N",
                                            "erro_mc_media_A", "erro_mc_desvio_G", "erro_mc_desvio_N", "erro_mc_desvio_A"};
        for (size_t i = 0; i < nMv; i++) principais.insert("dGN" + nomeMv[i]);
        int narquivos = 0; // sem Monte Carlo (R = 0), os campos erro_mc_* nao existem
        for (const auto &c : campos) {
            if (!principais.count(c.first)) continue;
            const auto m = MediaPorElemento(gmesh, P, *c.second);
            TPZVec<REAL> elData(gmesh->NElements());
            for (int64_t e = 0; e < gmesh->NElements(); e++) elData[e] = m[e];
            std::ofstream f(saida + "/campos_elementos_" + c.first + ".vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, f, elData);
            narquivos++;
        }
        res << "VTK: " << campos.size() << " campos em campos_pontos.vtk e campos_elementos.vtk; " << narquivos
            << " arquivos campos_elementos_<campo>.vtk (" << Fmt(Segundos(tv0), 1) << " s)\n";
    }
    res << "Saidas em " << saida << ": pontos.csv, modos.csv, varredura.csv, resumo.txt" << (gravaVTK ? ", *.vtk" : "")
        << "; tempo total " << Fmt(Segundos(tInicio), 1) << " s\n";
    std::cout << res.str() << std::flush;
    std::ofstream(saida + "/resumo.txt") << resumo << res.str();
    return 0;
}
