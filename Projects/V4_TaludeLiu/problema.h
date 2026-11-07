#pragma once
// V4_TaludeLiu: talude nao saturado de 5 m (1:2) sob chuva permanente (q/k_s = 0.01) de Liu, Han, Huang, Huang e Cheng
// (2024), Comput. Geotech. 165:105938: perfis de poropressao nas secoes X-X e Y-Y (Fig. 3), FS (1,525 por Bishop) e
// Monte Carlo com campos lognormais de c', phi' e k_s sem rotacao (alfa = 0 e 90 graus), cenarios 1 a 3 (Tabela 3).
//
// Modelo (ver README.md): variante iii (esqueleto rigido, ar a pressao atmosferica: so a equacao de Richards), o regime
// permanente atingido por um transiente longo com chuva constante; depois o FS por reducao de resistencia (SRM, com as
// poropressoes congeladas e nova tentativa do F que falhou), Mohr-Coulomb associado, chi constante = tan(phi_b)/tan(phi')
// (phi_b = 2/3 phi'), peso especifico constante. Retencao de Fredlund-Xing do artigo (a = 5 kPa, n = 2, m = 1,
// theta_s = 0.4); a permeabilidade FXH (eq. 9) nao existe na biblioteca e foi trocada por k_r = S^p com p ajustado a
// FXH (scripts/ajuste_kr.py). Campos de Karhunen-Loeve com o nucleo exponencial separavel (KLSeparavel.h) no lugar do
// gaussiano do artigo, com a mesma escala de flutuacao: L = sqrt(pi) theta / 2.
//
// Uso: main.cpp (o contrato do Monte Carlo, uma realizacao por execucao) e testes.cpp. As funcoes do RahardjoDriver.h nao
// sao inline: este cabecalho so pode ser incluido num unico .cpp por executavel.
#include "RahardjoDriver.h" // UnsatParams, UnsatMaterial, ReadGmshMesh, CreateSlopeCompMesh, DeleteCompMesh, PostProcessVTK
#include "KLSeparavel.h"    // Karhunen-Loeve com covariancia exponencial separavel, lognormal e gerador normal

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

#ifndef VERIF_DIR
#define VERIF_DIR "."
#endif

namespace v4 {

// ------------------------------------------------------------------------------------------------------------------
// ambiente e utilidades
inline double EnvD(const char *k, double d) { const char *v = std::getenv(k); return v && *v ? std::atof(v) : d; }
inline int EnvI(const char *k, int d) { const char *v = std::getenv(k); return v && *v ? std::atoi(v) : d; }
inline std::string EnvS(const char *k, const std::string &d) { const char *v = std::getenv(k); return v && *v ? std::string(v) : d; }
inline std::string Num(double v, int p = 10) { std::ostringstream o; o << std::setprecision(p) << v; return o.str(); }
using Relogio = std::chrono::steady_clock;
inline double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }
/// nivel de mensagens na tela (0: so o resumo; 1: passos e FS; 2: cada solucao do SRM)
inline int &Verboso() { static int v = 1; return v; }

// geometria da Fig. 1 (m): crista H(10,10), pe G(20,5), meia altura da face I(15,7.5); secoes X-X (x = 10) e Y-Y (x = 15)
constexpr double kXCrista = 10., kZCrista = 10., kXPe = 20., kZPe = 5., kXSecaoX = 10., kXSecaoY = 15.;
/// cota do terreno em x (crista z = 10 ate x = 10, face 1:2 ate x = 20, pe z = 5)
inline double CotaTerreno(double x) {
    if (x <= kXCrista) return kZCrista;
    if (x >= kXPe) return kZPe;
    return kZCrista - 0.5 * (x - kXCrista);
}
/// modo de ruptura pela profundidade D_f da superficie de ruptura abaixo da crista (Fig. 9, Cheng et al. 2017):
/// talude (0 a 5 m), pe (5 a 6 m) e base (6 a 10 m); os limites ficam na classe mais rasa
inline std::string ModoRuptura(double Df) {
    if (!(Df >= 0.)) return "";
    return Df <= 5. ? "talude" : Df <= 6. ? "pe" : "base";
}

// ------------------------------------------------------------------------------------------------------------------
/// Parametros de uma execucao (lidos do ambiente por DoAmbiente; os valores padrao sao os do artigo e as escolhas do
/// README.md)
struct Config {
    // ---- caso: det (medias), s1 (so k_s aleatorio), s2 (c' e phi'), s3 (os tres); angulo do acamamento 0 ou 90 graus
    std::string caso = "det";
    int angulo = 0;
    int real = 0, campos = 1;
    uint64_t semente = 20231122; // padrao fixo (a data de publicacao online do artigo)
    // ---- solo (Tabela 1 e sec. 3.1): c' = 10 kPa, phi' = 15 graus, phi_b = 2/3 phi' = 10 graus, gamma = 20 kN/m3
    double c = 10., phiDeg = 15., phibDeg = 10., gamma = 20.;
    // E e nu nao sao dados (equilibrio limite): os de Griffiths e Fenton (2004) para esta familia de taludes
    double E = 1.e5, nu = 0.3;
    // ---- hidraulica (sec. 3.1): k_s = 1e-5 m/s, theta_s = 0.4, theta_r = 0, Fredlund-Xing a = 5 kPa, n = 2, m = 1
    double ks = 1.e-5, poros = 0.4, aFX = 5., nFX = 2., mFX = 1.;
    // k_r = S^p ajustado a FXH (eq. 9) com psi_aev = 2.58 kPa (tangente no ponto de inflexao): scripts/ajuste_kr.py
    double pKr = 5.57;
    // ---- chuva e contornos: q/k_s = 0.01 (q = 1e-7 m/s), carga total de 6 m em BC e de 1 m em ED
    double qks = 0.01, hEsq = 6., hDir = 1.;
    int chuvaVertical = 0; ///< 0: q normal a superficie por comprimento real (unit flux do SEEP/W); 1: por projecao horizontal
    // ---- transiente ate o regime permanente: dt cresce de DT0_H por FATOR_DT ate DT_MAX_H; para quando a maior
    // variacao de p_w num passo com dt >= DT_PERM_H fica abaixo de TOL_PERM (kPa); T_MAX_H limita o tempo simulado
    double dt0_h = 1., fatorDt = 2., dtMax_h = 24. * 365., dtPerm_h = 24. * 30., tolPerm = 1.e-3, tMax_h = 24. * 365. * 50.;
    int maxCortes = 8;
    double sucIni = -1.; ///< succao maxima da condicao inicial (kPa); <= 0: o patamar k_r(psi) = q/k_s do modelo
    std::string pwPermanente; ///< PW_PERMANENTE: le o regime permanente de um pw.csv (k_s constante: det e s2)
    // ---- campos (Tabela 1): COV 0.3 (c'), 0.2 (phi'), 0.6 (k_s); theta_1 = 30 m, theta_2 = 3 m; rho_c'phi' = -0.5 no
    // artigo, aqui 0 (independentes) por padrao
    double covC = 0.3, covPhi = 0.2, covK = 0.6, theta1 = 30., theta2 = 3., rho = 0.;
    double Lx = -1., Ly = -1.; ///< > 0 trocam os comprimentos do nucleo exponencial (m)
    int klM = 2000;
    double klFrac = 1.;
    // ---- fator de seguranca (SRM): estado geostatico com c'/F0, tan(phi')/F0; F pela lista 0.7, 0.8, ... e bissecao
    double F0 = 0.6, fsTol = 0.01;
    double tolU = 1.e-7, tolP = 1.e-6;
    int maxIt = 40, maxItFrozen = 25;
    // massa que desliza: pontos com |du| >= limiar * max|du| (du: incremento de deslocamento entre os dois ultimos
    // estados convergidos do SRM)
    double massaLimiar = 0.5; // o padrao do SlopeRandomField: o meio da banda de cisalhamento
    // ---- malha e saidas
    std::string malha;        ///< padrao: malhas/talude_liu_h05.msh (0,5 m, o tamanho do artigo)
    int pwCsv = 1, camposCsv = 0, mecCsv = 0, vtk = 0;

    /// comprimentos do nucleo exponencial com a escala de flutuacao do gaussiano: delta = sqrt(pi) theta = 2 L; com
    /// alfa = 90 graus os eixos maior e menor trocam (o eixo maior na vertical)
    void Comprimentos(double &lx, double &ly) const {
        const double l1 = std::sqrt(M_PI) * theta1 / 2., l2 = std::sqrt(M_PI) * theta2 / 2.;
        lx = angulo == 90 ? l2 : l1;
        ly = angulo == 90 ? l1 : l2;
        if (Lx > 0.) lx = Lx;
        if (Ly > 0.) ly = Ly;
    }
    bool CampoK() const { return campos && (caso == "s1" || caso == "s3"); }
    bool CampoResist() const { return campos && (caso == "s2" || caso == "s3"); }
    double Chi() const { return std::tan(phibDeg * M_PI / 180.) / std::tan(phiDeg * M_PI / 180.); }

    static Config DoAmbiente() {
        Config cf;
        cf.caso = EnvS("CASO", cf.caso);
        if (cf.caso != "det" && cf.caso != "s1" && cf.caso != "s2" && cf.caso != "s3")
            throw std::invalid_argument("CASO deve ser det, s1, s2 ou s3");
        cf.angulo = EnvI("ANGULO", cf.angulo);
        if (cf.angulo != 0 && cf.angulo != 90) throw std::invalid_argument("ANGULO deve ser 0 ou 90 (sem rotacao)");
        cf.real = EnvI("REAL", cf.real);
        cf.campos = cf.caso == "det" ? 0 : EnvI("CAMPOS", cf.campos);
        cf.semente = (uint64_t)EnvD("SEMENTE", (double)cf.semente);
        cf.c = EnvD("COES", cf.c);
        cf.phiDeg = EnvD("PHI_DEG", cf.phiDeg);
        cf.phibDeg = EnvD("PHIB_DEG", cf.phiDeg * 2. / 3.);
        cf.gamma = EnvD("GAMMA", cf.gamma);
        cf.E = EnvD("E", cf.E);
        cf.nu = EnvD("NU", cf.nu);
        cf.ks = EnvD("KS", cf.ks);
        cf.poros = EnvD("POROS", cf.poros);
        cf.aFX = EnvD("A_FX", cf.aFX);
        cf.nFX = EnvD("N_FX", cf.nFX);
        cf.mFX = EnvD("M_FX", cf.mFX);
        cf.pKr = EnvD("P_KR", cf.pKr);
        cf.qks = EnvD("Q_KS", cf.qks);
        cf.hEsq = EnvD("H_ESQ", cf.hEsq);
        cf.hDir = EnvD("H_DIR", cf.hDir);
        cf.chuvaVertical = EnvI("CHUVA_VERTICAL", cf.chuvaVertical);
        cf.dt0_h = EnvD("DT0_H", cf.dt0_h);
        cf.fatorDt = EnvD("FATOR_DT", cf.fatorDt);
        cf.dtMax_h = EnvD("DT_MAX_H", cf.dtMax_h);
        cf.dtPerm_h = EnvD("DT_PERM_H", cf.dtPerm_h);
        cf.tolPerm = EnvD("TOL_PERM", cf.tolPerm);
        cf.tMax_h = EnvD("T_MAX_H", cf.tMax_h);
        cf.maxCortes = EnvI("MAX_CORTES", cf.maxCortes);
        cf.sucIni = EnvD("SUC_INI", cf.sucIni);
        cf.pwPermanente = EnvS("PW_PERMANENTE", "");
        cf.covC = EnvD("COV_C", cf.covC);
        cf.covPhi = EnvD("COV_PHI", cf.covPhi);
        cf.covK = EnvD("COV_K", cf.covK);
        cf.theta1 = EnvD("THETA1", cf.theta1);
        cf.theta2 = EnvD("THETA2", cf.theta2);
        cf.rho = EnvD("RHO", cf.rho);
        if (!(cf.rho > -1. && cf.rho < 1.)) throw std::invalid_argument("RHO deve estar em (-1, 1)");
        cf.Lx = EnvD("LX", cf.Lx);
        cf.Ly = EnvD("LY", cf.Ly);
        cf.klM = EnvI("KL_M", cf.klM);
        cf.klFrac = EnvD("KL_FRAC", cf.klFrac);
        cf.F0 = EnvD("FS_F0", cf.F0);
        cf.fsTol = EnvD("FS_TOL", cf.fsTol);
        cf.tolU = EnvD("TOLU", cf.tolU);
        cf.tolP = EnvD("TOLP", cf.tolP);
        cf.maxIt = EnvI("MAXIT", cf.maxIt);
        cf.maxItFrozen = EnvI("MAXIT_FROZEN", cf.maxItFrozen);
        cf.massaLimiar = EnvD("MASSA_LIMIAR", cf.massaLimiar);
        cf.malha = EnvS("MALHA", cf.malha);
        cf.pwCsv = EnvI("PW_CSV", cf.pwCsv);
        cf.camposCsv = EnvI("CAMPOS_CSV", cf.camposCsv);
        cf.mecCsv = EnvI("MEC_CSV", cf.mecCsv);
        cf.vtk = EnvI("VTK", cf.vtk);
        return cf;
    }
};

/// caminho da malha: MALHA (se dada) ou malhas/talude_liu_h05.msh da pasta do projeto; um nome sem '/' e procurado
/// na pasta malhas/ do projeto
inline std::string CaminhoMalha(const std::string &m) {
    if (m.empty()) return std::string(VERIF_DIR) + "/malhas/talude_liu_h05.msh";
    if (m.find('/') == std::string::npos && !std::filesystem::exists(m)) return std::string(VERIF_DIR) + "/malhas/" + m;
    return m;
}

/// Parametros do material (UnsatParams do RahardjoDriver.h)
inline UnsatParams Parametros(const Config &cf) {
    UnsatParams p;
    p.E = cf.E, p.nu = cf.nu, p.c = cf.c, p.phi = cf.phiDeg * M_PI / 180., p.gamma = cf.gamma;
    p.mixture = false; // peso especifico constante (sec. 3.1)
    p.n = cf.poros, p.ks = cf.ks;
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
    p.hydro.fA = cf.aFX, p.hydro.fN = cf.nFX, p.hydro.fM = cf.mFX;
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
    p.hydro.fPower = cf.pKr;
    p.hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // sem efeito: p_a = 0 imposto (variante iii)
    p.bishop = TUnsat::EChiConstant, p.chi = cf.Chi();     // tan(phi_b) = chi tan(phi'); chi = 1 onde p_w > 0
    p.coupled = false; // variante iii: esqueleto rigido no transiente
    return p;
}

/// succao do patamar de infiltracao permanente com gradiente unitario: k_r(psi) = q/k_s (bissecao em log psi)
inline double SuccaoPatamar(const TPZUnsaturatedHydraulics &h, double qks) {
    double lo = 1.e-4, hi = 1.e4;
    for (int it = 0; it < 200; ++it) {
        const double m = std::sqrt(lo * hi);
        const auto st = h.Evaluate(-m, 0.);
        (st.krw > qks ? lo : hi) = m;
    }
    return std::sqrt(lo * hi);
}

// ------------------------------------------------------------------------------------------------------------------
/// uma solucao da busca do SRM: um dos 4 passos de carga da gravidade com c'/F0, tan(phi')/F0 ("gravidade") ou um F
/// tentado depois ("srm")
struct TentativaSRM {
    const char *fase;
    REAL g, F;
    bool conv, retry;
    int iter;
    double seg;
};

/// FactorSRM do RahardjoDriver.h com a nova tentativa do F que falhou (como no SlopeRandomField do branch
/// UnsaturatedSoil): parte do estado sem tensoes com as pressoes de x congeladas, estado geostatico com F0 em 4
/// incrementos, F pela lista Fs e bissecao ate tol; o F que falhou e tentado de novo a partir do ultimo convergido (o
/// incremento agora e < tol). Fparada: o F em que o FactorSRM pararia sem a nova tentativa; Xant: o estado aceito
/// antes do ultimo (o ultimo, com F = FS, fica em an.X()). Retorna -1 se o estado geostatico nao convergir (FS < F0).
inline REAL FactorSRMRetry(TPZUnsatAnalysis &an, const TPZFMatrix<STATE> &x, REAL F0, REAL tol,
                           const std::vector<REAL> &Fs, REAL &Fparada, std::vector<TentativaSRM> &tent,
                           TPZFMatrix<STATE> &Xant) {
    an.ResetMemory(); // mantem as propriedades de cada ponto (os campos)
    TPZFMatrix<STATE> &X = an.X();
    X = x;
    for (int64_t i = 0; i < an.NEq(); i++)
        if (!an.IsPressure(i)) X(i, 0) = 0.;
    an.Restore();
    const auto fixed = an.FrozenFixed(X);
    an.SetStrengthReductionFactor(F0);
    an.SetGravityScalesSuction(true);
    Fparada = -1.;
    Xant = X;
    const int n = 4;
    for (int k = 1; k <= n; k++) {
        an.SetGravityFactor(REAL(k) / n);
        const auto t0 = Relogio::now();
        const auto r = an.Solve(0., fixed, {}, true);
        tent.push_back({"gravidade", REAL(k) / n, F0, r.conv, false, r.iter, Segundos(t0)});
        if (!r.conv) return -1.;
        Xant = an.X();
        an.Accept();
    }
    REAL Fok = F0, Ffail = -1.;
    size_t next = 0;
    while (true) {
        REAL F;
        bool retry = false;
        if (Ffail > 0.) {
            if (Ffail - Fok < tol) {
                if (Fparada < 0.) Fparada = Fok;
                F = Ffail;
                retry = true;
            } else F = 0.5 * (Fok + Ffail);
        } else {
            while (next < Fs.size() && Fs[next] <= Fok) next++;
            if (next == Fs.size()) break;
            F = Fs[next];
        }
        an.SetStrengthReductionFactor(F);
        const auto t0 = Relogio::now();
        const auto r = an.Solve(0., fixed, {}, true);
        tent.push_back({"srm", 1., F, r.conv, retry, r.iter, Segundos(t0)});
        if (Verboso() >= 2)
            std::cout << "  SRM: F = " << F << (r.conv ? " convergiu" : " nao convergiu") << " (" << r.iter << " it.)"
                      << std::endl;
        if (r.conv) {
            Xant = an.X();
            an.Accept();
            Fok = F;
            if (retry) Ffail = -1.;
        } else {
            an.Restore();
            if (retry) break;
            Ffail = F;
        }
    }
    if (Fparada < 0.) Fparada = Fok;
    an.SetStrengthReductionFactor(1.);
    return Fok;
}

// ------------------------------------------------------------------------------------------------------------------
/// regime permanente: resultado do transiente
struct InfoPermanente {
    bool ok = false, lido = false;
    double t_h = 0., dpw = 0., psiPatamar = 0., sucIni = 0.;
    int passos = 0, cortes = 0, newton = 0, faceNos = 0;
    // balanco de agua do ultimo passo (m3 por metro, divididos por dt: m2/s): chuva que entrou, laterais (entrada > 0),
    // face de percolacao, variacao do armazenamento e o residuo
    double chuva = 0., esq = 0., dir = 0., face = 0., dV = 0., residuo = 0.;
    double pwMin = 0., pwMax = 0., seg = 0.;
};

/// fator de seguranca e mecanismo
struct InfoFS {
    double FS = -1., Fparada = -1., seg = 0.;
    int tentativas = 0, newton = 0;
    double uMax = -1., duMax = 0.;
    double massaArea = 0., massaXc = 0., massaYc = 0., massaYmin = 0., massaXmin = 0., massaXmax = 0., Df = -1.;
    double massaC = 0., massaPhi = 0., massaK = 0.;
    int massaPts = 0;
    std::vector<std::pair<double, double>> areasLimiar; ///< (limiar, area) para a sensibilidade
    std::string modo;
};

/// fator de seguranca por SRM com arc-length (FactorSRMArc): o FS e o resumo do caminho de equilibrio
struct InfoFSArco {
    double FS = -1., seg = 0.;
    SaidaArco saida; ///< incrementos, cortes de l, iteracoes de Newton e por que o caminho terminou
};

/// O talude: malha, analise e as rotinas do caso (contornos, campos, permanente, secoes e FS)
class Talude {
public:
    Config cf;
    UnsatParams p;
    TPZGeoMesh *gmesh = nullptr;
    TPZMultiphysicsCompMesh *cmesh = nullptr;
    std::unique_ptr<TPZUnsatAnalysis> an;
    TUnsat *solo = nullptr;
    std::map<int64_t, REAL> borda;           ///< nos da superficie (AHIGF) e comprimento tributario
    std::map<int64_t, REAL> contorno;        ///< equacoes impostas no transiente (u, p_a e p_w das laterais)
    std::set<int64_t> eqEsq, eqDir;          ///< equacoes de p_w impostas em BC e ED
    std::set<int64_t> saida;                 ///< face de percolacao (nos da superficie com p_w = 0)
    // pontos de integracao do solo
    std::vector<int64_t> ips, ipel;
    std::vector<double> xs, ys, area;
    double areaSolo = 0.;
    int klMUsado = 0;
    double klFracUsada = 0., LxUsado = 0., LyUsado = 0.;
    std::vector<TentativaSRM> tentativas;

    explicit Talude(const Config &c) : cf(c) {
        p = Parametros(cf);
        gmesh = ReadGmshMesh(CaminhoMalha(cf.malha), {{"Solo", matSolo}},
                             {{"Base", bcBase}, {"Direita", bcDireita}, {"Superficie", bcSuperficie}, {"Esquerda", bcEsquerda}});
        cmesh = CreateSlopeCompMesh(gmesh, 2, p);
        an = std::make_unique<TPZUnsatAnalysis>(cmesh, std::map<int, std::array<bool, 2>>{
            {bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
        an->SetTolerances(cf.tolU, cf.tolP, cf.maxIt, cf.maxItFrozen);
        solo = dynamic_cast<TUnsat *>(cmesh->FindMaterial(matSolo));
        if (!solo) throw std::runtime_error("material do solo nao encontrado");
        borda = an->BoundaryNodes(bcSuperficie, cf.chuvaVertical != 0);
        if (borda.empty()) throw std::runtime_error("superficie sem nos");
        // contornos da percolacao (variante iii): todos os deslocamentos fixos, p_a = 0 em todos os nos, carga total
        // hEsq em BC (x = 0, z <= hEsq) e hDir em ED (x = 30, z <= hDir); AB, FE e a base sao impermeaveis (fluxo
        // natural nulo)
        contorno = an->AllDisplacements(an->X());
        for (const int id : {bcEsquerda, bcDireita}) {
            const REAL h = id == bcEsquerda ? cf.hEsq : cf.hDir;
            for (const auto &b : an->BoundaryNodes(id)) {
                const REAL y = an->Coord(b.first, 1);
                if (y <= h + 1.e-9) {
                    const int64_t eq = an->WaterEq().at(b.first);
                    contorno[eq] = p.gw * (h - y);
                    (id == bcEsquerda ? eqEsq : eqDir).insert(eq);
                }
            }
        }
        for (const auto &e : an->AirEq()) contorno[e.second] = 0.;
        ColetaPontos();
    }
    ~Talude() {
        an.reset();
        if (cmesh) DeleteCompMesh(cmesh);
        delete gmesh;
    }
    Talude(const Talude &) = delete;
    Talude &operator=(const Talude &) = delete;

    /// propriedades do ponto i (as do material se o ponto nao tem campo); phi em graus
    double CPonto(size_t i) const {
        const auto &mp = (*solo->GetMemory())[ips[i]].m_elastoplastic_state.fmatprop;
        return mp.size() >= 2 ? mp[0] : p.c;
    }
    double PhiPonto(size_t i) const {
        const auto &mp = (*solo->GetMemory())[ips[i]].m_elastoplastic_state.fmatprop;
        return (mp.size() >= 2 ? mp[1] : p.phi) * 180. / M_PI;
    }
    double KPonto(size_t i) const { return (*solo->GetMemory())[ips[i]].fKMult; }

    /// Campos de Karhunen-Loeve no retangulo que contem a malha, nos pontos de integracao: tres sequencias normais
    /// independentes por realizacao (c', phi' e k_s, sementes derivadas de SEMENTE, REAL e o campo: o campo de k_s da
    /// realizacao n e o mesmo nos cenarios 1 e 3); phi' = rho Z_c + sqrt(1 - rho^2) Z_2 (correlacao cruzada das normais,
    /// como as eqs. 3-5 do artigo); lognormais com as medias e COV da Tabela 1
    void AplicaCampos() {
        if (!cf.CampoK() && !cf.CampoResist()) return;
        REAL x0 = 1.e30, x1 = -1.e30, y0 = 1.e30, y1 = -1.e30;
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) {
            TPZManVector<REAL, 3> c(3);
            gmesh->NodeVec()[i].GetCoordinates(c);
            x0 = std::min(x0, c[0]); x1 = std::max(x1, c[0]); y0 = std::min(y0, c[1]); y1 = std::max(y1, c[1]);
        }
        cf.Comprimentos(LxUsado, LyUsado);
        const klsep::KL2D kl(x0, y0, x1 - x0, y1 - y0, LxUsado, LyUsado, cf.klFrac, cf.klM);
        const auto P = kl.PreparaPontos(xs, ys);
        auto xi = [&](int campo) {
            return klsep::Normal(cf.semente * 0x9E3779B97F4A7C15ULL + 3ULL * (uint64_t)cf.real + campo).Vetor(kl.M());
        };
        auto mem = solo->GetMemory();
        if (cf.CampoResist()) {
            const auto zc = kl.Campo(P, xi(0)), z2 = kl.Campo(P, xi(1));
            const klsep::Lognormal lc(cf.c, cf.covC), lphi(cf.phiDeg, cf.covPhi);
            const double r2 = std::sqrt(1. - cf.rho * cf.rho);
            for (size_t i = 0; i < ips.size(); ++i) {
                const REAL c = lc(zc[i]), phi = lphi(cf.rho * zc[i] + r2 * z2[i]) * M_PI / 180.;
                (*mem)[ips[i]].m_elastoplastic_state.fmatprop = {c, phi, phi}; // {c, phi, psi}, psi = phi (associado)
            }
        }
        if (cf.CampoK()) {
            const auto zk = kl.Campo(P, xi(2));
            const klsep::Lognormal lk(1., cf.covK);
            for (size_t i = 0; i < ips.size(); ++i) (*mem)[ips[i]].fKMult = lk(zk[i]);
        }
        klMUsado = kl.M();
        klFracUsada = kl.fracao;
    }

    /// condicao inicial: hidrostatica com o lencol reto de (0, hEsq) a (30, hDir), succao limitada a sucIni
    double CondicaoInicial() {
        double suc = cf.sucIni > 0. ? cf.sucIni : SuccaoPatamar(p.hydro, cf.qks);
        TPZFMatrix<STATE> x(cmesh->NEquations(), 1, 0.);
        const REAL xr = 30.;
        for (const auto &eq : an->WaterEq()) {
            const REAL xx = an->Coord(eq.first, 0), y = an->Coord(eq.first, 1);
            const REAL hf = cf.hEsq + (cf.hDir - cf.hEsq) * xx / xr;
            x(eq.second, 0) = std::max(p.gw * (hf - y), -suc);
        }
        for (const auto &c : contorno) x(c.first, 0) = c.second;
        an->X() = x;
        an->Restore();
        saida.clear();
        return suc;
    }

    /// Transiente com chuva constante q = (q/k_s) k_s ate o regime permanente (ver Config). passos.csv: um passo
    /// aceito por linha, com o balanco de agua (volumes por metro no passo)
    InfoPermanente Permanente(std::ostream *fpas = nullptr) {
        const auto t0 = Relogio::now();
        InfoPermanente info;
        info.psiPatamar = SuccaoPatamar(p.hydro, cf.qks);
        info.sucIni = CondicaoInicial();
        an->SetStrengthReductionFactor(1.);
        const REAL hora = 3600., q = cf.qks * cf.ks, tolP = 1.e-6, tolQ = 1.e-10;
        if (fpas)
            *fpas << "passo,t_h,dt_h,chuva_m3,esq_m3,dir_m3,face_m3,dV_m3,residuo_m3,face_nos,its_face,newton,dpw_max_kPa,"
                     "pw_min_kPa,pw_max_kPa,seg\n"
                  << std::setprecision(10);
        // um passo de dt com a face de percolacao: o no da superficie com p_w > 0 fica com p_w = 0 enquanto a entrada
        // exigida for menor que a chuva; em caso de falha, volta ao estado do inicio do passo
        auto passo = [&](REAL dt, double &dpw) -> bool {
            const auto r0 = Relogio::now();
            auto tentativa = saida;
            for (const auto &b : borda)
                if (an->X().GetVal(an->WaterEq().at(b.first), 0) > tolP) tentativa.insert(b.first);
            int newtonPasso = 0;
            for (int it = 0; it < 40; it++) {
                auto fixed = contorno;
                std::map<int64_t, REAL> fq;
                for (const auto &b : borda) {
                    const int64_t eq = an->WaterEq().at(b.first);
                    if (tentativa.count(b.first)) fixed[eq] = 0.;
                    else fq[eq] = q * b.second;
                }
                const auto sr = an->Solve(dt, fixed, fq, false);
                newtonPasso += sr.iter;
                info.newton += sr.iter;
                if (!sr.conv) break;
                const auto R = an->FullResidual();
                const TPZFMatrix<STATE> sol = cmesh->Solution();
                auto nova = tentativa;
                for (const auto &b : borda) {
                    const int64_t eq = an->WaterEq().at(b.first);
                    if (tentativa.count(b.first)) {
                        if (R.GetVal(eq, 0) / (dt * b.second) > q + tolQ) nova.erase(b.first);
                    } else if (sol.GetVal(eq, 0) > tolP) nova.insert(b.first);
                }
                if (nova != tentativa) {
                    tentativa = nova;
                    continue;
                }
                REAL somaR = 0., chuva = 0., face = 0., esq = 0., dir = 0., pwmin = 1.e300, pwmax = -1.e300;
                dpw = 0.;
                for (const auto &e : an->WaterEq()) {
                    somaR += R.GetVal(e.second, 0);
                    const REAL pw = sol.GetVal(e.second, 0);
                    pwmin = std::min(pwmin, pw);
                    pwmax = std::max(pwmax, pw);
                    dpw = std::max(dpw, std::fabs(pw - an->X().GetVal(e.second, 0)));
                }
                for (const auto &b : borda) {
                    if (tentativa.count(b.first)) face += R.GetVal(an->WaterEq().at(b.first), 0);
                    else chuva += q * b.second * dt;
                }
                for (auto e : eqEsq) esq += R.GetVal(e, 0);
                for (auto e : eqDir) dir += R.GetVal(e, 0);
                an->Accept();
                saida = tentativa;
                info.passos++;
                info.chuva = chuva / dt, info.esq = esq / dt, info.dir = dir / dt, info.face = face / dt;
                info.dV = somaR / dt;
                info.residuo = (somaR - chuva - face - esq - dir) / dt;
                info.pwMin = pwmin, info.pwMax = pwmax, info.faceNos = (int)saida.size();
                if (fpas)
                    *fpas << info.passos << ',' << (info.t_h + dt / hora) << ',' << dt / hora << ',' << chuva << ',' << esq
                          << ',' << dir << ',' << face << ',' << somaR << ',' << somaR - chuva - face - esq - dir << ','
                          << saida.size() << ',' << it + 1 << ',' << newtonPasso << ',' << dpw << ',' << pwmin << ','
                          << pwmax << ',' << Segundos(r0) << '\n';
                return true;
            }
            an->Restore();
            return false;
        };
        REAL dt = cf.dt0_h * hora;
        int cortesSeguidos = 0;
        while (info.t_h < cf.tMax_h) {
            double dpw = 0.;
            if (!passo(dt, dpw)) {
                info.cortes++;
                if (++cortesSeguidos > cf.maxCortes) break;
                dt *= 0.5;
                continue;
            }
            cortesSeguidos = 0;
            info.t_h += dt / hora;
            info.dpw = dpw;
            if (Verboso() >= 2)
                std::cout << "  passo " << info.passos << ": t = " << info.t_h << " h, dt = " << dt / hora
                          << " h, max|dp_w| = " << dpw << " kPa" << std::endl;
            if (dt >= cf.dtPerm_h * hora * (1. - 1.e-12) && dpw < cf.tolPerm) {
                info.ok = true;
                break;
            }
            dt = std::min(dt * cf.fatorDt, cf.dtMax_h * hora);
        }
        info.seg = Segundos(t0);
        if (fpas) fpas->flush();
        return info;
    }

    /// grava p_w e p_a de todos os nos (pw.csv) do estado convergido
    void GravaPw(const std::string &arq) {
        std::ofstream f(arq);
        f << "no,x,y,pw_kPa,pa_kPa,S,superficie,face\n" << std::setprecision(12);
        for (const auto &e : an->WaterEq()) {
            const REAL pw = an->X().GetVal(e.second, 0), pa = an->X().GetVal(an->AirEq().at(e.first), 0);
            f << e.first << ',' << an->Coord(e.first, 0) << ',' << an->Coord(e.first, 1) << ',' << pw << ',' << pa << ','
              << solo->Hydraulics().Evaluate(pw, pa).S << ',' << borda.count(e.first) << ',' << saida.count(e.first)
              << '\n';
        }
    }
    /// le o regime permanente de um pw.csv da mesma malha (k_s constante: o caso deterministico e o cenario 2)
    bool LePw(const std::string &arq) {
        std::ifstream f(arq);
        if (!f) return false;
        std::string l;
        std::getline(f, l);
        TPZFMatrix<STATE> x(an->X());
        size_t n = 0;
        while (std::getline(f, l)) {
            std::stringstream ss(l);
            std::string v[8];
            for (int k = 0; k < 8 && std::getline(ss, v[k], ','); ++k) {}
            const int64_t no = std::atoll(v[0].c_str());
            const auto it = an->WaterEq().find(no);
            if (it == an->WaterEq().end()) return false;
            if (std::fabs(std::atof(v[1].c_str()) - an->Coord(no, 0)) > 1.e-6 ||
                std::fabs(std::atof(v[2].c_str()) - an->Coord(no, 1)) > 1.e-6) return false; // outra malha
            x(it->second, 0) = std::atof(v[3].c_str());
            x(an->AirEq().at(no), 0) = std::atof(v[4].c_str());
            if (std::atoi(v[7].c_str())) saida.insert(no);
            n++;
        }
        if (n != an->WaterEq().size()) return false;
        an->X() = x; // os deslocamentos ficam nulos (esqueleto rigido)
        an->Restore();
        return true;
    }

    /// p_w (P1) interpolado no ponto (x, z) da malha geometrica; NaN fora do dominio
    double PwEm(double x, double z) const {
        for (int64_t i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *g = gmesh->Element(i);
            if (!g || g->Dimension() != 2 || g->MaterialId() != matSolo) continue;
            double X[3], Y[3], P[3];
            for (int k = 0; k < 3; ++k) {
                X[k] = an->Coord(g->NodeIndex(k), 0);
                Y[k] = an->Coord(g->NodeIndex(k), 1);
                P[k] = an->X().GetVal(an->WaterEq().at(g->NodeIndex(k)), 0);
            }
            const double det = (Y[1] - Y[2]) * (X[0] - X[2]) + (X[2] - X[1]) * (Y[0] - Y[2]);
            const double l0 = ((Y[1] - Y[2]) * (x - X[2]) + (X[2] - X[1]) * (z - Y[2])) / det;
            const double l1 = ((Y[2] - Y[0]) * (x - X[2]) + (X[0] - X[2]) * (z - Y[2])) / det;
            const double l2 = 1. - l0 - l1, e = -1.e-9;
            if (l0 >= e && l1 >= e && l2 >= e) return l0 * P[0] + l1 * P[1] + l2 * P[2];
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    /// Fator de seguranca (SRM) com as poropressoes do estado convergido congeladas, e o mecanismo: a massa que desliza
    /// (pontos com |du| >= MASSA_LIMIAR max|du|), a sua area, o ponto mais baixo e D_f = 10 - z_min (Fig. 9). Volta ao
    /// estado de antes no fim (solucao e memoria).
    InfoFS FS(std::ostream *fmec = nullptr) {
        const auto t0 = Relogio::now();
        InfoFS r;
        const auto estado = an->X();
        const auto memoria = an->SaveMemory();
        std::vector<REAL> Fs;
        for (REAL f = 0.7; f < 3.05; f += 0.1) Fs.push_back(std::round(f * 10.) / 10.);
        Fs.push_back(3.5);
        Fs.push_back(4.);
        TPZFMatrix<STATE> Xant;
        tentativas.clear();
        r.FS = FactorSRMRetry(*an, estado, cf.F0, cf.fsTol, Fs, r.Fparada, tentativas, Xant);
        r.tentativas = (int)tentativas.size();
        for (const auto &t : tentativas) r.newton += t.iter;
        if (r.FS > 0.) {
            // |du| nos vertices e, em cada ponto, a media dos 3 vertices do seu triangulo
            std::map<int64_t, double> du;
            r.uMax = 0.;
            for (const auto &e : an->DispEq()) {
                du[e.first] = std::hypot(an->NodalDisplacement(an->X(), e.first, 0) - an->NodalDisplacement(Xant, e.first, 0),
                                         an->NodalDisplacement(an->X(), e.first, 1) - an->NodalDisplacement(Xant, e.first, 1));
                r.uMax = std::max(r.uMax, (double)std::hypot(an->NodalDisplacement(an->X(), e.first, 0),
                                                             an->NodalDisplacement(an->X(), e.first, 1)));
            }
            std::vector<double> dp(ips.size(), 0.);
            for (size_t i = 0; i < ips.size(); ++i) {
                TPZGeoEl *g = gmesh->Element(ipel[i]);
                double s = 0.;
                for (int k = 0; k < 3; ++k) {
                    const auto it = du.find(g->NodeIndex(k));
                    if (it != du.end()) s += it->second;
                }
                dp[i] = s / 3.;
                r.duMax = std::max(r.duMax, dp[i]);
            }
            if (r.duMax > 0.) {
                for (double lim : {0.1, 0.2, 0.3, 0.5}) {
                    double a = 0.;
                    for (size_t i = 0; i < ips.size(); ++i) if (dp[i] >= lim * r.duMax) a += area[i];
                    r.areasLimiar.push_back({lim, a});
                }
                r.massaYmin = r.massaXmin = 1.e300;
                r.massaXmax = -1.e300;
                for (size_t i = 0; i < ips.size(); ++i) {
                    if (dp[i] < cf.massaLimiar * r.duMax) continue;
                    const double a = area[i];
                    r.massaPts++;
                    r.massaArea += a;
                    r.massaXc += a * xs[i]; r.massaYc += a * ys[i];
                    r.massaC += a * CPonto(i); r.massaPhi += a * PhiPonto(i); r.massaK += a * KPonto(i);
                    r.massaYmin = std::min(r.massaYmin, ys[i]);
                    r.massaXmin = std::min(r.massaXmin, xs[i]);
                    r.massaXmax = std::max(r.massaXmax, xs[i]);
                    if (fmec) *fmec << ips[i] << ',' << xs[i] << ',' << ys[i] << ',' << area[i] << ',' << dp[i] << '\n';
                }
                r.massaXc /= r.massaArea; r.massaYc /= r.massaArea;
                r.massaC /= r.massaArea; r.massaPhi /= r.massaArea; r.massaK /= r.massaArea;
                r.Df = kZCrista - r.massaYmin;
                r.modo = ModoRuptura(r.Df);
            }
            if (cf.vtk) {
                an->SetStrengthReductionFactor(r.FS);
                PostProcessVTK(*an, "ruptura", 0);
            }
        }
        an->SetStrengthReductionFactor(1.);
        an->SetGravityFactor(1.);
        an->SetGravityScalesSuction(true);
        an->RestoreMemory(memoria);
        an->X() = estado;
        an->Restore();
        solo->SetFrozenPressures(false);
        r.seg = Segundos(t0);
        return r;
    }

    /// FS por SRM com arc-length (FactorSRMArc de Comum/RahardjoDriver.h, o codigo do LearningMohrCoulomb2 do branch
    /// UnsaturatedSoil, commit 2f92365, com os seus parametros padrao, ParamsSRM) e as poropressoes do estado convergido
    /// congeladas, como em FS(), que usa a bissecao. O FS e o maior F do caminho de equilibrio. Sem o mecanismo. Volta ao
    /// estado de antes no fim (solucao e memoria). tol > 0 troca a tolerancia da parada (padrao 0,001); criterio >= 0 troca
    /// o criterio de parada (ParamsArco::criterio; padrao 0, dF < tol; 1, (dF / l) s < tol).
    InfoFSArco FSArco(double tol = 0., int criterio = -1) {
        const auto t0 = Relogio::now();
        InfoFSArco r;
        const auto estado = an->X();
        const auto memoria = an->SaveMemory();
        ParamsArco par = ParamsSRM();
        if (tol > 0.) par.tol = tol;
        if (criterio >= 0) par.criterio = criterio;
        r.FS = FactorSRMArc(*an, estado, nullptr, par, &r.saida);
        an->SetStrengthReductionFactor(1.);
        an->SetGravityFactor(1.);
        an->SetGravityScalesSuction(true);
        an->RestoreMemory(memoria);
        an->X() = estado;
        an->Restore();
        solo->SetFrozenPressures(false);
        r.seg = Segundos(t0);
        return r;
    }

private:
    void ColetaPontos() {
        for (TPZCompEl *cel : cmesh->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
            TPZManVector<int64_t, 10> idx;
            cel->GetMemoryIndices(idx);
            const TPZIntPoints &rule = cel->GetIntegrationRule();
            if (idx.size() != rule.NPoints()) DebugStop();
            TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
            TPZFNMatrix<9, REAL> jac, axes, jacinv;
            REAL w, detjac;
            for (int i = 0; i < rule.NPoints(); ++i) {
                rule.Point(i, qsi, w);
                cel->Reference()->X(qsi, x);
                cel->Reference()->Jacobian(qsi, jac, axes, detjac, jacinv);
                ips.push_back(idx[i]);
                ipel.push_back(cel->Reference()->Index());
                xs.push_back(x[0]);
                ys.push_back(x[1]);
                area.push_back(w * std::fabs(detjac));
                areaSolo += area.back();
            }
        }
    }
};

/// cotas das secoes da Fig. 3 (a cada 0,5 m, da base ao terreno)
inline std::vector<double> CotasSecao(double x) {
    std::vector<double> z;
    for (double v = 0.; v <= CotaTerreno(x) + 1.e-9; v += 0.5) z.push_back(v);
    return z;
}

/// estatistica de um campo nos pontos, ponderada pela area
struct Estat { double media = 0., dp = 0., min = 1.e300, max = -1.e300; };
template <class F> Estat Estatistica(const Talude &t, F f) {
    Estat e;
    double sa = 0., s1 = 0., s2 = 0.;
    for (size_t i = 0; i < t.ips.size(); ++i) {
        const double v = f(i), a = t.area[i];
        sa += a; s1 += a * v; s2 += a * v * v;
        e.min = std::min(e.min, v); e.max = std::max(e.max, v);
    }
    e.media = s1 / sa;
    e.dp = std::sqrt(std::max(0., s2 / sa - e.media * e.media));
    return e;
}

/// resultado de uma execucao (o que vai para resumo.csv e resultado.txt)
struct Resultado {
    InfoPermanente perm;
    InfoFS fs;
    std::vector<std::array<double, 3>> secoes; ///< (x, z, p_w)
    std::string situacao = "ok";
    double seg = 0.;
    int codigo = 0;
};

/// Roda uma execucao (uma realizacao) na pasta dir: os CSV, resumo.csv e por ultimo resultado.txt (um .tmp renomeado)
inline Resultado RodaExecucao(const Config &cf, const std::string &dir = ".") {
    const auto inicio = Relogio::now();
    Resultado res;
    namespace sfs = std::filesystem;
    sfs::create_directories(dir);
    auto arq = [&](const std::string &n) { return (sfs::path(dir) / n).string(); };
    std::remove(arq("resultado.txt").c_str());
    Talude t(cf);
    t.AplicaCampos();
    const Estat ec = Estatistica(t, [&](size_t i) { return t.CPonto(i); });
    const Estat ephi = Estatistica(t, [&](size_t i) { return t.PhiPonto(i); });
    const Estat ek = Estatistica(t, [&](size_t i) { return t.KPonto(i); });
    if (cf.camposCsv) {
        std::ofstream fc(arq("campos.csv"));
        fc << "ip,x,y,area,c,phi_deg,k_mult\n" << std::setprecision(7);
        for (size_t i = 0; i < t.ips.size(); ++i)
            fc << t.ips[i] << ',' << t.xs[i] << ',' << t.ys[i] << ',' << t.area[i] << ',' << t.CPonto(i) << ','
               << t.PhiPonto(i) << ',' << t.KPonto(i) << '\n';
    }
    if (Verboso())
        std::cout << "V4: caso " << cf.caso << ", angulo " << cf.angulo << ", real " << cf.real << ", campos "
                  << (cf.CampoK() || cf.CampoResist()) << ", neq " << t.an->NEq() << ", pontos " << t.ips.size()
                  << (t.klMUsado ? ", KL M = " + std::to_string(t.klMUsado) + " (" + Num(100. * t.klFracUsada, 4) + "% da variancia)" : "")
                  << std::endl;

    // ---- regime permanente (transiente longo com chuva constante, ou lido de PW_PERMANENTE quando k_s e constante)
    {
        std::ofstream fpas(arq("passos.csv"));
        bool lido = false;
        if (!cf.pwPermanente.empty()) {
            if (cf.CampoK()) std::cerr << "PW_PERMANENTE ignorado: k_s aleatorio neste caso\n";
            else if (!(lido = t.LePw(cf.pwPermanente))) std::cerr << "PW_PERMANENTE: " << cf.pwPermanente
                                                                 << " nao e um pw.csv desta malha; o transiente roda\n";
        }
        if (lido) {
            res.perm.ok = res.perm.lido = true;
            res.perm.psiPatamar = SuccaoPatamar(t.p.hydro, cf.qks);
            fpas << "lido_de\n" << cf.pwPermanente << '\n';
        } else res.perm = t.Permanente(&fpas);
    }
    if (!res.perm.ok) res.situacao = "permanente_nao_atingido_t_h=" + Num(res.perm.t_h, 6);
    if (Verboso())
        std::cout << "permanente: " << (res.perm.lido ? "lido" : res.perm.ok ? "atingido" : "NAO atingido") << " em t = "
                  << res.perm.t_h / 24. << " dias (" << res.perm.passos << " passos, " << res.perm.cortes << " cortes, "
                  << res.perm.seg << " s); max|dp_w| no ultimo passo " << res.perm.dpw << " kPa; chuva "
                  << res.perm.chuva << ", esquerda " << res.perm.esq << ", direita " << res.perm.dir << " m2/s"
                  << std::endl;
    if (cf.pwCsv) t.GravaPw(arq("pw.csv"));
    {
        std::ofstream fsec(arq("secoes.csv"));
        fsec << "secao,x,z,pw_kPa\n" << std::setprecision(8);
        for (double x : {kXSecaoX, kXSecaoY})
            for (double z : CotasSecao(x)) {
                const double pw = t.PwEm(x, z);
                res.secoes.push_back({x, z, pw});
                fsec << (x == kXSecaoX ? "X-X" : "Y-Y") << ',' << x << ',' << z << ',' << pw << '\n';
            }
    }

    // ---- fator de seguranca com as poropressoes do permanente congeladas
    {
        std::unique_ptr<std::ofstream> fmec;
        if (cf.mecCsv) {
            fmec = std::make_unique<std::ofstream>(arq("mecanismo.csv"));
            *fmec << "ip,x,y,area,du_m\n" << std::setprecision(7);
        }
        res.fs = t.FS(fmec.get());
        std::ofstream ften(arq("fs_tentativas.csv"));
        ften << "ordem,fase,g,F,convergiu,newton,nova_tentativa,seg\n" << std::setprecision(10);
        for (size_t k = 0; k < t.tentativas.size(); ++k) {
            const auto &q = t.tentativas[k];
            ften << k + 1 << ',' << q.fase << ',' << q.g << ',' << q.F << ',' << q.conv << ',' << q.iter << ',' << q.retry
                 << ',' << q.seg << '\n';
        }
    }
    const InfoFS &f = res.fs;
    if (!(f.FS > 0.) && res.situacao == "ok") res.situacao = "fs_abaixo_de_F0";
    {
        std::ofstream ffs(arq("fs.csv"));
        ffs << "t_h,FS_SRM,FS_parada,tentativas,newton_fs,u_max_m,du_max_m,massa_pts,massa_area_m2,massa_xc,massa_yc,"
               "massa_xmin,massa_xmax,massa_ymin,Df_m,modo,massa_c,massa_phi_deg,massa_k_mult,area_lim10_m2,"
               "area_lim20_m2,area_lim30_m2,area_lim50_m2,seg\n"
            << std::setprecision(10);
        auto a = [&](size_t k) { return k < f.areasLimiar.size() ? Num(f.areasLimiar[k].second) : std::string(); };
        const bool m = f.massaPts > 0;
        auto se = [&](bool ok, double v) { return ok ? Num(v) : std::string(); };
        ffs << res.perm.t_h << ',' << f.FS << ',' << f.Fparada << ',' << f.tentativas << ',' << f.newton << ','
            << se(f.FS > 0., f.uMax) << ',' << se(f.FS > 0., f.duMax) << ',' << f.massaPts << ',' << se(m, f.massaArea)
            << ',' << se(m, f.massaXc) << ',' << se(m, f.massaYc) << ',' << se(m, f.massaXmin) << ','
            << se(m, f.massaXmax) << ',' << se(m, f.massaYmin) << ',' << se(m, f.Df) << ',' << f.modo << ','
            << se(m, f.massaC) << ',' << se(m, f.massaPhi) << ',' << se(m, f.massaK) << ',' << a(0) << ',' << a(1) << ','
            << a(2) << ',' << a(3) << ',' << f.seg << '\n';
    }
    if (Verboso())
        std::cout << "FS (SRM) = " << f.FS << " (" << f.tentativas << " solucoes, " << f.newton << " iteracoes, " << f.seg
                  << " s); massa que desliza " << f.massaArea << " m2, D_f = " << f.Df << " m (" << f.modo << ")"
                  << std::endl;

    // ---- resumo.csv e resultado.txt
    res.seg = Segundos(inicio);
    auto pwSec = [&](double x, double z) {
        for (const auto &s : res.secoes) if (s[0] == x && std::fabs(s[1] - z) < 1.e-9) return s[2];
        return std::numeric_limits<double>::quiet_NaN();
    };
    double lx = 0., ly = 0.;
    cf.Comprimentos(lx, ly);
    std::vector<std::pair<std::string, std::string>> r = {
        {"real", std::to_string(cf.real)}, {"caso", cf.caso}, {"angulo", std::to_string(cf.angulo)},
        {"campos", std::to_string(int(cf.CampoK() || cf.CampoResist()))}, {"semente", std::to_string(cf.semente)},
        {"cov_c", Num(cf.covC)}, {"cov_phi", Num(cf.covPhi)}, {"cov_k", Num(cf.covK)}, {"rho", Num(cf.rho)},
        {"Lx", Num(lx)}, {"Ly", Num(ly)}, {"KL_M", std::to_string(t.klMUsado)}, {"KL_variancia", Num(t.klFracUsada)},
        {"p_kr", Num(cf.pKr)}, {"q_ks", Num(cf.qks)}, {"chi", Num(cf.Chi())}, {"area_solo_m2", Num(t.areaSolo)},
        {"c_med", Num(ec.media)}, {"c_dp", Num(ec.dp)}, {"phi_med", Num(ephi.media)}, {"phi_dp", Num(ephi.dp)},
        {"k_mult_med", Num(ek.media)}, {"k_mult_dp", Num(ek.dp)},
        {"t_perm_h", Num(res.perm.t_h)}, {"perm_ok", std::to_string(int(res.perm.ok))},
        {"perm_lido", std::to_string(int(res.perm.lido))}, {"dpw_ultimo_kPa", Num(res.perm.dpw)},
        {"psi_patamar_kPa", Num(res.perm.psiPatamar)}, {"q_chuva_m2s", Num(res.perm.chuva)},
        {"q_esq_m2s", Num(res.perm.esq)}, {"q_dir_m2s", Num(res.perm.dir)}, {"q_face_m2s", Num(res.perm.face)},
        {"pw_XX_z10", Num(pwSec(kXSecaoX, 10.))}, {"pw_XX_z5", Num(pwSec(kXSecaoX, 5.))},
        {"pw_XX_z0", Num(pwSec(kXSecaoX, 0.))}, {"pw_YY_z7.5", Num(pwSec(kXSecaoY, 7.5))},
        {"pw_YY_z0", Num(pwSec(kXSecaoY, 0.))},
        {"FS_SRM", Num(f.FS)}, {"FS_min", Num(f.FS)}, {"FS_parada", Num(f.Fparada)},
        {"massa_area_m2", f.massaPts ? Num(f.massaArea) : ""}, {"Df_m", f.massaPts ? Num(f.Df) : ""}, {"modo", f.modo},
        {"massa_xc", f.massaPts ? Num(f.massaXc) : ""}, {"massa_yc", f.massaPts ? Num(f.massaYc) : ""},
        {"situacao", res.situacao}, {"passos", std::to_string(res.perm.passos)},
        {"cortes", std::to_string(res.perm.cortes)}, {"newton", std::to_string(res.perm.newton)},
        {"newton_fs", std::to_string(f.newton)}, {"seg_perm", Num(res.perm.seg)}, {"seg_fs", Num(f.seg)},
        {"seg", Num(res.seg)}, {"versao", EnvS("V4_VERSAO", "-")}};
    for (auto &kv : r) { std::replace(kv.second.begin(), kv.second.end(), ',', ';'); std::replace(kv.second.begin(), kv.second.end(), ' ', '_'); }
    {
        std::ofstream o(arq("resumo.csv.tmp"));
        for (size_t i = 0; i < r.size(); ++i) o << (i ? "," : "") << r[i].first;
        o << '\n';
        for (size_t i = 0; i < r.size(); ++i) o << (i ? "," : "") << r[i].second;
        o << '\n';
    }
    std::rename(arq("resumo.csv.tmp").c_str(), arq("resumo.csv").c_str());
    {
        std::ofstream o(arq("resultado.txt.tmp"));
        o << "real=" << cf.real << " caso=" << cf.caso << " angulo=" << cf.angulo << " FS_min=" << Num(f.FS)
          << " massa_area_m2=" << (f.massaPts ? Num(f.massaArea) : "") << " Df_m=" << (f.massaPts ? Num(f.Df) : "")
          << " modo=" << f.modo << " situacao=" << res.situacao << " t_perm_h=" << Num(res.perm.t_h)
          << " passos=" << res.perm.passos << " seg=" << Num(res.seg) << '\n';
    }
    std::rename(arq("resultado.txt.tmp").c_str(), arq("resultado.txt").c_str());
    res.codigo = res.situacao == "ok" ? 0 : 1;
    if (Verboso())
        std::cout << "real " << cf.real << ": FS = " << f.FS << ", " << res.situacao << ", " << res.seg << " s" << std::endl;
    return res;
}

} // namespace v4
