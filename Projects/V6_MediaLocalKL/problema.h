#pragma once
// V6_MediaLocalKL: o problema de verificacao do gerador de Karhunen-Loeve (Comum/KLSeparavel.h) pela media do campo ao
// longo de um arco circular: o desvio padrao amostral da media no arco contra a reducao de variancia analitica
// (Eqs. 3, 14 e 15) e, onde o arco pode ser reconstruido, contra as Tabelas 5 e 6. Sem elementos finitos.
// Referencia: Yuan, Hu, Mei, Wang & Wang (2021), Slope reliability analysis based on curvilinear local averaging of a
// 2-D random field, Comput. Geotech. 137:104247.
//
// Os casos sao os taludes das Tabelas 5 a 8 (Fig. 9: beta = 40 graus com H = 10 a 50 m e H = 50 m com beta = 30 a 50
// graus) e o arco da Fig. 15. O campo KL cobre o retangulo do dominio da Fig. 9 (como num RFEM do talude), com o nucleo
// exp(-|dx| / Lx - |dy| / Ly) e, pela Eq. 13 do artigo, Lx = delta1 / 2 = 1,5 m e Ly = delta2 / 2 = 0,5 m.
// Os arcos das Tabelas 5 e 6 nao sao impressos: usamos o circulo critico pelo pe de Bishop simplificado (bishop.h),
// com gamma calibrado no FS = 0,9938 do talude H = 50 m, beta = 40 graus (Tabela 7), e o circulo digitalizado da Fig. 15.
#include "arco.h"
#include "bishop.h"
#include "klarco.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace v6 {

// ---------------------------------------------------------------------------------------------------------- dados
constexpr double kSigmaC = 8.68, kSigmaPhi = 2.99;        ///< desvios no ponto (Tabela 2)
constexpr double kDelta1 = 3.0, kDelta2 = 1.0;            ///< distancias de correlacao (Sec. 3.3, Eq. 13)
constexpr double kFSCalibracao = 0.9938;                  ///< Tabela 7, H = 50 m, beta = 40 graus
constexpr double kVarEq20C = 1.323, kVarEq20Phi = 0.158;  ///< variancias reduzidas da Eq. 20 (H = 50 m, beta = 40)
constexpr double kFig15A = -46.02, kFig15B = 116.06;      ///< centro do circulo da Fig. 15 (digitalizado, aproximado)
constexpr double kFig15PontosA = -47.25, kFig15PontosB = 112.42; ///< circulo pelos nos pretos da Fig. 15, pelo pe

struct Caso {
    std::string nome;
    Talude t;
    double scL = 0., sphiL = 0.; ///< Tabelas 5 e 6
    double FS = 0.;              ///< Tabelas 7 e 8
    bool fig15 = false;          ///< o arco e o da Fig. 15 (senao, o circulo critico de Bishop)
};

inline const std::vector<Caso> &Casos() {
    static const std::vector<Caso> c = {
        {"b40_h10", {10., 40.}, 2.46, 0.901, 2.210, false},  {"b40_h20", {20., 40.}, 1.88, 0.599, 1.546, false},
        {"b40_h30", {30., 40.}, 1.53, 0.571, 1.300, false},  {"b40_h40", {40., 40.}, 1.356, 0.462, 1.167, false},
        {"b40_h50", {50., 40.}, 1.23, 0.413, 0.9938, false}, {"h50_b30", {50., 30.}, 0.971, 0.344, 1.318, false},
        {"h50_b35", {50., 35.}, 1.12, 0.355, 1.191, false},  {"h50_b45", {50., 45.}, 1.28, 0.412, 0.8752, false},
        {"h50_b50", {50., 50.}, 1.29, 0.432, 0.8023, false}, {"fig15", {50., 45.}, 1.28, 0.412, 0.8752, true},
    };
    return c;
}

inline const Caso &BuscaCaso(const std::string &nome) {
    for (const auto &c : Casos())
        if (c.nome == nome) return c;
    throw std::invalid_argument("caso desconhecido: " + nome + " (b40_h10 ... b40_h50, h50_b30 ... h50_b50, fig15)");
}

/// gamma (kN/m3) com o qual o FS minimo de Bishop do talude H = 50 m, beta = 40 graus e 0,9938 (calculado uma vez)
inline double GammaCalibrado() {
    static const double g = CalibraGamma(Talude{50., 40.}, SoloMC{}, kFSCalibracao);
    return g;
}

struct ArcoCaso {
    Arco arco;
    double FSbishop = 0., gamma = 0.;
};

/// O arco do caso: o da Fig. 15 ou o circulo critico de Bishop com o gamma calibrado
inline ArcoCaso ArcoDoCaso(const Caso &c) {
    SoloMC s;
    s.gamma = GammaCalibrado();
    if (c.fig15) return {Arco::PeloPe(kFig15A, kFig15B, c.t.H), FSBishop(c.t, s, kFig15A, kFig15B), s.gamma};
    const auto cc = BuscaCirculoCritico(c.t, s);
    return {Arco::PeloPe(cc.a, cc.b, c.t.H), cc.FS, s.gamma};
}

// ---------------------------------------------------------------------------------------------- configuracao
inline std::string Ambiente(const char *n, const std::string &padrao) {
    const char *v = std::getenv(n);
    return (v && *v) ? std::string(v) : padrao;
}

struct Config {
    std::string caso = "b40_h10";   ///< CASO: um nome de Casos() ou "tabelas" (so as Eqs. 14/15 e Bishop, sem KL)
    std::string convencao = "eq13"; ///< CONVENCAO: eq13 (L = delta / 2, padrao) ou delta (L = delta)
    double fracao = 0.90;           ///< FRACAO: fracao de variancia da KL (Eq. 29 do artigo: >= 90%)
    int Mmax = 400000;              ///< MKL: limite de termos (o que vier primeiro)
    long N = 10000;                 ///< N: realizacoes
    uint64_t semente = 20210607;    ///< SEMENTE
    long nConferir = 20;            ///< NCONF: realizacoes em que a media no arco tambem e feita por KL2D::Campo
    double d1 = kDelta1, d2 = kDelta2; ///< DELTA1, DELTA2 (m)
    int verbose = 1;                ///< VERBOSE

    Correlacao Rho() const {
        return convencao == "delta" ? Correlacao::CIgualUmSobreDelta(d1, d2) : Correlacao::Eq13(d1, d2);
    }
    Correlacao RhoOutra() const {
        return convencao == "delta" ? Correlacao::Eq13(d1, d2) : Correlacao::CIgualUmSobreDelta(d1, d2);
    }
    static Config DoAmbiente() {
        Config c;
        c.caso = Ambiente("CASO", c.caso);
        c.convencao = Ambiente("CONVENCAO", c.convencao);
        c.fracao = std::stod(Ambiente("FRACAO", std::to_string(c.fracao)));
        c.Mmax = std::stoi(Ambiente("MKL", std::to_string(c.Mmax)));
        c.N = std::stol(Ambiente("N", std::to_string(c.N)));
        c.semente = std::stoull(Ambiente("SEMENTE", std::to_string(c.semente)));
        c.nConferir = std::stol(Ambiente("NCONF", std::to_string(c.nConferir)));
        c.d1 = std::stod(Ambiente("DELTA1", std::to_string(c.d1)));
        c.d2 = std::stod(Ambiente("DELTA2", std::to_string(c.d2)));
        c.verbose = std::stoi(Ambiente("VERBOSE", std::to_string(c.verbose)));
        if (c.convencao != "eq13" && c.convencao != "delta") throw std::invalid_argument("CONVENCAO: eq13 ou delta");
        return c;
    }
};

// ------------------------------------------------------------------------------------------------- sondas
/// Pontos onde o campo e guardado em cada realizacao: P0 = P(T/2) (meio do arco), defasagens horizontais e verticais
/// de 0,25 a 3 L, duas diagonais, tres pontos adiante no arco e quatro pontos do arco (variancia no ponto)
struct Sondas {
    std::vector<double> x, y;
    std::vector<std::string> nome;
    std::vector<std::pair<int, int>> pares; ///< (0, i): correlacao entre P0 e a sonda i
};

inline Sondas CriaSondas(const Arco &A, const Correlacao &rho) {
    Sondas s;
    double x0, y0;
    A.Ponto(0.5 * A.T, x0, y0);
    auto add = [&](double x, double y, const std::string &n) { s.x.push_back(x); s.y.push_back(y); s.nome.push_back(n); };
    add(x0, y0, "P0");
    for (double f : {0.25, 0.5, 1., 2., 3.}) add(x0 + f * rho.Lx, y0, "dx=" + std::to_string(f).substr(0, 4) + "Lx");
    for (double f : {0.25, 0.5, 1., 2., 3.}) add(x0, y0 + f * rho.Ly, "dy=" + std::to_string(f).substr(0, 4) + "Ly");
    add(x0 + 0.5 * rho.Lx, y0 + 0.5 * rho.Ly, "diag 0.5L");
    add(x0 + rho.Lx, y0 + rho.Ly, "diag 1L");
    for (double dl : {0.5, 1., 2.}) {
        double x, y;
        A.Ponto(0.5 * A.T + dl, x, y);
        add(x, y, "arco dl=" + std::to_string(dl).substr(0, 3) + "m");
    }
    const int npares = (int)s.x.size();
    for (double f : {0.125, 0.25, 0.75, 0.875}) {
        double x, y;
        A.Ponto(f * A.T, x, y);
        add(x, y, "arco l=" + std::to_string(f).substr(0, 5) + "T");
    }
    for (int i = 1; i < npares; ++i) s.pares.push_back({0, i});
    return s;
}

// ---------------------------------------------------------------------------------------------- um caso com KL
struct LinhaConvergencia {
    int M = 0;
    double fracao = 0., G2 = 0., varP0 = 0., rhoLx = 0., rhoLy = 0.;
};

struct ResultadoCaso {
    Caso caso;
    ArcoCaso ac;
    Correlacao rho;
    double G2 = 0., G2impressa = 0., G2outra = 0.; ///< Eq. 14; Eq. 15 impressa; a outra convencao de delta
    int M = 0;
    double fracao = 0., G2KL = 0., wmax = 0.;
    int nq = 0;
    Sondas sondas;
    std::vector<double> varKL;           ///< variancia exata do campo truncado em cada sonda
    std::vector<double> rhoNucleo, rhoKL; ///< por par
    std::vector<LinhaConvergencia> conv;
    ResultadoMC mc;
    double segKL = 0., segMC = 0., seg = 0.;
};

/// Constroi a KL no retangulo do talude, as medias dos termos no arco, as variancias exatas, a convergencia com M e o
/// Monte Carlo de N realizacoes (N = 0: sem Monte Carlo)
inline ResultadoCaso RodaCaso(const Config &cf, std::ostream *log = nullptr) {
    const auto t0 = std::chrono::steady_clock::now();
    auto seg = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    ResultadoCaso R;
    R.caso = BuscaCaso(cf.caso);
    R.ac = ArcoDoCaso(R.caso);
    R.rho = cf.Rho();
    const Arco &A = R.ac.arco;
    const Talude &t = R.caso.t;
    if (!ArcoDentro(A, t)) throw std::runtime_error("o arco do caso " + cf.caso + " sai do talude");
    R.G2 = FatorVariancia(A, R.rho);
    R.G2impressa = FatorVarianciaEq15Impressa(A, R.rho);
    R.G2outra = FatorVariancia(A, cf.RhoOutra());
    if (log)
        *log << "V6 " << cf.caso << ": arco centro (" << A.a << ", " << A.b << "), r = " << A.r << " m, T = " << A.T
             << " m, FS Bishop = " << R.ac.FSbishop << " (gamma = " << R.ac.gamma << "), Gamma = " << std::sqrt(R.G2)
             << "\n";
    // KL no retangulo do dominio da Fig. 9
    const klsep::KL2D kl(t.X0(), t.Y0(), t.X1() - t.X0(), t.Y1() - t.Y0(), R.rho.Lx, R.rho.Ly, cf.fracao, cf.Mmax);
    R.M = kl.M();
    R.fracao = kl.fracao;
    R.wmax = FrequenciaMaxima(kl);
    const QuadArco q = CriaQuadArcoParaKL(A, kl);
    R.nq = q.N();
    const auto F = MediasDosTermos(kl, q);
    const auto V = VarianciaMediaAcumulada(kl, F);
    R.G2KL = V.back();
    // sondas: variancias e correlacoes exatas do campo truncado
    R.sondas = CriaSondas(A, R.rho);
    const auto phi = TermosNosPontos(kl, R.sondas.x, R.sondas.y);
    for (size_t i = 0; i < R.sondas.x.size(); ++i) R.varKL.push_back(CovarianciaTruncada(kl, phi[i], phi[i]));
    for (const auto &p : R.sondas.pares) {
        const int i = p.first, j = p.second;
        R.rhoNucleo.push_back(R.rho(R.sondas.x[j] - R.sondas.x[i], R.sondas.y[j] - R.sondas.y[i]));
        R.rhoKL.push_back(CovarianciaTruncada(kl, phi[i], phi[j]) / std::sqrt(R.varKL[i] * R.varKL[j]));
    }
    // convergencia com o numero de termos (os primeiros m termos sao a KL de m termos)
    {
        std::vector<int> Ms;
        for (double m = 10.; m < R.M; m *= std::sqrt(10.)) Ms.push_back((int)std::lround(m));
        Ms.push_back(R.M);
        const int iLx = 3, iLy = 8; // sondas dx = 1 Lx e dy = 1 Ly
        double soma = 0., c00 = 0., c0x = 0., cxx = 0., c0y = 0., cyy = 0.;
        size_t im = 0;
        const double area = kl.Dx * kl.Dy;
        for (int k = 0; k < R.M && im < Ms.size(); ++k) {
            const double l = kl.raiz[k] * kl.raiz[k];
            soma += l;
            c00 += l * phi[0][k] * phi[0][k];
            c0x += l * phi[0][k] * phi[iLx][k];
            cxx += l * phi[iLx][k] * phi[iLx][k];
            c0y += l * phi[0][k] * phi[iLy][k];
            cyy += l * phi[iLy][k] * phi[iLy][k];
            if (k + 1 == Ms[im]) {
                R.conv.push_back({k + 1, soma / area, V[k], c00, c0x / std::sqrt(c00 * cxx), c0y / std::sqrt(c00 * cyy)});
                ++im;
            }
        }
    }
    R.segKL = seg();
    if (log)
        *log << "V6 " << cf.caso << ": KL com M = " << R.M << " termos (fracao " << R.fracao << "), " << R.nq
             << " pontos no arco, Gamma_KL = " << std::sqrt(R.G2KL) << " (" << R.segKL << " s)\n";
    if (cf.N > 0) {
        R.mc = MonteCarloArco(kl, F, q, cf.N, cf.semente, R.sondas.x, R.sondas.y, R.sondas.pares, cf.nConferir);
        R.segMC = seg() - R.segKL;
        if (log)
            *log << "V6 " << cf.caso << ": Monte Carlo com " << cf.N << " realizacoes: dp da media no arco = "
                 << R.mc.media.DP() << " (+- " << EPdoDP(R.mc.media.DP(), cf.N) << "), Gamma = " << std::sqrt(R.G2)
                 << " (" << R.segMC << " s)\n";
    }
    R.seg = seg();
    return R;
}

// --------------------------------------------------------------------------------------- tabelas (sem KL)
struct LinhaTabela {
    Caso caso;
    ArcoCaso ac;
    double G = 0., Gimpressa = 0., Goutra = 0.; ///< sqrt do fator de variancia: Eq. 14 (eq13), Eq. 15 impressa, L = delta
};

/// Os arcos (Bishop ou Fig. 15) e os fatores de reducao dos 10 casos, para comparar com as Tabelas 5 a 8
inline std::vector<LinhaTabela> CalculaTabelas(const Config &cf) {
    std::vector<LinhaTabela> L;
    for (const auto &c : Casos()) {
        LinhaTabela l;
        l.caso = c;
        l.ac = ArcoDoCaso(c);
        l.G = std::sqrt(FatorVariancia(l.ac.arco, cf.Rho()));
        l.Gimpressa = std::sqrt(FatorVarianciaEq15Impressa(l.ac.arco, cf.Rho()));
        l.Goutra = std::sqrt(FatorVariancia(l.ac.arco, cf.RhoOutra()));
        L.push_back(l);
    }
    return L;
}

// ------------------------------------------------------------------------------------------------- saidas
inline void GravaResultado(const std::string &linha) { // resultado.txt por ultimo (um .tmp renomeado)
    {
        std::ofstream f("resultado.txt.tmp");
        f << linha << "\n";
    }
    std::rename("resultado.txt.tmp", "resultado.txt");
}

inline void GravaTabelas(const std::vector<LinhaTabela> &L, const Config &cf, double seg) {
    std::ofstream f("tabelas.csv");
    f << std::setprecision(8);
    f << "caso,H_m,beta_graus,arco,a_m,b_m,r_m,T_m,alpha_rad,theta_rad,gamma_kNm3,FS_bishop,FS_artigo,Gamma,"
         "sigma_cL,sigma_cL_artigo,sigma_phiL,sigma_phiL_artigo,Gamma_c_artigo,Gamma_phi_artigo,Gamma_eq15_impressa,"
         "Gamma_L_igual_delta\n";
    for (const auto &l : L) {
        const Arco &A = l.ac.arco;
        f << l.caso.nome << "," << l.caso.t.H << "," << l.caso.t.beta << "," << (l.caso.fig15 ? "fig15" : "bishop") << ","
          << A.a << "," << A.b << "," << A.r << "," << A.T << "," << A.alpha << "," << A.theta << "," << l.ac.gamma
          << "," << l.ac.FSbishop << "," << l.caso.FS << "," << l.G << "," << l.G * kSigmaC << "," << l.caso.scL << ","
          << l.G * kSigmaPhi << "," << l.caso.sphiL << "," << l.caso.scL / kSigmaC << "," << l.caso.sphiL / kSigmaPhi
          << "," << l.Gimpressa << "," << l.Goutra << "\n";
    }
    std::ostringstream r;
    r << "caso=tabelas convencao=" << cf.convencao << " gamma=" << GammaCalibrado() << " n=" << L.size()
      << " situacao=ok seg=" << seg;
    {
        std::ofstream s("resumo.csv");
        s << "caso,convencao,gamma_kNm3,n,seg\ntabelas," << cf.convencao << "," << GammaCalibrado() << "," << L.size()
          << "," << seg << "\n";
    }
    GravaResultado(r.str());
}

inline void GravaCaso(const ResultadoCaso &R, const Config &cf) {
    const Arco &A = R.ac.arco;
    const double G = std::sqrt(R.G2), GKL = std::sqrt(R.G2KL);
    const bool mc = R.mc.N > 0;
    const double sd = mc ? R.mc.media.DP() : 0., ep = mc ? EPdoDP(sd, R.mc.N) : 0.;
    {
        std::ofstream f("convergencia.csv");
        f << std::setprecision(10) << "M,fracao,Gamma_KL,Gamma_KL_sobre_Gamma,var_P0,rho_KL_dx_Lx,rho_KL_dy_Ly,rho_nucleo_1L\n";
        for (const auto &c : R.conv)
            f << c.M << "," << c.fracao << "," << std::sqrt(c.G2) << "," << std::sqrt(c.G2) / G << "," << c.varP0 << ","
              << c.rhoLx << "," << c.rhoLy << "," << std::exp(-1.) << "\n";
    }
    {
        std::ofstream f("sondas.csv");
        f << std::setprecision(10) << "sonda,x_m,y_m,var_KL,var_MC,ep_var_MC,media_MC\n";
        for (size_t i = 0; i < R.sondas.x.size(); ++i) {
            f << R.sondas.nome[i] << "," << R.sondas.x[i] << "," << R.sondas.y[i] << "," << R.varKL[i];
            if (mc) f << "," << R.mc.sonda[i].Var() << "," << EPdaVar(R.varKL[i], R.mc.N) << "," << R.mc.sonda[i].media;
            else f << ",,,";
            f << "\n";
        }
    }
    {
        std::ofstream f("correlacao.csv");
        f << std::setprecision(10) << "par,dx_m,dy_m,rho_nucleo,rho_KL,rho_MC,ep_rho_MC\n";
        for (size_t p = 0; p < R.sondas.pares.size(); ++p) {
            const int i = R.sondas.pares[p].first, j = R.sondas.pares[p].second;
            f << R.sondas.nome[j] << "," << R.sondas.x[j] - R.sondas.x[i] << "," << R.sondas.y[j] - R.sondas.y[i] << ","
              << R.rhoNucleo[p] << "," << R.rhoKL[p];
            if (mc) f << "," << R.mc.par[p].Rho() << "," << EPdaCorrelacao(R.rhoKL[p], R.mc.N);
            else f << ",,";
            f << "\n";
        }
    }
    std::ostringstream cab, lin;
    cab << "caso,convencao,Lx_m,Ly_m,H_m,beta_graus,a_m,b_m,r_m,T_m,FS_bishop,gamma_kNm3,M,fracao,pontos_arco,"
           "Gamma_eq14,Gamma_KL,Gamma_eq15_impressa,Gamma_outra_convencao,N,semente,dp_MC,ep_dp_MC,media_MC,"
           "dp_MC_sobre_Gamma,dp_MC_sobre_Gamma_KL,var_P0_KL,var_P0_MC,max_dif_Campo,sigma_cL,sigma_cL_artigo,"
           "sigma_phiL,sigma_phiL_artigo,assimetria_MC,excesso_curtose_MC,fracao_fora_1.96,seg_KL,seg_MC,seg";
    lin << std::setprecision(10) << R.caso.nome << "," << cf.convencao << "," << R.rho.Lx << "," << R.rho.Ly << ","
        << R.caso.t.H << "," << R.caso.t.beta << "," << A.a << "," << A.b << "," << A.r << "," << A.T << ","
        << R.ac.FSbishop << "," << R.ac.gamma << "," << R.M << "," << R.fracao << "," << R.nq << "," << G << "," << GKL
        << "," << std::sqrt(R.G2impressa) << "," << std::sqrt(R.G2outra) << "," << R.mc.N << "," << cf.semente << ","
        << sd << "," << ep << "," << (mc ? R.mc.media.media : 0.) << "," << (mc ? sd / G : 0.) << ","
        << (mc ? sd / GKL : 0.) << "," << R.varKL[0] << "," << (mc ? R.mc.sonda[0].Var() : 0.) << ","
        << R.mc.maxDifCampo << "," << G * kSigmaC << "," << R.caso.scL << "," << G * kSigmaPhi << "," << R.caso.sphiL
        << "," << (mc ? R.mc.Assimetria() : 0.) << "," << (mc ? R.mc.ExcessoCurtose() : 0.) << ","
        << (mc ? (double)R.mc.nFora / R.mc.N : 0.) << "," << R.segKL << "," << R.segMC << "," << R.seg;
    {
        std::ofstream f("resumo.csv");
        f << cab.str() << "\n" << lin.str() << "\n";
    }
    std::ostringstream r;
    r << std::setprecision(8) << "caso=" << R.caso.nome << " convencao=" << cf.convencao << " M=" << R.M
      << " fracao=" << R.fracao << " Gamma=" << G << " Gamma_KL=" << GKL << " N=" << R.mc.N << " dp_MC=" << sd
      << " ep_dp_MC=" << ep << " situacao=ok seg=" << R.seg;
    GravaResultado(r.str());
}

} // namespace v6
