#pragma once
// V5_EncostaDuasCamadas: encosta em duas camadas de Abbasov, Fahs, Younes, Nowamooz, Maloy e Toussaint (2024),
// Comput. Geotech. 167:106102 (o benchmark de Lu et al. 2012 / Shao et al. 2015): 10 anos de chuva a 600 mm/ano a partir
// do estado hidrostatico (nivel d'agua em z = 17 m), depois 20 mm/h por 15 h; fluxo de Richards nas duas camadas,
// chi = S_e; campo de tensoes linear elastico e o fator de seguranca local LFS de Lu et al. (2012) (eq. 8) na zona de
// investigacao: a media do LFS e a area com LFS < 1 (eqs. 11 a 13).
//
// Modelo (ver README.md):
//  - material TPZMatUnsaturatedUPwPa com p_a = 0 imposto em todos os nos (so a equacao de Richards) e esqueleto rigido
//    no fluxo (todos os deslocamentos impostos): o fluxo e o de um modelo de Richards num dominio fixo, como no artigo;
//  - coeficiente de Biot alpha = 0: as poropressoes nao entram no equilibrio, que fica div(sigma) + gamma(theta) g = 0
//    com sigma = D eps (tensao TOTAL elastica, eq. 5 do artigo; abordagem "hillslope cube"), peso da mistura
//    (1 - theta_s) rho_s g + theta rho_w g com a umidade do instante; depois sigma' = sigma - chi p_w I (eq. 7) so no LFS;
//  - Mohr-Coulomb com coesao de 1e9 kPa (o calculo fica elastico); c' e phi' do artigo entram so no LFS;
//  - retencao: Brooks-Corey (eq. 3) nao existe na biblioteca; van Genuchten ajustado a ela em 0 < h <= 1,5 m
//    (scripts/ajuste_vg.py); condutividade: a de Brooks-Corey exatamente, k_rw = S_e^(l + 2 + 2/n_BC) (eq. 4);
//  - contornos (Fig. 2): base e lateral esquerda impermeaveis; lateral direita com carga hidrostatica (nivel 17 m)
//    abaixo de z = 17 m e face de percolacao acima; superficie com a chuva (fluxo vertical por area horizontal, o
//    "fluxo x cos(inclinacao)" do artigo) e face de percolacao quando satura; apoios de rolete na base e nas laterais.
//
// Uso: main.cpp (o caso completo, configuravel pelo ambiente) e testes.cpp (Catch2).
#include "UnsatAnalysis.h" // UnsatModel.h (UnsatParams, UnsatMaterial, UnsatMesh com um material por id) e TPZUnsatAnalysis
#include "TPZGmshReader.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef VERIF_DIR
#define VERIF_DIR "."
#endif

namespace v5 {

// ------------------------------------------------------------------------------------------------------------------
// ambiente e utilidades
inline double EnvD(const char *k, double d) { const char *v = std::getenv(k); return v && *v ? std::atof(v) : d; }
inline int EnvI(const char *k, int d) { const char *v = std::getenv(k); return v && *v ? std::atoi(v) : d; }
inline std::string EnvS(const char *k, const std::string &d) { const char *v = std::getenv(k); return v && *v ? std::string(v) : d; }
inline std::vector<double> EnvLista(const char *k, const std::vector<double> &d) {
    const char *v = std::getenv(k);
    if (!v || !*v) return d;
    std::vector<double> r;
    std::stringstream ss(v);
    std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) r.push_back(std::atof(t.c_str()));
    return r;
}
inline std::string Num(double v, int p = 10) { std::ostringstream o; o << std::setprecision(p) << v; return o.str(); }
using Relogio = std::chrono::steady_clock;
inline double Segundos(Relogio::time_point t0) { return std::chrono::duration<double>(Relogio::now() - t0).count(); }
/// nivel de mensagens na tela (0: nada; 1: etapas e saidas; 2: cada passo)
inline int &Verboso() { static int v = 1; return v; }

constexpr double kHora = 3600., kDia = 86400., kAno = 365.25 * 86400.;

// ------------------------------------------------------------------------------------------------------------------
// geometria (Fig. 2; ver malhas/encosta.geo): dominio 42 m x 25 m, crista z = 25 m ate x = 19 m, talude de 6 m em 15 m
// (21,8 graus; o rotulo "23 graus" da Fig. 2 nao fecha com os 22 m e os 113 m2 da zona), pe em x = 34 m, z = 19 m;
// camada superior de 2 m (espessura vertical); zona de investigacao 15 <= x <= 37 m, z >= 17 m
constexpr double kLarg = 42., kZCrista = 25., kXCrista = 19., kXPe = 34., kZPe = 19., kEspSup = 2.;
constexpr double kXZona0 = 15., kXZona1 = 37., kZZona0 = 17.;
/// cota do terreno em x
inline double CotaTerreno(double x) {
    if (x <= kXCrista) return kZCrista;
    if (x >= kXPe) return kZPe;
    return kZCrista - (kZCrista - kZPe) * (x - kXCrista) / (kXPe - kXCrista);
}
inline double CotaInterface(double x) { return CotaTerreno(x) - kEspSup; }
/// ponto (centroide de um elemento) na zona de investigacao
inline bool NaZonaEncosta(double x, double z) {
    return x >= kXZona0 - 1.e-9 && x <= kXZona1 + 1.e-9 && z >= kZZona0 - 1.e-9;
}

/// ids dos dominios (grupos fisicos "Superior" e "Inferior"); contornos: os EBnd do UnsatModel.h (base -1, direita -2,
/// superficie -3, esquerda -4)
enum ECamada { kSuperior = 1, kInferior = 2 };

// ------------------------------------------------------------------------------------------------------------------
/// Uma camada: dados de Brooks-Corey da Tabela 1, o ajuste de van Genuchten e a coesao efetiva
struct Camada {
    std::string nome;
    double thetaS = 0.4, thetaR = 0.05; ///< teor de umidade saturado (= porosidade) e residual
    double ksCmH = 1.;                  ///< condutividade saturada (cm/h)
    double alphaBC = 0.05, nBC = 0.3;   ///< Brooks-Corey: alpha (1/cm) e n
    double lBC = 1.;                    ///< l de Brooks-Corey (o artigo da um valor so, 1)
    double aVG = 0.1, bVG = 1.5;        ///< van Genuchten ajustado: a (1/kPa) e b (m = 1 - 1/b)
    double c = 5.;                      ///< coesao efetiva (kPa), so no LFS
    double Ks() const { return ksCmH / 100. / 3600.; }     ///< m/s
    double PotenciaK() const { return lBC + 2. + 2. / nBC; } ///< expoente de K = K_s S_e^p (eq. 4)
    double Swr() const { return thetaR / thetaS; }          ///< S_wr de van Genuchten (S_w = theta / theta_s)
    /// S_e de Brooks-Corey para a succao h (m de agua), eq. (3)
    double SeBC(double h) const {
        const double x = alphaBC * 100. * h;
        return x > 1. ? std::pow(x, -nBC) : 1.;
    }
    /// S_e de van Genuchten para a succao p_c (kPa)
    double SeVG(double pc) const {
        if (pc <= 0.) return 1.;
        return std::pow(1. + std::pow(aVG * pc, bVG), -(1. - 1. / bVG));
    }
    double Theta(double Se) const { return thetaR + (thetaS - thetaR) * Se; }
};

// ajuste de van Genuchten a Brooks-Corey em 0 < h <= 1,5 m (scripts/ajuste_vg.py; RMSE de S_e 0,023 e 0,011)
constexpr double kAVGSup = 0.33500, kBVGSup = 1.4762; // camada superior (a em 1/kPa)
constexpr double kAVGInf = 0.05630, kBVGInf = 1.5664; // camada inferior

/// Parametros de uma execucao (CasoPadrao e DoAmbiente; os valores padrao sao os do artigo e as escolhas do README.md)
struct Config {
    /// tabela1: Tabela 1 (mapas das Figs. 4 a 6: nu = 0,35, c' = 3 / 6 kPa, K_s(l) = 0,06 cm/h);
    /// tabela4: valores deterministicos da Tabela 4 (Fig. 14 e sec. 4.6: nu = 0,3, c' = 4,5 kPa, K_s(l) = 0,05 cm/h);
    /// fig9: Fig. 9 (K_s(l) = 0,06 cm/h, nu = 0,35, c' = 5 kPa)
    std::string caso = "tabela1";
    Camada sup, inf;
    double E = 1.e4, nu = 0.35, phiDeg = 35.; ///< E = 10 MPa (kPa); E nao muda as tensoes (so forcas de volume)
    double rhoS = 2636., g = 9.81, gw = 9.81; ///< kg/m3, m/s2, kN/m3
    double Kw = 2.2e6;                         ///< modulo volumetrico da agua (kPa): armazenamento na zona saturada
    // ---- chuva e tempos (sec. 2.3)
    double anos = 10., chuva1MmAno = 600., chuva2MmH = 20., horas = 15.;
    double nivel = 17.; ///< nivel d'agua inicial (hidrostatico) e carga na lateral direita abaixo de zFace
    double zFace = 17.; ///< lateral direita: carga abaixo, face de percolacao acima (leitura da Fig. 2)
    // ---- passos: etapa 1 com dt crescente (de DT0_H por FATOR_DT ate DT_MAX_D); etapa 2 com DT2_H fixo
    double dt0_h = 1., fatorDt = 1.5, dtMax_d = 30., dt2_h = 0.1;
    int maxCortes = 12;
    std::vector<double> saidasH = {0., 7.5, 15.};
    double tolU = 1.e-7, tolP = 1.e-6;
    int maxIt = 40;
    // ---- malha e arquivos
    std::string malha;       ///< MALHA: padrao malhas/encosta.msh do projeto (um nome sem '/' e procurado em malhas/)
    std::string estadoIni;   ///< ESTADO_INICIAL: le o estado do fim da etapa 1 (estado.csv) em vez de calcula-lo
    std::string gravaEstado; ///< GRAVA_ESTADO: grava o estado do fim da etapa 1
    int vtk = 0, csvPontos = 0, csvCampos = 1;

    double Chuva1() const { return chuva1MmAno / 1000. / kAno; } ///< m/s
    double Chuva2() const { return chuva2MmH / 1000. / kHora; }  ///< m/s
    double Phi() const { return phiDeg * M_PI / 180.; }
    const Camada &DaCamada(int id) const { return id == kSuperior ? sup : inf; }

    static Config CasoPadrao(const std::string &caso) {
        Config cf;
        cf.caso = caso;
        // Tabela 1: theta_s, theta_r, K_s (cm/h), alpha_BC (1/cm), n_BC, l_BC; c' das colunas de cada camada
        cf.sup = {"superior", 0.412, 0.041, 2.59, 0.068, 0.322, 1., kAVGSup, kBVGSup, 3.};
        cf.inf = {"inferior", 0.385, 0.09, 0.06, 0.027, 0.131, 1., kAVGInf, kBVGInf, 6.};
        if (caso == "tabela1") {
        } else if (caso == "tabela4") { // Tabela 4 (e Tabela 3 para o resto)
            cf.inf.ksCmH = 0.05, cf.nu = 0.3, cf.phiDeg = 35., cf.sup.c = cf.inf.c = 4.5;
        } else if (caso == "fig9") {
            cf.inf.ksCmH = 0.06, cf.nu = 0.35, cf.phiDeg = 35., cf.sup.c = cf.inf.c = 5.;
        } else
            throw std::invalid_argument("CASO deve ser tabela1, tabela4 ou fig9");
        return cf;
    }

    /// o caso de CASO e as trocas do ambiente (nomes no README.md)
    static Config DoAmbiente() {
        Config cf = CasoPadrao(EnvS("CASO", "tabela1"));
        cf.nu = EnvD("NU", cf.nu);
        cf.phiDeg = EnvD("PHI_DEG", cf.phiDeg);
        cf.E = EnvD("E", cf.E);
        cf.sup.c = EnvD("C_SUP", EnvD("COES", cf.sup.c));
        cf.inf.c = EnvD("C_INF", EnvD("COES", cf.inf.c));
        cf.sup.ksCmH = EnvD("KS_SUP", cf.sup.ksCmH);
        cf.inf.ksCmH = EnvD("KS_INF", cf.inf.ksCmH);
        cf.sup.aVG = EnvD("A_VG_SUP", cf.sup.aVG);
        cf.sup.bVG = EnvD("B_VG_SUP", cf.sup.bVG);
        cf.inf.aVG = EnvD("A_VG_INF", cf.inf.aVG);
        cf.inf.bVG = EnvD("B_VG_INF", cf.inf.bVG);
        cf.rhoS = EnvD("RHO_S", cf.rhoS);
        cf.anos = EnvD("ANOS", cf.anos);
        cf.chuva1MmAno = EnvD("CHUVA1_MM_ANO", cf.chuva1MmAno);
        cf.chuva2MmH = EnvD("CHUVA2_MM_H", cf.chuva2MmH);
        cf.horas = EnvD("HORAS", cf.horas);
        cf.nivel = EnvD("NIVEL", cf.nivel);
        cf.zFace = EnvD("Z_FACE", cf.zFace);
        cf.dt0_h = EnvD("DT0_H", cf.dt0_h);
        cf.fatorDt = EnvD("FATOR_DT", cf.fatorDt);
        cf.dtMax_d = EnvD("DT_MAX_D", cf.dtMax_d);
        cf.dt2_h = EnvD("DT2_H", cf.dt2_h);
        cf.maxCortes = EnvI("MAX_CORTES", cf.maxCortes);
        cf.saidasH = EnvLista("SAIDAS_H", cf.saidasH);
        cf.tolU = EnvD("TOLU", cf.tolU);
        cf.tolP = EnvD("TOLP", cf.tolP);
        cf.maxIt = EnvI("MAXIT", cf.maxIt);
        cf.malha = EnvS("MALHA", cf.malha);
        cf.estadoIni = EnvS("ESTADO_INICIAL", cf.estadoIni);
        cf.gravaEstado = EnvS("GRAVA_ESTADO", cf.gravaEstado);
        cf.vtk = EnvI("VTK", cf.vtk);
        cf.csvPontos = EnvI("PONTOS_CSV", cf.csvPontos);
        cf.csvCampos = EnvI("CAMPOS_CSV", cf.csvCampos);
        if (!(cf.nu > -1. && cf.nu < 0.5) || !(cf.dt0_h > 0.) || !(cf.dt2_h > 0.) || !(cf.fatorDt >= 1.))
            throw std::invalid_argument("NU, DT0_H, DT2_H ou FATOR_DT fora da faixa");
        std::sort(cf.saidasH.begin(), cf.saidasH.end());
        return cf;
    }
};

/// caminho da malha: o dado, ou malhas/<nome> do projeto (padrao encosta.msh)
inline std::string CaminhoMalha(const std::string &m) {
    if (m.empty()) return std::string(VERIF_DIR) + "/malhas/encosta.msh";
    if (m.find('/') == std::string::npos && !std::filesystem::exists(m)) return std::string(VERIF_DIR) + "/malhas/" + m;
    return m;
}

/// Le a malha do Gmsh (.msh 4.1) com os grupos fisicos Superior, Inferior, Base, Direita, Superficie e Esquerda
inline TPZGeoMesh *LeMalha(const std::string &arq) {
    if (!std::filesystem::exists(arq)) throw std::runtime_error("malha nao encontrada: " + arq);
    TPZManVector<std::map<std::string, int>, 4> nomes(4);
    nomes[2] = {{"Superior", kSuperior}, {"Inferior", kInferior}};
    nomes[1] = {{"Base", EBase}, {"Direita", ERight}, {"Superficie", ETop}, {"Esquerda", ELeft}};
    TPZGmshReader leitor;
    leitor.SetDimNamePhysical(nomes);
    TPZGeoMesh *g = leitor.GeometricGmshMesh(arq);
    for (TPZGeoEl *gel : g->ElementVec()) {
        if (!gel) continue;
        bool ok = false;
        for (auto &n : nomes[gel->Dimension()]) ok = ok || n.second == gel->MaterialId();
        if (!ok) throw std::runtime_error("elemento fora dos grupos fisicos em " + arq);
    }
    g->BuildConnectivity();
    return g;
}

/// Parametros do material de uma camada (UnsatParams do UnsatModel.h)
inline UnsatParams ParametrosCamada(const Config &cf, const Camada &L) {
    UnsatParams p;
    p.E = cf.E, p.nu = cf.nu, p.phi = cf.Phi();
    p.c = 1.e9;                         // resistencia muito alta: o calculo de tensoes e linear elastico
    p.mixture = true;                   // peso (1 - n) gamma_s + n S_w gamma_w = (1 - theta_s) rho_s g + theta rho_w g
    p.gammas = cf.rhoS * cf.g / 1000.;  // kN/m3
    p.bishop = TUnsat::EChiEffective;   // chi = S_e (so no LFS: com alpha = 0 p_B nao entra no equilibrio)
    p.n = L.thetaS, p.Kw = cf.Kw, p.gw = cf.gw, p.ks = L.Ks();
    p.alpha = 0., p.Ks = 0.;            // tensao total elastica (eq. 5): as poropressoes nao entram no equilibrio
    p.coupled = false, p.lumped = true; // fluxo de Richards num dominio fixo
    p.hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    p.hydro.fAVG = L.aVG, p.hydro.fBVG = L.bVG, p.hydro.fSwr = L.Swr();
    p.hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw; // K = K_s S_e^p, a de Brooks-Corey
    p.hydro.fPower = L.PotenciaK();
    p.hydro.fKra = TPZUnsaturatedHydraulics::EMualemAir; // sem efeito: p_a = 0 imposto
    return p;
}

/// LFS de Lu et al. (2012), eq. (8): cos(phi') (c' + sigma'_I tan(phi')) / sigma_II, com o centro sigma'_I e o raio
/// sigma_II do circulo de Mohr no plano (compressao positiva). sxx, syy, sxy: tensao TOTAL com tracao positiva (a da
/// memoria do material); p_w > 0 em compressao; sigma' = sigma - chi p_w I (eq. 7) com chi = S_e
inline double LFSdeTensao(double sxx, double syy, double sxy, double pw, double Se, double c, double phi,
                          double *centro = nullptr, double *raio = nullptr) {
    const double a = -sxx, b = -syy, t = -sxy; // compressao positiva
    const double sI = 0.5 * (a + b) - Se * pw;
    const double sII = std::sqrt(0.25 * (a - b) * (a - b) + t * t);
    if (centro) *centro = sI;
    if (raio) *raio = sII;
    return std::cos(phi) * (c + sI * std::tan(phi)) / std::max(sII, 1.e-12);
}

// ------------------------------------------------------------------------------------------------------------------
/// balanco de agua e convergencia de um passo aceito (volumes por metro de encosta)
struct InfoPasso {
    double dt = 0., chuva = 0., face = 0., dir = 0., dV = 0., residuo = 0.;
    double pwMin = 0., pwMax = 0., dpwMax = 0., seg = 0.;
    int newton = 0, itsFace = 0, faceNos = 0;
};

/// metricas do LFS na zona de investigacao (eqs. 11 a 13)
struct Metricas {
    double t_h = 0., lfsMedio = 0., aLFS = 0., areaZona = 0., lfsMin = 0., lfsMedioLim = 0., lfsMedioLim3 = 0.;
    double aLFSsup = 0., aLFSinf = 0., aLFS3 = 0.; ///< area com LFS < 1 em cada camada e com LFS > 3
    double thetaSupMin = 0., thetaSupMax = 0., thetaInfMin = 0., thetaInfMax = 0., sigsMax = 0.;
};

/// O modelo: malha, analise, contornos, pontos de integracao, as etapas de chuva e o LFS
class Encosta {
public:
    Config cf;
    TPZGeoMesh *gmesh = nullptr;
    TPZMultiphysicsCompMesh *mphys = nullptr;
    std::unique_ptr<TPZUnsatAnalysis> an;
    std::map<int, TUnsat *> mat;          ///< id da camada -> material
    std::map<int64_t, REAL> chuvaNo;      ///< nos da superficie -> comprimento tributario horizontal (m)
    std::map<int64_t, REAL> compNo;       ///< candidatos a face de percolacao -> comprimento tributario real (m)
    std::map<int64_t, REAL> contorno;     ///< equacoes impostas no fluxo (u, p_a e p_w da carga na lateral direita)
    std::set<int64_t> eqDir;              ///< equacoes de p_w com carga imposta na lateral direita
    std::set<int64_t> face;               ///< nos com p_w = 0 (face de percolacao ativa)
    /// ponto de integracao: camada, memoria, elemento geometrico, coordenadas, peso (w det J) e se o elemento esta na zona
    struct Ponto {
        int camada;
        int64_t mem, el;
        double x, z, w;
        bool zona;
    };
    std::vector<Ponto> pontos;
    /// estado mecanico de cada ponto do ultimo calculo de tensoes (tensao total, tracao positiva) e o LFS
    struct EstadoPonto { double sxx = 0., syy = 0., sxy = 0., szz = 0., pw = 0., S = 1., Se = 1., lfs = 0.; };
    std::vector<EstadoPonto> estado;
    TPZFMatrix<STATE> xMec; ///< solucao do ultimo calculo de tensoes (deslocamentos)
    double areaZona = 0.;
    double tempo = 0.;      ///< tempo desde o inicio da etapa (s)
    std::function<bool(double, double)> naZona;
    // contadores
    int passos = 0, cortes = 0, newtonTotal = 0;
    InfoPasso ultimo;

    Encosta(const Config &c, TPZGeoMesh *g, std::function<bool(double, double)> zona = NaZonaEncosta)
        : cf(c), gmesh(g), naZona(std::move(zona)) {
        mphys = UnsatMesh(gmesh, {{kSuperior, ParametrosCamada(cf, cf.sup)}, {kInferior, ParametrosCamada(cf, cf.inf)}});
        // apoios de rolete (Fig. 2): u_y = 0 na base, u_x = 0 nas laterais
        an = std::make_unique<TPZUnsatAnalysis>(mphys, std::map<int, std::array<bool, 2>>{
            {EBase, {false, true}}, {ELeft, {true, false}}, {ERight, {true, false}}});
        an->SetTolerances(cf.tolU, cf.tolP, cf.maxIt, 25);
        for (int id : {kSuperior, kInferior}) {
            auto *m = dynamic_cast<TUnsat *>(mphys->FindMaterial(id));
            if (m) mat[id] = m;
        }
        if (mat.empty()) throw std::runtime_error("malha sem as camadas Superior/Inferior");
        // superficie: chuva pela projecao horizontal; candidatos a face: superficie e lateral direita acima de zFace
        chuvaNo = an->BoundaryNodes(ETop, true);
        for (const auto &b : an->BoundaryNodes(ETop, false)) compNo[b.first] += b.second;
        for (const auto &b : an->BoundaryNodes(ERight, false))
            if (an->Coord(b.first, 1) > cf.zFace + 1.e-9) compNo[b.first] += b.second;
        // contornos do fluxo: esqueleto rigido, p_a = 0, carga hidrostatica na lateral direita abaixo de zFace
        contorno = an->AllDisplacements(an->X());
        for (const auto &e : an->AirEq()) contorno[e.second] = 0.;
        for (const auto &b : an->BoundaryNodes(ERight, false)) {
            const REAL z = an->Coord(b.first, 1);
            if (z <= cf.zFace + 1.e-9) {
                const int64_t eq = an->WaterEq().at(b.first);
                contorno[eq] = cf.gw * std::max(cf.nivel - z, REAL(0.));
                eqDir.insert(eq);
                compNo.erase(b.first);
            }
        }
        ColetaPontos();
    }
    ~Encosta() {
        an.reset();
        if (mphys) {
            TPZVec<TPZCompMesh *> mv = mphys->MeshVector();
            delete mphys;
            for (auto *m : mv) delete m;
        }
        delete gmesh;
    }
    Encosta(const Encosta &) = delete;
    Encosta &operator=(const Encosta &) = delete;

    const Camada &CamadaDe(int id) const { return cf.DaCamada(id); }
    const TPZUnsaturatedHydraulics &Hidraulica(int id) const { return mat.at(id)->Hydraulics(); }
    /// p_w do estado convergido num no
    REAL Pw(int64_t no) const { return an->X().GetVal(an->WaterEq().at(no), 0); }
    /// theta e S_e de uma camada com p_w (p_a = 0)
    double ThetaDe(int id, double pw) const {
        REAL S, G;
        Hidraulica(id).Retention(-pw, S, G);
        return CamadaDe(id).thetaS * S;
    }
    double SeDe(int id, double pw) const {
        REAL S, G, Se, dSe;
        Hidraulica(id).Retention(-pw, S, G);
        Hidraulica(id).Effective(S, Se, dSe);
        return Se;
    }

    /// pontos de integracao das camadas: indice de memoria, coordenadas, peso e se o elemento esta na zona
    void ColetaPontos() {
        pontos.clear();
        areaZona = 0.;
        for (TPZCompEl *cel : mphys->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 2) continue;
            TPZGeoEl *gel = cel->Reference();
            const int id = gel->MaterialId();
            if (!mat.count(id)) continue;
            TPZManVector<int64_t, 10> idx;
            cel->GetMemoryIndices(idx);
            const TPZIntPoints &rule = cel->GetIntegrationRule();
            if (idx.size() != rule.NPoints()) throw std::runtime_error("pontos de integracao e memoria diferentes");
            double xc = 0., zc = 0.;
            for (int k = 0; k < 3; ++k) xc += gel->NodePtr(k)->Coord(0) / 3., zc += gel->NodePtr(k)->Coord(1) / 3.;
            const bool z = naZona(xc, zc);
            TPZManVector<REAL, 3> qsi(2, 0.), x(3, 0.);
            TPZFNMatrix<9, REAL> jac, axes, jacinv;
            REAL w, detjac;
            for (int i = 0; i < rule.NPoints(); ++i) {
                rule.Point(i, qsi, w);
                gel->X(qsi, x);
                gel->Jacobian(qsi, jac, axes, detjac, jacinv);
                pontos.push_back({id, idx[i], gel->Index(), x[0], x[1], w * std::fabs(detjac), z});
                if (z) areaZona += w * std::fabs(detjac);
            }
        }
        estado.assign(pontos.size(), EstadoPonto());
    }

    /// carrega um campo de p_w (no -> valor) como estado convergido, com p_a = 0, u = 0 e os contornos do fluxo
    void CarregaEstado(const std::function<double(int64_t)> &pwNo) {
        TPZFMatrix<STATE> x(mphys->NEquations(), 1, 0.);
        for (const auto &e : an->WaterEq()) x(e.second, 0) = pwNo(e.first);
        for (const auto &c : contorno) x(c.first, 0) = c.second;
        an->X() = x;
        an->Restore();
        an->Accept(); // a memoria guarda as pressoes do passo anterior (armazenamento concentrado nos vertices)
        face.clear();
        for (const auto &c : compNo)
            if (Pw(c.first) >= -1.e-9) face.insert(c.first);
    }
    /// condicao inicial (sec. 2.3): hidrostatica com o nivel d'agua em z = nivel
    void CondicaoInicial() {
        CarregaEstado([&](int64_t no) { return cf.gw * (cf.nivel - an->Coord(no, 1)); });
        tempo = 0.;
    }

    /// Um passo de dt com chuva (m/s, por area horizontal) e as faces de percolacao (superficie e lateral direita acima
    /// de zFace): o candidato com p_w > 0 fica com p_w = 0 enquanto a entrada exigida nele nao passar da chuva que ele
    /// recebe; se falhar, volta ao estado do inicio do passo
    bool Passo(REAL dt, REAL chuva) {
        const auto t0 = Relogio::now();
        const REAL tolPw = 1.e-6, tolQ = 1.e-10; // kPa e m/s
        auto tentativa = face;
        for (const auto &c : compNo)
            if (Pw(c.first) > tolPw) tentativa.insert(c.first);
        int newton = 0;
        for (int it = 0; it < 30; it++) {
            auto fixed = contorno;
            std::map<int64_t, REAL> fq;
            for (const auto &c : compNo) {
                const int64_t eq = an->WaterEq().at(c.first);
                if (tentativa.count(c.first)) fixed[eq] = 0.;
            }
            for (const auto &b : chuvaNo) {
                const int64_t eq = an->WaterEq().at(b.first);
                if (!fixed.count(eq) && b.second > 0.) fq[eq] = chuva * b.second; // entrada (m3/s por metro)
            }
            const auto sr = an->Solve(dt, fixed, fq, false);
            newton += sr.iter;
            if (!sr.conv) break;
            const auto R = an->FullResidual(); // R_w/dt: entrada de agua nos nos com p_w imposto
            const TPZFMatrix<STATE> sol = mphys->Solution();
            auto nova = tentativa;
            for (const auto &c : compNo) {
                const int64_t eq = an->WaterEq().at(c.first);
                const REAL qno = chuvaNo.count(c.first) ? chuva * chuvaNo.at(c.first) : 0.;
                if (tentativa.count(c.first)) {
                    if (R.GetVal(eq, 0) / dt > qno + tolQ * c.second) nova.erase(c.first);
                } else if (sol.GetVal(eq, 0) > tolPw) nova.insert(c.first);
            }
            if (nova != tentativa) {
                tentativa = nova;
                continue;
            }
            InfoPasso info;
            info.dt = dt, info.newton = newton, info.itsFace = it + 1;
            double somaR = 0.;
            info.pwMin = 1.e300, info.pwMax = -1.e300;
            for (const auto &e : an->WaterEq()) {
                somaR += R.GetVal(e.second, 0);
                const REAL pw = sol.GetVal(e.second, 0);
                info.pwMin = std::min<double>(info.pwMin, pw);
                info.pwMax = std::max<double>(info.pwMax, pw);
                info.dpwMax = std::max<double>(info.dpwMax, std::fabs(pw - an->X().GetVal(e.second, 0)));
            }
            for (const auto &b : chuvaNo)
                if (!tentativa.count(b.first)) info.chuva += chuva * b.second * dt;
            for (int64_t no : tentativa) info.face += R.GetVal(an->WaterEq().at(no), 0);
            for (int64_t eq : eqDir) info.dir += R.GetVal(eq, 0);
            info.dV = somaR;
            info.residuo = somaR - info.chuva - info.face - info.dir;
            an->Accept();
            face = tentativa;
            info.faceNos = (int)face.size();
            info.seg = Segundos(t0);
            ultimo = info;
            passos++;
            newtonTotal += newton;
            return true;
        }
        an->Restore();
        newtonTotal += newton;
        return false;
    }

    /// avanca de dt dividindo ao meio o passo que nao converge (ate maxCortes niveis)
    bool Avanca(REAL dt, REAL chuva, std::ostream *fpas, int etapa, int nivel = 0) {
        if (Passo(dt, chuva)) {
            tempo += dt;
            Registra(fpas, etapa);
            return true;
        }
        if (nivel >= cf.maxCortes) return false;
        cortes++;
        return Avanca(0.5 * dt, chuva, fpas, etapa, nivel + 1) && Avanca(0.5 * dt, chuva, fpas, etapa, nivel + 1);
    }
    static void CabecalhoPassos(std::ostream &f) {
        f << "etapa,passo,t_h,t_dias,dt_h,newton,its_face,face_nos,chuva_m3,face_m3,dir_m3,dV_m3,residuo_m3,dpw_max_kPa,"
             "pw_min_kPa,pw_max_kPa,seg\n";
    }
    void Registra(std::ostream *fpas, int etapa) const {
        if (Verboso() >= 2)
            std::cout << "  etapa " << etapa << " passo " << passos << ": t = " << tempo / kHora << " h (" << tempo / kDia
                      << " d), dt = " << ultimo.dt / kHora << " h, newton = " << ultimo.newton << ", face = "
                      << ultimo.faceNos << ", max|dp_w| = " << ultimo.dpwMax << " kPa" << std::endl;
        if (!fpas) return;
        const InfoPasso &i = ultimo;
        *fpas << etapa << ',' << passos << ',' << tempo / kHora << ',' << tempo / kDia << ',' << i.dt / kHora << ','
              << i.newton << ',' << i.itsFace << ',' << i.faceNos << ',' << i.chuva << ',' << i.face << ',' << i.dir << ','
              << i.dV << ',' << i.residuo << ',' << i.dpwMax << ',' << i.pwMin << ',' << i.pwMax << ',' << i.seg << '\n';
    }

    /// Etapa 1: ANOS anos de chuva CHUVA1_MM_ANO a partir do estado carregado, com dt crescente (de DT0_H, vezes
    /// FATOR_DT se o passo levou ate 8 iteracoes de Newton, ate DT_MAX_D; vezes 0,7 se levou mais de 15)
    bool Etapa1(std::ostream *fpas = nullptr) {
        const REAL T = cf.anos * kAno, q = cf.Chuva1();
        REAL dt = cf.dt0_h * kHora;
        tempo = 0.;
        while (tempo < T * (1. - 1.e-12)) {
            const REAL d = std::min(dt, T - tempo);
            if (!Avanca(d, q, fpas, 1)) return false;
            const int nw = ultimo.newton;
            dt = nw <= 8 ? dt * cf.fatorDt : nw > 15 ? dt * 0.7 : dt;
            dt = std::min(dt, REAL(cf.dtMax_d * kDia));
        }
        return true;
    }

    /// Etapa 2: chuva CHUVA2_MM_H com passo DT2_H ate cada instante de SAIDAS_H (o ultimo passo de cada intervalo e
    /// encurtado); saida(t_h) e chamada em cada instante (t_h = 0: o fim da etapa 1)
    bool Etapa2(const std::function<void(double)> &saida, std::ostream *fpas = nullptr) {
        const REAL q = cf.Chuva2(), dt = cf.dt2_h * kHora;
        tempo = 0.;
        for (double th : cf.saidasH) {
            const REAL alvo = th * kHora;
            while (tempo < alvo * (1. - 1.e-12) - 1.e-9) {
                const REAL d = std::min(dt, alvo - tempo);
                if (!Avanca(d, q, fpas, 2)) return false;
            }
            if (saida) saida(th);
        }
        return true;
    }

    /// Tensoes do estado convergido: equilibrio elastico com as pressoes congeladas (so o peso, que depende da
    /// umidade), a partir do estado sem tensoes; guarda a tensao total, p_w e S de cada ponto em estado e o LFS com os
    /// c' e phi' do Config; o estado do fluxo (solucao e memoria) volta ao que era
    bool Tensoes() {
        const TPZFMatrix<STATE> X0 = an->X();
        const auto mem0 = an->SaveMemory();
        an->ResetMemory();
        TPZFMatrix<STATE> &X = an->X();
        for (int64_t i = 0; i < an->NEq(); ++i)
            if (!an->IsPressure(i)) X(i, 0) = 0.;
        an->Restore();
        const auto fixed = an->FrozenFixed(X);
        an->SetGravityFactor(1.);
        an->SetStrengthReductionFactor(1.);
        an->SetGravityScalesSuction(true);
        const auto r = an->Solve(0., fixed, {}, true);
        bool ok = r.conv;
        if (ok) {
            an->Accept();
            xMec = an->X();
            for (size_t i = 0; i < pontos.size(); ++i) {
                const TUnsat *m = mat.at(pontos[i].camada);
                const TPZUnsaturatedMem &mm = (*const_cast<TUnsat *>(m)->GetMemory())[pontos[i].mem];
                EstadoPonto &e = estado[i];
                e.sxx = mm.m_sigma.XX(), e.syy = mm.m_sigma.YY(), e.sxy = mm.m_sigma.XY(), e.szz = mm.m_sigma.ZZ();
                e.pw = mm.fPw, e.S = mm.fS;
                REAL Se, dSe;
                m->Hydraulics().Effective(mm.fS, Se, dSe);
                e.Se = Se;
            }
            CalculaLFS();
        } else
            std::cerr << "Tensoes: o equilibrio elastico nao convergiu (" << r.iter << " iteracoes)" << std::endl;
        an->RestoreMemory(mem0);
        an->X() = X0;
        an->Restore();
        for (auto *m : an->Materials()) m->SetFrozenPressures(false);
        return ok;
    }

    /// LFS de cada ponto com o estado mecanico guardado e os c' (por camada) e phi' do Config
    void CalculaLFS() { CalculaLFS(cf.sup.c, cf.inf.c, cf.Phi()); }
    /// LFS de cada ponto com c' de cada camada e phi' (rad) dados
    void CalculaLFS(double cSup, double cInf, double phi) {
        for (size_t i = 0; i < pontos.size(); ++i) {
            EstadoPonto &e = estado[i];
            e.lfs = LFSdeTensao(e.sxx, e.syy, e.sxy, e.pw, e.Se, pontos[i].camada == kSuperior ? cSup : cInf, phi);
        }
    }
    /// troca o coeficiente de Poisson do modelo elastico das duas camadas (o proximo calculo de tensoes o usa; o
    /// fluxo nao depende dele)
    void TrocaNu(double nu) {
        for (auto &km : mat) {
            UnsatParams p = ParametrosCamada(cf, CamadaDe(km.first));
            p.nu = nu;
            km.second->SetPlasticModel(MohrCoulomb(p));
        }
    }

    /// media do LFS e area com LFS < 1 na zona (eqs. 11 a 13), com os pontos de integracao (os elementos estao todos
    /// dentro ou fora da zona: a malha tem as arestas da zona); lfsMedioLim e lfsMedioLim3: a media com o LFS
    /// limitado a 10 e a 3 (o topo da escala de cores das Figs. 6 e 9 do artigo; ver README.md)
    Metricas CalculaMetricas(double t_h) const {
        Metricas m;
        m.t_h = t_h;
        double s = 0., sl = 0., s3 = 0., a = 0.;
        m.lfsMin = 1.e300;
        m.thetaSupMin = m.thetaInfMin = 1.e300;
        m.thetaSupMax = m.thetaInfMax = -1.e300;
        for (size_t i = 0; i < pontos.size(); ++i) {
            const Ponto &p = pontos[i];
            if (!p.zona) continue;
            const EstadoPonto &e = estado[i];
            a += p.w;
            s += p.w * e.lfs;
            sl += p.w * std::min(e.lfs, 10.);
            s3 += p.w * std::min(e.lfs, 3.);
            m.lfsMin = std::min(m.lfsMin, e.lfs);
            if (e.lfs < 1.) {
                m.aLFS += p.w;
                (p.camada == kSuperior ? m.aLFSsup : m.aLFSinf) += p.w;
            }
            if (e.lfs > 3.) m.aLFS3 += p.w;
            const double th = CamadaDe(p.camada).thetaS * e.S;
            if (p.camada == kSuperior) m.thetaSupMin = std::min(m.thetaSupMin, th), m.thetaSupMax = std::max(m.thetaSupMax, th);
            else m.thetaInfMin = std::min(m.thetaInfMin, th), m.thetaInfMax = std::max(m.thetaInfMax, th);
            if (e.pw < 0.) m.sigsMax = std::max(m.sigsMax, -e.Se * e.pw);
        }
        m.areaZona = a;
        m.lfsMedio = a > 0. ? s / a : 0.;
        m.lfsMedioLim = a > 0. ? sl / a : 0.;
        m.lfsMedioLim3 = a > 0. ? s3 / a : 0.;
        return m;
    }

    // ------------------------------------------------------------------ amostragem (perfis)
    /// elemento 2D que contem (x, z) e as coordenadas baricentricas; -1 se fora
    int64_t Localiza(double x, double z, double l[3]) const {
        for (int64_t i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *g = gmesh->Element(i);
            if (!g || g->Dimension() != 2 || g->HasSubElement()) continue;
            const double x0 = g->NodePtr(0)->Coord(0), z0 = g->NodePtr(0)->Coord(1);
            const double x1 = g->NodePtr(1)->Coord(0), z1 = g->NodePtr(1)->Coord(1);
            const double x2 = g->NodePtr(2)->Coord(0), z2 = g->NodePtr(2)->Coord(1);
            const double det = (x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0);
            const double l1 = ((x - x0) * (z2 - z0) - (x2 - x0) * (z - z0)) / det;
            const double l2 = ((x1 - x0) * (z - z0) - (x - x0) * (z1 - z0)) / det;
            const double l0 = 1. - l1 - l2, tol = -1.e-9;
            if (l0 >= tol && l1 >= tol && l2 >= tol) {
                l[0] = l0, l[1] = l1, l[2] = l2;
                return i;
            }
        }
        return -1;
    }
    /// valor de um ponto amostrado: p_w (P1), theta, S_e e tensao de succao S_e |p_w| da camada do elemento, e o LFS
    /// medio dos pontos de integracao do elemento (do ultimo calculo de tensoes)
    struct Amostra { int64_t el = -1; int camada = 0; double pw = 0., theta = 0., Se = 0., sigs = 0., lfs = 0.; };
    std::map<int64_t, std::pair<double, int>> lfsEl; ///< LFS medio e numero de pontos por elemento
    void AtualizaLFSElementos() {
        lfsEl.clear();
        for (size_t i = 0; i < pontos.size(); ++i) {
            auto &v = lfsEl[pontos[i].el];
            v.first += estado[i].lfs;
            v.second++;
        }
    }
    Amostra Amostrar(double x, double z) const {
        Amostra a;
        double l[3];
        a.el = Localiza(x, z, l);
        if (a.el < 0) return a;
        TPZGeoEl *g = gmesh->Element(a.el);
        a.camada = g->MaterialId();
        for (int k = 0; k < 3; ++k) a.pw += l[k] * Pw(g->NodeIndex(k));
        a.theta = ThetaDe(a.camada, a.pw);
        a.Se = SeDe(a.camada, a.pw);
        a.sigs = a.pw < 0. ? -a.Se * a.pw : 0.;
        auto it = lfsEl.find(a.el);
        if (it != lfsEl.end() && it->second.second > 0) a.lfs = it->second.first / it->second.second;
        return a;
    }
    /// cota do nivel d'agua numa vertical x: subindo da base, a primeira troca de p_w >= 0 para p_w < 0
    double NivelAgua(double x, double dz = 0.02) const {
        const double zt = CotaTerreno(x);
        double zant = 0., pant = Amostrar(x, 1.e-6).pw;
        for (double z = dz; z <= zt + 1.e-9; z += dz) {
            const double p = Amostrar(x, std::min(z, zt - 1.e-6)).pw;
            if (pant >= 0. && p < 0.) return zant + (z - zant) * pant / (pant - p);
            zant = z, pant = p;
        }
        return pant >= 0. ? zt : -1.;
    }
    /// profundidade da frente de umedecimento numa vertical x (camada superior): descendo da superficie, a primeira
    /// profundidade com theta < thetaFrente; a espessura da camada se ela toda passa de thetaFrente
    double Frente(double x, double thetaFrente = 0.32, double dz = 0.01) const {
        const double zt = CotaTerreno(x);
        double dant = 0., tant = Amostrar(x, zt - 1.e-6).theta;
        for (double d = dz; d <= kEspSup + 1.e-9; d += dz) {
            const double th = Amostrar(x, zt - std::min(d, kEspSup - 1.e-6)).theta;
            if (tant >= thetaFrente && th < thetaFrente) return dant + (d - dant) * (tant - thetaFrente) / (tant - th);
            dant = d, tant = th;
        }
        return tant >= thetaFrente ? kEspSup : 0.;
    }

    // ------------------------------------------------------------------ arquivos
    /// estado (p_w nos nos): no, x, z, pw_kPa
    void GravaEstado(const std::string &arq) const {
        std::ofstream f(arq);
        f << "no,x,z,pw_kPa\n" << std::setprecision(15);
        for (const auto &e : an->WaterEq())
            f << e.first << ',' << an->Coord(e.first, 0) << ',' << an->Coord(e.first, 1) << ',' << Pw(e.first) << '\n';
    }
    void LeEstado(const std::string &arq) {
        std::ifstream f(arq);
        if (!f) throw std::runtime_error("estado nao encontrado: " + arq);
        std::string l;
        std::getline(f, l);
        std::map<int64_t, double> pw;
        while (std::getline(f, l)) {
            std::stringstream ss(l);
            std::string c[4];
            for (auto &s : c) std::getline(ss, s, ',');
            const int64_t no = std::atoll(c[0].c_str());
            if (!an->WaterEq().count(no) || std::fabs(an->Coord(no, 0) - std::atof(c[1].c_str())) > 1.e-6 ||
                std::fabs(an->Coord(no, 1) - std::atof(c[2].c_str())) > 1.e-6)
                throw std::runtime_error("estado de outra malha: " + arq);
            pw[no] = std::atof(c[3].c_str());
        }
        if (pw.size() != an->WaterEq().size()) throw std::runtime_error("estado incompleto: " + arq);
        CarregaEstado([&](int64_t no) { return pw.at(no); });
    }

    /// malha (uma vez): nos.csv (no, x, z) e elementos.csv (el, camada, n0, n1, n2, zona)
    void GravaMalha(const std::string &pasta = ".") const {
        std::ofstream fn(pasta + "/nos.csv"), fe(pasta + "/elementos.csv");
        fn << "no,x,z\n" << std::setprecision(12);
        for (int64_t i = 0; i < gmesh->NNodes(); ++i)
            fn << i << ',' << gmesh->NodeVec()[i].Coord(0) << ',' << gmesh->NodeVec()[i].Coord(1) << '\n';
        fe << "el,camada,n0,n1,n2,zona\n";
        for (int64_t i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *g = gmesh->Element(i);
            if (!g || g->Dimension() != 2) continue;
            double xc = 0., zc = 0.;
            for (int k = 0; k < 3; ++k) xc += g->NodePtr(k)->Coord(0) / 3., zc += g->NodePtr(k)->Coord(1) / 3.;
            fe << i << ',' << g->MaterialId() << ',' << g->NodeIndex(0) << ',' << g->NodeIndex(1) << ','
               << g->NodeIndex(2) << ',' << (naZona(xc, zc) ? 1 : 0) << '\n';
        }
    }
    /// campos de um instante por elemento (para os mapas): theta e tensao de succao nos tres vertices (com a camada do
    /// elemento), LFS medio e minimo dos pontos, tensoes totais medias
    void GravaCampos(const std::string &arq) const {
        std::map<int64_t, std::array<double, 6>> el; // LFS soma, min, sxx, syy, sxy, n
        for (size_t i = 0; i < pontos.size(); ++i) {
            auto it = el.find(pontos[i].el);
            if (it == el.end()) it = el.emplace(pontos[i].el, std::array<double, 6>{0., 1.e300, 0., 0., 0., 0.}).first;
            auto &v = it->second;
            const EstadoPonto &e = estado[i];
            v[0] += e.lfs, v[1] = std::min(v[1], e.lfs), v[2] += e.sxx, v[3] += e.syy, v[4] += e.sxy, v[5] += 1.;
        }
        std::ofstream f(arq);
        f << "el,camada,pw0,pw1,pw2,theta0,theta1,theta2,sigs0,sigs1,sigs2,lfs,lfs_min,sxx,syy,sxy\n" << std::setprecision(8);
        for (const auto &kv : el) {
            TPZGeoEl *g = gmesh->Element(kv.first);
            const int id = g->MaterialId();
            const auto &v = kv.second;
            double pw[3], th[3], sg[3];
            for (int k = 0; k < 3; ++k) {
                pw[k] = Pw(g->NodeIndex(k));
                th[k] = ThetaDe(id, pw[k]);
                sg[k] = pw[k] < 0. ? -SeDe(id, pw[k]) * pw[k] : 0.;
            }
            f << kv.first << ',' << id << ',' << pw[0] << ',' << pw[1] << ',' << pw[2] << ',' << th[0] << ',' << th[1] << ','
              << th[2] << ',' << sg[0] << ',' << sg[1] << ',' << sg[2] << ',' << v[0] / v[5] << ',' << v[1] << ','
              << v[2] / v[5] << ',' << v[3] / v[5] << ',' << v[4] / v[5] << '\n';
        }
    }
    /// pontos de integracao de um instante
    void GravaPontos(const std::string &arq) const {
        std::ofstream f(arq);
        f << "x,z,w,camada,zona,pw_kPa,S,Se,theta,sigs_kPa,sxx,syy,sxy,szz,lfs\n" << std::setprecision(8);
        for (size_t i = 0; i < pontos.size(); ++i) {
            const Ponto &p = pontos[i];
            const EstadoPonto &e = estado[i];
            f << p.x << ',' << p.z << ',' << p.w << ',' << p.camada << ',' << p.zona << ',' << e.pw << ',' << e.S << ','
              << e.Se << ',' << CamadaDe(p.camada).thetaS * e.S << ',' << (e.pw < 0. ? -e.Se * e.pw : 0.) << ',' << e.sxx
              << ',' << e.syy << ',' << e.sxy << ',' << e.szz << ',' << e.lfs << '\n';
        }
    }
    /// VTK legado: p_w, theta (media das camadas nos nos da interface) e deslocamento nos nos; LFS, theta e camada
    /// por elemento
    void GravaVTK(const std::string &arq) const {
        std::vector<int64_t> els;
        for (int64_t i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *g = gmesh->Element(i);
            if (g && g->Dimension() == 2) els.push_back(i);
        }
        std::map<int64_t, std::array<double, 2>> lfs; // soma, n
        for (size_t i = 0; i < pontos.size(); ++i) lfs[pontos[i].el][0] += estado[i].lfs, lfs[pontos[i].el][1] += 1.;
        std::vector<double> thNo(gmesh->NNodes(), 0.), nNo(gmesh->NNodes(), 0.);
        for (int64_t i : els) {
            TPZGeoEl *g = gmesh->Element(i);
            for (int k = 0; k < 3; ++k) {
                const int64_t n = g->NodeIndex(k);
                thNo[n] += ThetaDe(g->MaterialId(), Pw(n));
                nNo[n] += 1.;
            }
        }
        std::ofstream v(arq);
        v << "# vtk DataFile Version 3.0\nV5_EncostaDuasCamadas\nASCII\nDATASET UNSTRUCTURED_GRID\n" << std::setprecision(8);
        v << "POINTS " << gmesh->NNodes() << " double\n";
        for (int64_t i = 0; i < gmesh->NNodes(); ++i)
            v << gmesh->NodeVec()[i].Coord(0) << ' ' << gmesh->NodeVec()[i].Coord(1) << " 0\n";
        v << "CELLS " << els.size() << ' ' << 4 * els.size() << '\n';
        for (int64_t i : els) {
            TPZGeoEl *g = gmesh->Element(i);
            v << "3 " << g->NodeIndex(0) << ' ' << g->NodeIndex(1) << ' ' << g->NodeIndex(2) << '\n';
        }
        v << "CELL_TYPES " << els.size() << '\n';
        for (size_t i = 0; i < els.size(); ++i) v << "5\n";
        v << "POINT_DATA " << gmesh->NNodes() << "\nSCALARS pw_kPa double 1\nLOOKUP_TABLE default\n";
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) v << (an->WaterEq().count(i) ? Pw(i) : 0.) << '\n';
        v << "SCALARS theta double 1\nLOOKUP_TABLE default\n";
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) v << (nNo[i] > 0. ? thNo[i] / nNo[i] : 0.) << '\n';
        v << "VECTORS deslocamento_m double\n";
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) {
            const bool tem = an->DispEq().count(i) && xMec.Rows() == an->NEq();
            v << (tem ? an->NodalDisplacement(xMec, i, 0) : 0.) << ' ' << (tem ? an->NodalDisplacement(xMec, i, 1) : 0.)
              << " 0\n";
        }
        v << "CELL_DATA " << els.size() << "\nSCALARS LFS double 1\nLOOKUP_TABLE default\n";
        for (int64_t i : els) {
            auto it = lfs.find(i);
            v << (it != lfs.end() && it->second[1] > 0. ? it->second[0] / it->second[1] : 0.) << '\n';
        }
        v << "SCALARS camada int 1\nLOOKUP_TABLE default\n";
        for (int64_t i : els) v << gmesh->Element(i)->MaterialId() << '\n';
    }
};

// ------------------------------------------------------------------------------------------------------------------
/// secoes verticais dos perfis (x absoluto): crista (x_zona = 2 m), meio do talude (10 m), pe do talude (17 m) e
/// patamar do pe (20,5 m); x_zona = x - 15
inline std::vector<double> SecoesPerfis() { return {17., 25., 32., 35.5}; }

/// resultados de um instante usados no resumo e nos testes
struct Instante {
    Metricas m;
    double nivelCrista = 0., frenteCrista = 0., frenteMeio = 0.;
    double thetaSupCrista = 0., thetaInfCrista = 0., sigsSupCrista = 0., sigsInfCrista = 0.;
    double thetaPe = 0.; ///< theta medio da camada superior no patamar do pe (x = 35,5 m)
};

/// Os valores de um instante: metricas do LFS, nivel d'agua e frente na crista e no meio do talude, theta e tensao de
/// succao no meio da camada superior (1 m de profundidade) e a 1 m abaixo da interface na crista (x = 17 m)
inline Instante AvaliaInstante(Encosta &e, double t_h) {
    Instante r;
    r.m = e.CalculaMetricas(t_h);
    e.AtualizaLFSElementos();
    r.nivelCrista = e.NivelAgua(17.);
    r.frenteCrista = e.Frente(17.);
    r.frenteMeio = e.Frente(25.);
    const auto s1 = e.Amostrar(17., kZCrista - 1.), s2 = e.Amostrar(17., kZCrista - kEspSup - 1.);
    r.thetaSupCrista = s1.theta, r.sigsSupCrista = s1.sigs;
    r.thetaInfCrista = s2.theta, r.sigsInfCrista = s2.sigs;
    double s = 0.;
    int n = 0;
    for (double z = kZPe - 0.05; z > kZPe - kEspSup; z -= 0.1, ++n) s += e.Amostrar(35.5, z).theta;
    r.thetaPe = n ? s / n : 0.;
    return r;
}

/// Sensibilidade das metricas do LFS num instante, com o mesmo fluxo: para cada nu de nus, um novo calculo de tensoes;
/// para cada c' (o mesmo nas duas camadas) de cs e cada phi' (graus) de phis, so o LFS com as tensoes do nu do caso.
/// Grava em f as linhas t_h, parametro, valor, LFS_medio, LFS_medio_lim3, A_LFS_m2; no fim o estado mecanico e o LFS
/// voltam aos do caso. Retorna false se algum calculo de tensoes falhou.
inline bool Sensibilidade(Encosta &e, double t_h, const std::vector<double> &nus, const std::vector<double> &cs,
                          const std::vector<double> &phis, std::ostream &f) {
    bool ok = true;
    auto linha = [&](const char *par, double v) {
        const Metricas m = e.CalculaMetricas(t_h);
        f << t_h << ',' << par << ',' << v << ',' << m.lfsMedio << ',' << m.lfsMedioLim3 << ',' << m.aLFS << '\n';
    };
    for (double nu : nus) {
        e.TrocaNu(nu);
        if (e.Tensoes()) linha("nu", nu);
        else ok = false;
    }
    if (!nus.empty()) {
        e.TrocaNu(e.cf.nu);
        ok = e.Tensoes() && ok;
    }
    for (double c : cs) {
        e.CalculaLFS(c, c, e.cf.Phi());
        linha("c", c);
    }
    for (double p : phis) {
        e.CalculaLFS(e.cf.sup.c, e.cf.inf.c, p * M_PI / 180.);
        linha("phi", p);
    }
    e.CalculaLFS();
    f.flush();
    return ok;
}

/// perfis de um instante nas secoes de SecoesPerfis (passo de 0,05 m da superficie ate a base da zona, z = 17 m)
inline void GravaPerfis(Encosta &e, double t_h, std::ostream &f, bool cabecalho) {
    if (cabecalho) f << "t_h,x,x_zona,z,z_zona,prof,camada,pw_kPa,theta,Se,sigs_kPa,lfs\n";
    e.AtualizaLFSElementos();
    f << std::setprecision(8);
    for (double x : SecoesPerfis()) {
        const double zt = CotaTerreno(x);
        for (double d = 0.; zt - d >= kZZona0 - 1.e-9; d += 0.05) {
            const double z = std::max(std::min(zt - d, zt - 1.e-6), kZZona0 + 1.e-6);
            const auto a = e.Amostrar(x, z);
            if (a.el < 0) continue;
            f << t_h << ',' << x << ',' << x - kXZona0 << ',' << z << ',' << z - kZZona0 << ',' << d << ',' << a.camada << ','
              << a.pw << ',' << a.theta << ',' << a.Se << ',' << a.sigs << ',' << a.lfs << '\n';
        }
    }
}

} // namespace v5
