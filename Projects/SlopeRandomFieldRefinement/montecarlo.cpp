// Monte Carlo de um talude sob chuva com campos aleatorios lognormais de c', phi' e k_s. Padrao: o talude de Rahardjo
// et al. (2007) (RandomField: Karhunen-Loeve analitico; campos da Tabela 2 de Vargas Ceron et al. 2025), talude da
// Fig. 2 do Rahardjo (H = 10 m a 45 graus) na TriGMesh do main.cpp. As opcoes da segunda lista abaixo trocam o problema
// (ex.: o caso A, Hu et al. 2024, nos cluster/casos_hu2024.txt). Em cada realizacao:
//  1. transiente da chuva (RainTransient: esqueleto rigido, ar passivo) na malha do transiente, com o k de cada ponto do
//     campo, ate o ultimo instante pedido; guarda p_w em cada instante do fator. Com k fixo (covk = 0) o transiente e o
//     mesmo em todas as realizacoes: e calculado uma vez por processo (ou lido de um arquivo, pressoes=);
//  2. em cada instante, o fator de seguranca (SRM: c' e tan(phi') de cada ponto divididos por F; ou GIM) com as pressoes
//     levadas para a malha do fator (PointLocator + SetFrozenPressures), com refinamento adaptativo em h feito do zero
//     (instantes de 'adapta') ou na malha refinada do instante adaptativo mais proximo (os outros);
//  3. descritores do mecanismo (a faixa dos pontos com sqrt(J2(eps_p)) >= 0,1 do maior) e do estado hidraulico.
// Uso: SlopeRainMonteCarlo [chave=valor ...], com (padrao entre parenteses):
//   n (10), primeira (0)   realizacoes primeira..primeira+n-1
//   det (0)                1: antes, o caso deterministico (as medias, sem campos), como a realizacao -1
//   valores ()             c',phi',k: a deterministica com esses valores em todo o solo (c' em kPa, phi' em graus,
//                          k / k_s), em vez das medias; ex.: valores=10,26,0.5 (o FS de um ponto do FORM)
//   solo (f50)             f10, f50 ou f100 (Tabela 1 do Rahardjo)
//   I (1 k_s do solo)      intensidade da chuva (mm/h); chuva (24): duracao da chuva (h)
//   hw (5)                 profundidade do nivel d'agua inicial abaixo do pe (m)
//   tempos (0,24)          instantes do fator (h, separados por virgula); o transiente vai ate o ultimo
//   adapta (os tempos)     instantes com refinamento do zero; os outros usam a malha do adaptativo mais proximo
//   metodo (srm)           srm ou gim (aumento do peso; a succao fica fixa)
//   srm_f0 (0.8)           F do peso no SRM, o inicio do caminho de F (1: o SRM antigo, com o buraco logo abaixo de FS = 1)
//   covc (0.3), covphi (0.1), covk (0.6)   CoV de c', phi' e k_s
//   s (1)                  escala das distancias de autocorrelacao: 20 s (horizontal) e 2 s (vertical), em m
//   alfa (1)               anisotropia k_h / k_v (k_v e o do campo)
//   transiente (tri2)      malha do transiente: triK, camH ou camqH (camHpD, camqHpD: camadas de H ate D m); uT (2):
//                          ordem de u no transiente
//   fator (tri1), nref (4) malha grossa do fator (triK ou camH) e niveis de refinamento
//   figuras ()             realizacoes com os campos gravados para figuras (pasta <saida>_real<n>; nelas todos os
//                          instantes sao adaptativos), ex.: 0,3,7
//   saida (mc_chuva)       prefixo dos arquivos
//   semente (2025)         semente dos campos (outra para uma amostra independente de outra rodada)
// O problema (padrao: o Rahardjo de antes, com os mesmos numeros):
//   geometria ()           os 6 pontos do talude, x,y alternados (base esquerda, base direita, terreno do pe, pe, crista,
//                          crista a esquerda); padrao: o do Rahardjo, 0,0,70,0,70,20,40,20,30,30,0,30. Caso A (Hu et al.
//                          2024): 0,0,70,0,70,10,40,10,30,20,0,20
//   na ()                  cota do nivel d'agua inicial no pe (m), no lugar de hw (hw = cota do pe - na)
//   na_beta (7)            inclinacao do nivel inicial (graus, subindo para dentro do macico; 0: horizontal)
//   succao_max (75)        teto da succao inicial (kPa): p_w = max(gamma_w (nivel - y), -succao_max)
//   esquerda (carga), direita (carga)   carga: p_w hidrostatica abaixo do nivel, sem fluxo acima; fechada: sem fluxo
//                          em toda a lateral (a base e sempre fechada)
//   carga_esq (), carga_dir ()          carga total (m, a cota do nivel) da lateral com carga, no lugar do nivel inicial
//   chuva_face (inclinada) inclinada: q = I por metro inclinado da superficie (normal a face, o "unit flux" do
//                          GeoStudio); vertical: I por metro de projecao horizontal (chuva vertical, a decisao D1)
//   ks ()                  k_s (m/s) no lugar do k_s do solo
//   vg ()                  theta_s,theta_r,alfa,n: retencao de van Genuchten (alfa em 1/kPa) com k_r de Mualem, no lugar
//                          de Fredlund-Xing com k_s Theta^4 (porosidade = theta_s; S_wr = theta_r / theta_s)
//   vg_tab (0)             eps[,fecho]: com vg=, o k_r de Mualem tabelado contra S_w (MualemTabelado.h) ate o teto S_e =
//                          1 - eps (1e-5 <= eps <= 0,1) e fechado ate k_r(1) = 1 por reta (padrao), c1 ou ippisch, no lugar
//                          do EMualem (com n < 2, dk_r/dpsi explode na saturacao e o Newton do transiente cicla); 0: o
//                          EMualem (desligado)
//   vg_tab_psi ()          psi_c[,fecho]: o mesmo, com eps = 1 - S_e(psi_c) (psi_c em kPa; ippisch: a entrada de ar psi_e;
//                          ex.: vg_tab_psi=0.01,ippisch). Um ou outro
//   chi (const)            const: chi = tan(phi^b) / tan(phi') com phi^b = 26 graus (o Rahardjo); se: chi = S_e
//                          (Vanapalli et al. 1996)
//   c (10), phi (26), gamma (20), E (20000), nu (0.3)   medias de c' (kPa) e phi' (graus), peso especifico (kN/m3), E (kPa)
//   rcphi (0)              correlacao cruzada dos campos gaussianos de c' e phi' (como no Cho 2010: -0,5 no caso A)
//   fx ()                  a,n,m,theta_s: Fredlund-Xing (a em kPa, C(psi) = 1) no lugar do solo (porosidade = theta_s);
//                          caso B (Liu et al. 2024): fx=5,2,1,0.4
//   kfredlund ()           psi_aev (kPa): k_r pela integral de Fredlund et al. (1994) (Eq. 9 de Liu et al. 2024), com a
//                          retencao com C(psi) e psi_r (psi_r=, 3000 kPa), tabelada contra S_w (TPZUnsaturatedHydraulics::
//                          ETabulated), no lugar de k_s Theta^4; k_r = min(1, integral / a do denominador desde psi_aev)
//   q0 (0)                 infiltracao do estado inicial (mm/h, na convencao de chuva_face=): o estado inicial passa a ser
//                          a percolacao permanente com esse fluxo na superficie (RainTransient::SteadyState), a partir do
//                          nivel inicial; a chuva comeca depois (t = 0). Caso B: q0=0.36 (1e-7 m/s), I=0, tempos=0
//   q0_cov (0)             q0 normal com esse CoV, sorteado por realizacao (caso C, Liu e Wang 2023: a chuva anual de
//                          1200 mm com CoV 0,1, q0 = 1200/8760 mm/h); truncado em 0
//   chi (const)            tambem ponto: chi = tan(phi^b(x)) / tan(phi'(x)) em cada ponto, com phi^b = phib= (graus) ou
//                          phib_razao= vezes phi'(x) (posto nas pressoes congeladas do fator: p_w chi(x) onde p_w < 0)
//   phib (26), phib_razao ()   phi^b (graus) de chi=const e de chi=ponto; ou phi^b = razao phi' (chi=ponto; caso B: 2/3)
//   dt (0.1)               passo de tempo do transiente (h)
//   dt_min (0)             > 0: um passo do transiente que nao converge e refeito em duas metades, ate dt_min (h); 0: o
//                          transiente para (o de antes)
//   newton ()              tol_p,maxit[,busca]: o Newton do transiente para com |dp| < tol_p (kPa; padrao do
//                          TPZUnsatAnalysis 1e-6) e ate maxit iteracoes (40); busca: metades da busca linear (4)
//   newton_trace (<arquivo>)  diagnostico de um transiente que para: o ultimo Solve que falhou e refeito com maxit - 2,
//                          maxit - 1 e maxit iteracoes, e o arquivo ganha um bloco com os 10 nos de maior mudanca de p_w na
//                          ultima iteracao (no, x, y, os tres p_w, superficie, face; RainTransient::SetNewtonTrace); a tela
//                          ganha uma linha por instante de tempos= atingido (passos, iteracoes, nos na face)
//   face_tol (0)           um no da face de percolacao so volta ao fluxo da chuva se R / (dt L) > q (1 + face_tol) + 1e-10
//                          (0: o criterio de antes)
// A chuva aleatoria (caso C):
//   chuva_aleatoria (0)    1: I e D sorteados por realizacao (Pareto generalizada com copula de Frank, Liu e Wang 2023,
//                          Tabela 2), com os instantes do fator 0 e D (o FS inicial e o do fim da chuva; adapta nos dois);
//                          <saida>_chuva_<faixa>.csv tem I, D e q0 de cada realizacao. A deterministica usa I= e chuva=
//   gpd (-0.04,3.44,5.13,0.43,5.15,0.30,-2.38)  xi_I,mu_I,sigma_I,xi_D,mu_D,sigma_D,omega (mm/h, h): as marginais e a copula;
//                          o artigo imprime sigma_D = -0,30 h, impossivel (a escala e positiva): o padrao usa |sigma_D|
//   crpc ()                I1,I2,... (mm/h, decrescentes): no lugar do Monte Carlo, a curva critica (CRPC, Liu e Wang 2023)
//                          de cada realizacao: para cada I, chuva constante ate crpc_dmax h; o FS em 6, 12, 24, 48, 96, 150
//                          e crpc_dmax h ate o primeiro abaixo de 1 e, entre os dois ultimos, bissecao nos instantes
//                          guardados de crpc_passo em crpc_passo h; D_C por interpolacao linear do FS. Se uma I nao derruba
//                          ate crpc_dmax, as menores nao sao feitas; FS(0) < 1: realizacao descartada (nada mais).
//                          Saida: <saida>_crpc_<faixa>.csv (cluster/crpc_liuwang.py integra P_1)
//   crpc_dmax (200), crpc_passo (1)
// Os campos:
//   M (2000)               termos da KL por campo (com nystrom=: os M primeiros da base, M <= os do arquivo)
//   nystrom (<arquivo>)    a base de Karhunen-Loeve de Nystrom do proprio talude (Projects/Cho2010/nystrom_base.py, com a
//                          geometria do talude) no lugar da analitica do retangulo; o nucleo vem do arquivo e tem de ter
//                          as distancias 20 s e 2 s (exponencial)
//   xi (<arquivo>)         os xi dados, sem sorteio (um ponto do FORM, Projects/Cho2010/form_kl.py prog=chuva): 3K numeros,
//                          os K de c', os K de phi' e os K de k (K <= M; os outros termos ficam em 0); a realizacao
//                          primeira usa esses xi
//   desloc_xi (<arquivo>)  amostragem por importancia com o deslocamento dado termo a termo: 3K numeros (c', phi', k, como
//                          o u* do form_kl.py); os outros termos ficam sem deslocamento. Use com outra semente=
//   gradiente (0)          K > 0: o gradiente dFS/dxi (SRM) ou dGamma/dxi (GIM) nos K primeiros termos de cada campo pela
//                          projecao modal da sensibilidade de limite superior (Mohr-Coulomb associado, tensoes efetivas
//                          do colapso, pressoes congeladas), no instante do menor fator: <saida>_grad_<real>.txt (o menor
//                          fator e os 3K numeros, c', phi' e k; os de k sao 0: o k so muda as pressoes, congeladas no
//                          fator)
// O transiente reaproveitado (so com k fixo: covk = 0 e o k de valores = 1):
//   grava_pressoes (<arquivo>)  grava p_w nos nos da malha do transiente em cada instante do fator, com os dados do
//                          transiente (geometria, contornos, solo, chuva, malha) para conferir na leitura
//   pressoes (<arquivo>)   le p_w desse arquivo em vez de calcular o transiente (os dados do transiente e a malha tem de ser
//                          os mesmos; os tempos pedidos tem de estar nele)
//   so_transiente (0)      1: so o transiente (com grava_pressoes=), sem fator
// Amostragem por importancia (AI; ImportanceSampling.h): o deslocamento da media dos xi sai de um Monte Carlo direto ja
// feito (o piloto) pela entropia cruzada, e cada realizacao grava o peso exato (log_peso): P_f = soma(peso * falha) / N
// nao tem vies (posproc/amostragem_importancia.py). A falha e FS < 1 em algum instante: no piloto vale o menor FS.
//   banda (<arquivo>), banda_t (todos)   o mecanismo.csv de um deterministico e o instante dele: os pontos com
//                          sqrt(J2) >= 0,1 do maior dao a direcao d dos xi (o gradiente da media do campo gaussiano na
//                          banda); com banda, o arquivo grava as projecoes xi . d (proj_c, proj_phi, proj_k)
//   piloto (<pasta>)       a pasta do piloto (as realizacoes 0..piloto_n-1 dos <prefixo>_<a>_<b>.csv); cada bloco calcula
//                          o mesmo deslocamento (<saida>_desloc_<faixa>.txt). Sai com o codigo 3 se o piloto tiver menos
//                          de 90% das realizacoes
//   piloto_n (1000), piloto_semente (2025), elite (0.1)   o conjunto elite: as falhas do piloto ou, se forem menos que a
//                          fracao elite, as de menor FS
//   ai (banda), modos (30) banda: s_c d por campo, com s_c a media de xi_c . d no conjunto elite (precisa de banda);
//                          modos: a media dos 'modos' primeiros modos de cada campo
// Saidas (<faixa> = primeira_ultima):
//   <saida>_<faixa>.csv       uma linha por realizacao e instante, com as colunas do mc_resumo.csv do rebaixamento
//                             (FS_min, situacao, banda_*, massa_*, c_med, phi_med, k_mult_med, parametros do caso) e as da
//                             chuva (nivel d'agua na crista, no meio da face e no pe, p_w na banda, nos na face, ...);
//                             niveis: o FS de cada nivel do refinamento
//   <saida>_parametros.txt    chave=valor do caso
//   <saida>_mapa_<faixa>.csv  por instante e celula de 1 m: quantas realizacoes tem a banda do mecanismo ali (n) e quantas
//                             dessas tem FS < 1 (n_falha); n_real: realizacoes somadas (somar os arquivos das partes)
//   <saida>_real<n>/          campos de uma realizacao: o formato do LearningMohrCoulomb2/animacao.py, campos.csv (c', phi',
//                             k nos nos) e mecanismo.csv (pontos com sqrt(J2) >= 1% do maior, como no rebaixamento)
//   <saida>_grad_<n>.txt      com gradiente=K
#include "AnimationFiles.h"
#include "HPMaterial.h"
#include "ImportanceSampling.h"
#include "LayerGMesh.h"
#include "MualemTabelado.h"
#include "PointLocator.h"
#include "RainTransient.h"
#include "RandomField.h"
#include "UnsatAnalysis.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "TPZSpStructMatrix.h"
#include "pzgeotriangle.h"
#include "pzstepsolver.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;
using TUnsatHP = HPMaterial<TUnsat>;
using Clock = std::chrono::steady_clock;
// propriedades num ponto (x, y): {c' (kPa), phi' (graus), k / k_s}
using Campo = std::function<std::vector<double>(double, double)>;
static double Seg(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// ids da TriGMesh do main
enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };

// o talude do problema: cantos no sentido anti-horario (base esquerda, base direita, terreno do pe, pe, crista, crista a
// esquerda), como no main.cpp. Padrao: o da Fig. 2 de Rahardjo et al. (2007), H = 10 m a 45 graus; a opcao geometria=
// troca (o LeOpcoes, antes de tudo)
std::vector<std::vector<REAL>> gTalude = {{0, 0}, {70, 0}, {70, 20}, {40, 20}, {30, 30}, {0, 30}};

// ---- dados do problema
struct Solo { // Tabela 1 do Rahardjo: o que muda de um solo para outro
    std::string nome;
    REAL a, ks; // parametro a de Fredlund-Xing (kPa) e k_s (m/s)
};
const std::vector<Solo> kSolos = {{"f10", 10., 1.e-4}, {"f50", 50., 1.e-5}, {"f100", 100., 1.e-6}};
struct Resistencia { // igual nos tres solos do Rahardjo (opcoes c=, phi=, gamma=, E=, nu=)
    REAL c = 10., phiDeg = 26.; // medias de c' (kPa) e phi' (graus)
    REAL gamma = 20., gammaw = 9.81, E = 2.e4, nu = 0.3;
};
struct Campos { // Tabela 2 do Ceron: CoV e distancias de autocorrelacao (m); k: multiplicador de k_s com media 1
    REAL covC = 0.3, covPhi = 0.1, covK = 0.6, Lx = 20., Ly = 2., rcphi = 0.;
    int M = 2000;
    uint64_t semente = 2025; // a realizacao n tem os mesmos numeros aleatorios em todos os casos
};
struct VanGenuchten { // opcao vg=: theta_s, theta_r, alfa (1/kPa), n
    bool usa = false;
    REAL ts = 0., tr = 0., a = 0., n = 0.;
};

// ---- opcoes (o cabecalho explica cada uma)
struct Opcoes {
    int n = 10, primeira = 0;
    bool det = false, gim = false;
    Solo solo = kSolos[1];
    REAL I = -1., chuva = 24., hw = 5., alfa = 1., s = 1.;
    std::vector<REAL> tempos = {0., 24.}, adapta;
    Campos campos;
    std::string malhaT = "tri2", malhaF = "tri1", saida = "mc_chuva";
    int uT = 2, nref = 4;
    std::set<int> figuras;
    std::vector<REAL> valores; // deterministica com c' (kPa), phi' (graus) e k / k_s dados (vazio: as medias)
    std::string banda;         // amostragem por importancia: o mecanismo.csv da banda e o instante dele
    REAL bandaT = -1., srmF0 = 0.8;
    ai::Opcoes ai;             // o deslocamento pelo piloto
    // o problema
    REAL betaNA = 7., succaoMax = 75.;
    bool esqCarga = true, dirCarga = true, chuvaVertical = false, chiSe = false;
    REAL cargaEsq = std::numeric_limits<REAL>::quiet_NaN(), cargaDir = std::numeric_limits<REAL>::quiet_NaN();
    VanGenuchten vg;
    REAL vgTab = 0.;                       // vg_tab: o eps do teto S_e = 1 - eps (0: o EMualem)
    EFechoVG vgFecho = EFechoVG::Reta;
    Resistencia res;
    // os campos
    std::string nystrom, xi, deslocXi;
    int gradiente = 0;
    // o transiente reaproveitado
    std::string pressoes, gravaPressoes;
    bool soTransiente = false;
    // casos B e C: Fredlund-Xing dado, k de Fredlund, estado inicial permanente, chi por ponto, chuva aleatoria e CRPC
    std::vector<REAL> fx;                  // a, n, m, theta_s (vazio: o solo)
    REAL kFredlund = -1., psiR = 3000.;    // psi_aev (kPa; < 0: k_s Theta^4) e psi_r
    REAL q0 = 0., q0Cov = 0.;              // infiltracao do estado inicial (mm/h) e o CoV dela
    bool chiPonto = false;
    REAL phibDeg = 26., phibRazao = -1.;   // phi^b (graus) ou phi^b / phi' (> 0)
    REAL dt = 0.1, dtMin = 0.;             // passo do transiente e o menor da subdivisao (h)
    std::vector<REAL> newton;              // tol_p, maxit do Newton do transiente (vazio: os do TPZUnsatAnalysis)
    std::string newtonTrace;               // arquivo do diagnostico do Newton (vazio: sem)
    REAL faceTol = 0.;                     // folga relativa da volta da face ao fluxo
    bool chuvaAleatoria = false;
    std::array<REAL, 7> gpd = {-0.04, 3.44, 5.13, 0.43, 5.15, 0.30, -2.38};
    std::vector<REAL> crpc;                // intensidades da CRPC (vazio: Monte Carlo)
    REAL crpcDmax = 200., crpcPasso = 1.;
};

std::vector<REAL> Lista(const std::string &s)
{
    std::vector<REAL> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(std::stod(item));
    return v;
}

[[noreturn]] void Erro(const std::string &m)
{
    std::fprintf(stderr, "%s\n", m.c_str());
    std::exit(1);
}

Opcoes LeOpcoes(int argc, char *argv[])
{
    const std::set<std::string> chaves = {
        "n",          "primeira",   "det",        "solo",       "I",          "chuva",     "hw",        "tempos",
        "adapta",     "metodo",     "covc",       "covphi",     "covk",       "s",         "alfa",      "transiente",
        "uT",         "fator",      "nref",       "figuras",    "saida",      "valores",   "semente",   "banda",
        "banda_t",    "piloto",     "piloto_n",   "piloto_semente", "ai",     "elite",     "modos",     "srm_f0",
        "geometria",  "na",         "na_beta",    "succao_max", "esquerda",   "direita",   "carga_esq", "carga_dir",
        "chuva_face", "ks",         "vg",         "chi",        "c",          "phi",       "gamma",     "E",
        "nu",         "rcphi",      "M",          "nystrom",    "xi",         "desloc_xi", "gradiente", "grava_pressoes",
        "pressoes",   "so_transiente", "fx",       "kfredlund", "psi_r",     "q0",         "q0_cov",    "phib",
        "phib_razao", "dt",         "dt_min",    "newton",         "chuva_aleatoria", "gpd",  "crpc",      "crpc_dmax", "crpc_passo",
        "vg_tab",     "vg_tab_psi", "newton_trace", "face_tol"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        if (eq == std::string::npos || !chaves.count(a.substr(0, eq)))
            Erro("opcao desconhecida: " + a + " (as opcoes estao no cabecalho do montecarlo.cpp)");
        kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    auto tem = [&kv](const char *k) { return kv.count(k) > 0; };
    Opcoes op;
    if (tem("geometria")) {
        const std::vector<REAL> v = Lista(kv["geometria"]);
        if (v.size() != 12) Erro("geometria: 12 numeros (os 6 pontos x,y do talude)");
        for (int i = 0; i < 6; i++) gTalude[i] = {v[2 * i], v[2 * i + 1]};
        const auto &p = gTalude;
        if (p[1][1] != p[0][1] || p[2][0] != p[1][0] || p[3][1] != p[2][1] || p[5][0] != p[0][0] || p[5][1] != p[4][1] ||
            !(p[0][0] < p[4][0] && p[4][0] < p[3][0] && p[3][0] < p[1][0] && p[0][1] < p[3][1] && p[3][1] < p[4][1]))
            Erro("geometria: os pontos tem de ser base esquerda, base direita, terreno do pe, pe, crista e crista a "
                 "esquerda (base e terrenos horizontais, laterais verticais, crista a esquerda e acima do pe)");
    }
    if (tem("n")) op.n = std::stoi(kv["n"]);
    if (tem("primeira")) op.primeira = std::stoi(kv["primeira"]);
    if (tem("det")) op.det = std::stoi(kv["det"]) != 0;
    if (tem("solo")) {
        auto it = std::find_if(kSolos.begin(), kSolos.end(), [&](const Solo &s) { return s.nome == kv["solo"]; });
        if (it == kSolos.end()) Erro("solo desconhecido: " + kv["solo"] + " (use f10, f50 ou f100)");
        op.solo = *it;
    }
    if (tem("ks")) op.solo.ks = std::stod(kv["ks"]);
    if (!(op.solo.ks > 0.)) Erro("ks: k_s > 0 (m/s)");
    op.I = tem("I") ? std::stod(kv["I"]) : op.solo.ks * 3.6e6; // padrao: 1 k_s (m/s -> mm/h)
    if (tem("chuva")) op.chuva = std::stod(kv["chuva"]);
    if (tem("hw") && tem("na")) Erro("hw e na: um ou outro");
    if (tem("hw")) op.hw = std::stod(kv["hw"]);
    if (tem("na")) op.hw = gTalude[3][1] - std::stod(kv["na"]);
    if (tem("na_beta")) op.betaNA = std::stod(kv["na_beta"]);
    if (tem("succao_max")) op.succaoMax = std::stod(kv["succao_max"]);
    auto lateral = [&](const char *k, bool &carga) {
        if (!tem(k)) return;
        if (kv[k] != "carga" && kv[k] != "fechada") Erro(std::string(k) + ": carga ou fechada");
        carga = kv[k] == "carga";
    };
    lateral("esquerda", op.esqCarga);
    lateral("direita", op.dirCarga);
    if (tem("carga_esq")) op.cargaEsq = std::stod(kv["carga_esq"]);
    if (tem("carga_dir")) op.cargaDir = std::stod(kv["carga_dir"]);
    if ((tem("carga_esq") && !op.esqCarga) || (tem("carga_dir") && !op.dirCarga))
        Erro("carga_esq/carga_dir so numa lateral com carga");
    if (tem("chuva_face")) {
        if (kv["chuva_face"] != "inclinada" && kv["chuva_face"] != "vertical") Erro("chuva_face: inclinada ou vertical");
        op.chuvaVertical = kv["chuva_face"] == "vertical";
    }
    if (tem("vg")) {
        const std::vector<REAL> v = Lista(kv["vg"]);
        if (v.size() != 4 || !(v[0] > v[1] && v[1] >= 0. && v[2] > 0. && v[3] > 1.))
            Erro("vg: theta_s,theta_r,alfa,n com theta_s > theta_r >= 0, alfa > 0 (1/kPa) e n > 1");
        op.vg = {true, v[0], v[1], v[2], v[3]};
        if (!tem("solo")) op.solo.nome = "vg"; // a coluna solo do CSV
    }
    if (tem("vg_tab") && tem("vg_tab_psi")) Erro("vg_tab e vg_tab_psi: um ou outro");
    for (const std::string k : {"vg_tab", "vg_tab_psi"}) { // k_r de Mualem tabelado: eps (ou psi_c) e o fecho
        if (!tem(k.c_str())) continue;
        const size_t vir = kv[k].find(',');
        const REAL x = std::stod(kv[k].substr(0, vir));
        if (k == "vg_tab" && x == 0.) break; // vg_tab=0: desligado
        if (!op.vg.usa) Erro(k + " precisa de vg=");
        const std::string fecho = vir == std::string::npos ? "reta" : kv[k].substr(vir + 1);
        if (fecho != "reta" && fecho != "c1" && fecho != "ippisch") Erro(k + ": o fecho e reta, c1 ou ippisch");
        op.vgFecho = fecho == "reta" ? EFechoVG::Reta : fecho == "c1" ? EFechoVG::C1 : EFechoVG::Ippisch;
        op.vgTab = k == "vg_tab" ? x : x > 0. ? EpsDePsi(op.vg.a, op.vg.n, x) : 0.;
        if (!(op.vgTab >= 1.e-5 && op.vgTab <= 0.1)) {
            char b[128];
            std::snprintf(b, sizeof b, "%s: eps = %g fora de [1e-5, 0.1]%s", k.c_str(), op.vgTab,
                          k == "vg_tab" ? "" : " (eps = 1 - S_e(psi_c))");
            Erro(b);
        }
    }
    if (tem("chi")) {
        if (kv["chi"] != "const" && kv["chi"] != "se" && kv["chi"] != "ponto") Erro("chi: const, se ou ponto");
        op.chiSe = kv["chi"] == "se";
        op.chiPonto = kv["chi"] == "ponto";
    }
    if (tem("phib")) op.phibDeg = std::stod(kv["phib"]);
    if (tem("phib_razao")) op.phibRazao = std::stod(kv["phib_razao"]);
    if (tem("phib_razao") && !(op.chiPonto && op.phibRazao > 0.)) Erro("phib_razao: > 0 e so com chi=ponto");
    if (tem("phib") && op.chiSe) Erro("phib nao vale com chi=se");
    if (tem("fx")) {
        op.fx = Lista(kv["fx"]);
        if (op.fx.size() != 4 || !(op.fx[0] > 0. && op.fx[1] > 0. && op.fx[2] > 0. && op.fx[3] > 0. && op.fx[3] < 1.))
            Erro("fx: a,n,m,theta_s (a em kPa, todos > 0, theta_s < 1)");
        if (tem("vg")) Erro("fx e vg: um ou outro");
        op.solo.a = op.fx[0];
        if (!tem("solo")) op.solo.nome = "fx";
    }
    if (tem("kfredlund")) op.kFredlund = std::stod(kv["kfredlund"]);
    if (tem("kfredlund") && (op.vg.usa || !(op.kFredlund > 0.))) Erro("kfredlund: psi_aev > 0 (kPa), so com Fredlund-Xing");
    if (tem("psi_r")) op.psiR = std::stod(kv["psi_r"]);
    if (tem("q0")) op.q0 = std::stod(kv["q0"]);
    if (tem("q0_cov")) op.q0Cov = std::stod(kv["q0_cov"]);
    if (op.q0 < 0. || op.q0Cov < 0. || (op.q0Cov > 0. && op.q0 == 0.)) Erro("q0 >= 0 (mm/h) e q0_cov >= 0 (so com q0 > 0)");
    if (tem("dt")) op.dt = std::stod(kv["dt"]);
    if (!(op.dt > 0.)) Erro("dt > 0 (h)");
    if (tem("dt_min")) op.dtMin = std::stod(kv["dt_min"]);
    if (tem("newton")) {
        op.newton = Lista(kv["newton"]);
        if (op.newton.size() < 2 || op.newton.size() > 3 || !(op.newton[0] > 0. && op.newton[1] >= 1.))
            Erro("newton: tol_p,maxit[,busca] (kPa, >= 1, metades da busca linear)");
    }
    if (tem("newton_trace")) op.newtonTrace = kv["newton_trace"];
    if (tem("face_tol")) op.faceTol = std::stod(kv["face_tol"]);
    if (!(op.faceTol >= 0.)) Erro("face_tol >= 0");
    if (tem("chuva_aleatoria")) op.chuvaAleatoria = std::stoi(kv["chuva_aleatoria"]) != 0;
    if (tem("gpd")) {
        const std::vector<REAL> v = Lista(kv["gpd"]);
        if (v.size() != 7 || !(v[2] > 0. && v[5] > 0.) || v[6] == 0.)
            Erro("gpd: xi_I,mu_I,sigma_I,xi_D,mu_D,sigma_D,omega (escalas > 0, omega != 0)");
        std::copy(v.begin(), v.end(), op.gpd.begin());
    }
    if (tem("crpc")) {
        op.crpc = Lista(kv["crpc"]);
        if (op.crpc.empty() || !std::is_sorted(op.crpc.rbegin(), op.crpc.rend()) || !(op.crpc.back() > 0.))
            Erro("crpc: intensidades > 0 em ordem decrescente (mm/h)");
    }
    if (tem("crpc_dmax")) op.crpcDmax = std::stod(kv["crpc_dmax"]);
    if (tem("crpc_passo")) op.crpcPasso = std::stod(kv["crpc_passo"]);
    if (!op.crpc.empty() && (op.chuvaAleatoria || tem("pressoes") || tem("grava_pressoes") || tem("gradiente")))
        Erro("crpc: sem chuva_aleatoria, pressoes, grava_pressoes e gradiente");
    if (op.chuvaAleatoria && (tem("pressoes") || tem("grava_pressoes")))
        Erro("chuva_aleatoria: o transiente muda a cada realizacao (sem pressoes e grava_pressoes)");
    if (tem("c")) op.res.c = std::stod(kv["c"]);
    if (tem("phi")) op.res.phiDeg = std::stod(kv["phi"]);
    if (tem("gamma")) op.res.gamma = std::stod(kv["gamma"]);
    if (tem("E")) op.res.E = std::stod(kv["E"]);
    if (tem("nu")) op.res.nu = std::stod(kv["nu"]);
    if (tem("tempos")) op.tempos = Lista(kv["tempos"]);
    std::sort(op.tempos.begin(), op.tempos.end());
    op.adapta = tem("adapta") ? Lista(kv["adapta"]) : op.tempos;
    for (REAL t : op.adapta)
        if (std::find(op.tempos.begin(), op.tempos.end(), t) == op.tempos.end())
            Erro("adapta: " + std::to_string(t) + " nao esta em tempos");
    if (op.adapta.empty()) Erro("adapta precisa de pelo menos um instante");
    if (tem("metodo")) {
        if (kv["metodo"] != "srm" && kv["metodo"] != "gim") Erro("metodo: srm ou gim");
        op.gim = kv["metodo"] == "gim";
    }
    if (tem("covc")) op.campos.covC = std::stod(kv["covc"]);
    if (tem("covphi")) op.campos.covPhi = std::stod(kv["covphi"]);
    if (tem("covk")) op.campos.covK = std::stod(kv["covk"]);
    if (tem("rcphi")) op.campos.rcphi = std::stod(kv["rcphi"]);
    if (!(std::abs(op.campos.rcphi) < 1.)) Erro("rcphi precisa estar entre -1 e 1");
    if (tem("M")) op.campos.M = std::stoi(kv["M"]);
    if (op.campos.M < 1) Erro("M >= 1");
    if (tem("s")) op.s = std::stod(kv["s"]);
    op.campos.Lx = 20. * op.s, op.campos.Ly = 2. * op.s;
    if (tem("alfa")) op.alfa = std::stod(kv["alfa"]);
    if (tem("transiente")) op.malhaT = kv["transiente"];
    if (tem("uT")) op.uT = std::stoi(kv["uT"]);
    if (tem("fator")) op.malhaF = kv["fator"];
    if (tem("nref")) op.nref = std::stoi(kv["nref"]);
    if (tem("figuras"))
        for (REAL v : Lista(kv["figuras"])) op.figuras.insert(int(v));
    if (tem("saida")) op.saida = kv["saida"];
    if (tem("semente")) op.campos.semente = std::stoull(kv["semente"]);
    if (tem("banda")) op.banda = kv["banda"];
    if (tem("banda_t")) op.bandaT = std::stod(kv["banda_t"]);
    if (tem("srm_f0")) op.srmF0 = std::stod(kv["srm_f0"]);
    op.ai.semente = tem("piloto_semente") ? std::stoull(kv["piloto_semente"]) : Campos{}.semente;
    if (tem("piloto")) op.ai.piloto = kv["piloto"];
    if (tem("piloto_n")) op.ai.n = std::stoi(kv["piloto_n"]);
    if (tem("ai")) op.ai.tipo = kv["ai"];
    if (tem("elite")) op.ai.elite = std::stod(kv["elite"]);
    if (tem("modos")) op.ai.modos = std::stoi(kv["modos"]);
    if (op.ai.tipo != "banda" && op.ai.tipo != "modos") Erro("ai=" + op.ai.tipo + ": banda ou modos");
    if (!op.ai.piloto.empty() && op.ai.tipo == "banda" && op.banda.empty())
        Erro("piloto com ai=banda precisa de banda=<mecanismo.csv>");
    if (tem("valores")) {
        op.valores = Lista(kv["valores"]);
        if (op.valores.size() != 3) Erro("valores: c,phi,k (c' em kPa, phi' em graus, k / k_s)");
    }
    if (op.malhaF.rfind("camq", 0) == 0)
        Erro("a malha do fator precisa de triangulos (a bisseccao do refinamento): use triK ou camH");
    if (tem("nystrom")) op.nystrom = kv["nystrom"];
    if (tem("xi")) op.xi = kv["xi"];
    if (tem("desloc_xi")) op.deslocXi = kv["desloc_xi"];
    if (!op.deslocXi.empty() && !op.ai.piloto.empty()) Erro("desloc_xi e piloto: um ou outro");
    if (tem("gradiente")) op.gradiente = std::stoi(kv["gradiente"]);
    if (tem("grava_pressoes")) op.gravaPressoes = kv["grava_pressoes"];
    if (tem("pressoes")) op.pressoes = kv["pressoes"];
    if (tem("so_transiente")) op.soTransiente = std::stoi(kv["so_transiente"]) != 0;
    if (op.soTransiente && op.gravaPressoes.empty()) Erro("so_transiente=1 precisa de grava_pressoes=<arquivo>");
    if (!op.pressoes.empty() && !op.gravaPressoes.empty()) Erro("pressoes e grava_pressoes: um ou outro");
    const bool kFixo = op.campos.covK == 0. && (op.valores.empty() || op.valores[2] == 1.) && op.q0Cov == 0.;
    if ((!op.pressoes.empty() || !op.gravaPressoes.empty()) && !kFixo)
        Erro("pressoes/grava_pressoes so com k fixo (covk=0 e, com valores, k = 1): o transiente tem de ser o mesmo "
             "em todas as realizacoes");
    return op;
}

// ---- malhas
// TriGMesh do main.cpp (a mesma funcao, sem gravar o VTK): passo H, cada trecho em round(comprimento / H) partes, cada
// retangulo em dois triangulos e a cunha da face num triangulo so; ref refinamentos uniformes
TPZGeoMesh *TriGMesh(const std::vector<std::vector<REAL>> &points, int ref)
{
    if (points.size() != 6) DebugStop();
    const REAL x0 = points[0][0], y0 = points[0][1], L = points[1][0], yPe = points[2][1];
    const REAL xPe = points[3][0], xCr = points[4][0], yCr = points[4][1];
    if (points[1][1] != y0 || points[2][0] != L || points[3][1] != yPe || points[5][0] != x0 || points[5][1] != yCr ||
        !(x0 < xCr && xCr < xPe && xPe < L && y0 < yPe && yPe < yCr))
        DebugStop();
    const REAL H = yCr - yPe;
    auto divide = [H](REAL a, REAL b) {
        const int n = std::max(1, int(std::lround((b - a) / H)));
        std::vector<REAL> v(n + 1);
        for (int k = 0; k <= n; k++) v[k] = a + (b - a) * k / n;
        return v;
    };
    std::vector<REAL> xs = divide(x0, xCr);
    const int iCr = int(xs.size()) - 1, iPe = iCr + 1;
    for (REAL x : divide(xPe, L)) xs.push_back(x);
    const std::vector<REAL> ys = divide(y0, yPe);
    const int nx = int(xs.size()) - 1, ny = int(ys.size()) - 1;
    auto no = [nx](int i, int j) { return int64_t(j) * (nx + 1) + i; };
    auto topo = [nx, ny](int i) { return int64_t(ny + 1) * (nx + 1) + i; };
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(topo(iCr) + 1);
    auto poeNo = [gmesh](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        gmesh->NodeVec()[id] = TPZGeoNode(id, c, *gmesh);
    };
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) poeNo(no(i, j), xs[i], ys[j]);
    for (int i = 0; i <= iCr; i++) poeNo(topo(i), xs[i], yCr);
    auto triangulo = [gmesh](int64_t a, int64_t b, int64_t c) {
        TPZManVector<int64_t, 3> nodes = {a, b, c};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, matSolo, *gmesh);
    };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            triangulo(no(i, j), no(i + 1, j), no(i, j + 1));
            triangulo(no(i + 1, j), no(i + 1, j + 1), no(i, j + 1));
        }
    for (int i = 0; i < iCr; i++) {
        triangulo(no(i, ny), no(i + 1, ny), topo(i));
        triangulo(no(i + 1, ny), topo(i + 1), topo(i));
    }
    triangulo(no(iCr, ny), no(iPe, ny), topo(iCr));
    auto linha = [gmesh](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> nodes = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, id, *gmesh);
    };
    for (int i = 0; i < nx; i++) linha(no(i, 0), no(i + 1, 0), base);
    for (int j = 0; j < ny; j++) linha(no(nx, j), no(nx, j + 1), right);
    for (int i = nx; i > iPe; i--) linha(no(i, ny), no(i - 1, ny), topright);
    for (int i = iCr; i > 0; i--) linha(topo(i), topo(i - 1), topleft);
    linha(topo(0), no(0, ny), left);
    for (int j = ny; j > 0; j--) linha(no(0, j), no(0, j - 1), left);
    linha(no(iPe, ny), topo(iCr), face);
    gmesh->BuildConnectivity();
    for (int d = 0; d < ref; d++) {
        const int64_t nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> sub;
        for (int64_t iel = 0; iel < nel; iel++)
            if (!gmesh->Element(iel)->HasSubElement()) gmesh->Element(iel)->Divide(sub);
    }
    return gmesh;
}

// Copia so com as folhas: o TPZUnsatAnalysis acha os nos da chuva pelos elementos de contorno da malha geometrica, e
// com os pais junto a chuva entraria duas vezes nos mesmos nos
TPZGeoMesh *MalhaPlana(TPZGeoMesh *gmesh)
{
    UnsatModel<TUnsatHP> tmp;
    tmp.SetGeoMesh(gmesh); // tmp fica dono da malha antiga e a apaga
    std::vector<int64_t> pai;
    return tmp.RefinedGeoMesh({}, pai);
}

// Malha em camadas do mesmo talude (LayerGMesh): passo h nas colunas da face e nas camadas ate a profundidade D (na
// coluna mais alta; D = 0: so a primeira camada), crescendo 1,3 vezes para longe do talude e para baixo, ate 10 m. Na face
// as colunas tem 0,75 h: a camada ali tem de 0,75 h a h na vertical, e coluna mais estreita que a camada nao deixa
// triangulo obtuso
TPZGeoMesh *Camadas(REAL h, REAL D, bool quad)
{
    const auto &p = gTalude;
    const REAL r = 1.3, hmax = 10., xCrista = p[4][0], xPe = p[3][0], xFim = p[1][0], H = p[4][1] - p[0][1];
    LayerGMeshData d;
    d.superficie = {{p[0][0], p[4][1]}, {xCrista, p[4][1]}, {xPe, p[3][1]}, {xFim, p[3][1]}};
    d.idSuperficie = {topleft, face, topright};
    d.idMaterial = matSolo, d.idBase = base, d.idEsquerda = left, d.idDireita = right;
    const std::vector<REAL> esq = GeometricSteps(xCrista - p[0][0], h, r, hmax), dir = GeometricSteps(xFim - xPe, h, r, hmax);
    const int nf = std::max(1, int(std::ceil((xPe - xCrista) / (0.75 * h))));
    d.xs = {p[0][0]};
    for (auto it = esq.rbegin(); it != esq.rend(); ++it) d.xs.push_back(d.xs.back() + *it); // fino perto da crista
    d.xs.back() = xCrista;
    for (int k = 1; k <= nf; k++) d.xs.push_back(xCrista + (xPe - xCrista) * k / nf);
    for (REAL q : dir) d.xs.push_back(d.xs.back() + q); // fino perto do pe
    d.xs.back() = xFim;
    d.eta = {1.}; // da superficie para a base: camadas de h ate D, depois crescendo
    const int nD = int(std::lround(D / h));
    for (int k = 0; k < nD; k++) d.eta.push_back(d.eta.back() - h / H);
    for (REAL q : nD ? GeometricSteps(H - nD * h, h * r, r, hmax) : GeometricSteps(H, h, r, hmax))
        d.eta.push_back(d.eta.back() - q / H);
    d.eta.back() = 0.;
    std::reverse(d.eta.begin(), d.eta.end());
    d.quadrilateros = quad;
    return LayerGMesh(d);
}

// Malha pelo nome: triK (TriGMesh com K divisoes, sem hierarquia), camH ou camqH (camadas, triangulos ou
// quadrilateros); camHpD e camqHpD: com as camadas de H ate a profundidade D
TPZGeoMesh *Malha(const std::string &nome)
{
    auto camadas = [](const std::string &s, bool quad) {
        const size_t p = s.find('p');
        return Camadas(std::stod(s.substr(0, p)), p == std::string::npos ? 0. : std::stod(s.substr(p + 1)), quad);
    };
    if (nome.rfind("camq", 0) == 0) return camadas(nome.substr(4), true);
    if (nome.rfind("cam", 0) == 0) return camadas(nome.substr(3), false);
    if (nome.rfind("tri", 0) == 0) return MalhaPlana(TriGMesh(gTalude, std::stoi(nome.substr(3))));
    Erro("malha desconhecida: " + nome + " (use triK, camH[pD] ou camqH[pD])");
}

// cota da superficie do talude em x
REAL Superficie(REAL x)
{
    const auto &p = gTalude;
    if (x <= p[4][0]) return p[4][1];
    if (x >= p[3][0]) return p[3][1];
    return p[4][1] + (p[3][1] - p[4][1]) * (x - p[4][0]) / (p[3][0] - p[4][0]);
}

// ---- k_r de Fredlund et al. (1994) (opcao kfredlund=; Liu et al. 2024, Eq. 9), tabelado contra S_w:
//   k_r(psi) = int_{ln psi}^{ln 1e6} [theta(e^y) - theta(psi)] theta'(e^y) / e^y dy
//            / int_{ln psi_aev}^{ln 1e6} [theta(e^y) - theta_s] theta'(e^y) / e^y dy,   limitado a 1,
// com theta = theta_s C(psi) / ln(e + (psi/a)^n)^m e C(psi) = 1 - ln(1 + psi/psi_r) / ln(1 + 1e6/psi_r) (Eqs. 7 e 8). A
// tabela vai em S_w da retencao do material (Fredlund-Xing com C = 1) no mesmo psi: o k e o de Eq. 9 em cada succao. As
// integrais saem da regra do trapezio em y com 40 000 intervalos, de ln(1e-4) a ln(1e6), acumuladas de cima para baixo
struct TabelaK {
    std::vector<REAL> S, lnK;
};
TabelaK gTabelaK; // preenchida no main com kfredlund=
TabelaVG gTabelaVG; // k_r de Mualem tabelado ate S_e = 1 - eps (MualemTabelado.h), preenchida no main com vg_tab=

TabelaK TabelaFredlund(REAL a, REAL n, REAL m, REAL ts, REAL psiR, REAL psiAev)
{
    auto C = [psiR](REAL p) { return 1. - std::log(1. + p / psiR) / std::log(1. + 1.e6 / psiR); };
    auto theta = [&](REAL p) { return ts * C(p) / std::pow(std::log(M_E + std::pow(p / a, n)), m); };
    auto dtheta = [&](REAL p) { const REAL h = 1.e-5 * p; return (theta(p + h) - theta(p - h)) / (2. * h); };
    const int N = 40000;
    const REAL y0 = std::log(1.e-4), y1 = std::log(1.e6), dy = (y1 - y0) / N;
    // A(y) = int_y^y1 theta theta' / e^y, B(y) = int_y^y1 theta' / e^y (acumuladas do fim)
    std::vector<REAL> A(N + 1, 0.), B(N + 1, 0.), fa(N + 1), fb(N + 1);
    for (int i = 0; i <= N; i++) {
        const REAL p = std::exp(y0 + i * dy), d = dtheta(p) / p;
        fa[i] = theta(p) * d, fb[i] = d;
    }
    for (int i = N - 1; i >= 0; i--) A[i] = A[i + 1] + 0.5 * dy * (fa[i] + fa[i + 1]), B[i] = B[i + 1] + 0.5 * dy * (fb[i] + fb[i + 1]);
    auto integral = [&](REAL p, REAL t) { // int_{ln p}^{y1} (theta - t) theta' / e^y, interpolado na grade
        const REAL x = std::clamp((std::log(p) - y0) / dy, REAL(0.), REAL(N));
        const int i = std::min(int(x), N - 1);
        const REAL f = x - i;
        return (1. - f) * (A[i] - t * B[i]) + f * (A[i + 1] - t * B[i + 1]);
    };
    const REAL den = integral(psiAev, ts);
    TabelaK T;
    for (int i = 0; i <= 400; i++) { // psi de 1e6 a 1e-3 kPa: S_w crescente
        const REAL p = std::pow(10., 6. - 9. * i / 400.);
        const REAL S = std::pow(std::log(M_E + std::pow(p / a, n)), -m); // a retencao do material (C = 1)
        const REAL kr = std::min(REAL(1.), std::max(integral(p, theta(p)) / den, REAL(1.e-30)));
        if (!T.S.empty() && !(S > T.S.back())) continue;
        T.S.push_back(S), T.lnK.push_back(std::log(kr));
    }
    if (T.S.back() < 1.) T.S.push_back(1.), T.lnK.push_back(T.lnK.back());
    return T;
}

// ---- material do solo (o MaterialRahardjo do rahardjo.cpp), com a anisotropia k_h / k_v; com vg=, van Genuchten-Mualem;
// com chi=se, chi = S_e
TUnsatHP *MaterialChuva(const RainData &r, const Opcoes &op)
{
    const Resistencia &s = op.res;
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    if (op.vg.usa) { // S_w = S_wr + (1 - S_wr) S_e com theta = theta_s S_w; k_r de Mualem
        h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        h.fAVG = op.vg.a, h.fBVG = op.vg.n, h.fSwr = op.vg.tr / op.vg.ts;
        h.fKrw = TPZUnsaturatedHydraulics::EMualem;
        if (op.vgTab > 0.) UsaTabela(h, gTabelaVG); // vg_tab=: Mualem tabelado ate S_e = 1 - eps
    } else {
        h.fRetention = TPZUnsaturatedHydraulics::EFredlundXing; // S_w(p_c) com C(psi) = 1
        h.fA = r.aFX, h.fN = r.nFX, h.fM = r.mFX;
        h.fKrw = TPZUnsaturatedHydraulics::EPowerLaw; // k_w = k_s Theta^p
        h.fPower = r.power;
        if (op.kFredlund > 0.) { // k_r de Fredlund et al. (1994), tabelado
            h.fKrw = TPZUnsaturatedHydraulics::ETabulated;
            h.fTabS = gTabelaK.S, h.fTabLnK = gTabelaK.lnK;
        }
    }
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // o ar fica passivo (p_a = 0)
    mat->SetHydraulics(h);
    if (op.chiSe) mat->SetBishop(TUnsat::EChiEffective); // chi = S_e (Vanapalli et al. 1996)
    else if (op.chiPonto) mat->SetBishop(TUnsat::EChiConstant, 1.); // chi(x) vai nas pressoes congeladas (CalculaFator)
    else mat->SetBishop(TUnsat::EChiConstant, std::tan(r.phib) / std::tan(phi)); // (u_a - u_w) tan(phi^b)
    mat->SetPorosity(r.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(r.Kw);
    mat->SetConductivity(r.ks, s.gammaw, 1.e-3 / 1.8e-5);
    if (op.alfa != 1.) mat->SetAnisotropicConductivity(op.alfa * r.ks, r.ks); // k_x = alfa k_s, k_y = k_s
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false); // esqueleto rigido na percolacao, como no SEEP/W do artigo
    return mat;
}

// c', phi' (fmatprop) e k (fKMult) de cada ponto de integracao, do campo
void AplicaCampos(UnsatModel<TUnsatHP> &m, const Campo &f)
{
    m.ForEachIntegrationPoint([&f](const TPZVec<REAL> &x, TPZUnsaturatedMem &mem) {
        const std::vector<double> v = f(x[0], x[1]);
        const REAL phi = v[1] * M_PI / 180.;
        mem.m_elastoplastic_state.fmatprop = {v[0], phi, phi}; // {c, phi, psi}: fluxo associado
        mem.fKMult = v[2];
    });
}

// u: base engastada, laterais presas so em x
void InsertBC(UnsatModel<TUnsatHP> &m)
{
    m.AddDirichletBC(base, TUnsat::EDirichletU, 0, {0, 1});
    m.AddDirichletBC(left, TUnsat::EDirichletU, 0, {0});
    m.AddDirichletBC(right, TUnsat::EDirichletU, 0, {0});
}

// ---- fator de seguranca
// Modelo e analise do fator numa malha (a analise aponta para o modelo: apagar a analise antes)
struct Fator {
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a; // declarada depois: destruida antes do modelo
    REAL t = -1.;                               // o instante do estado carregado (o ultimo calculado nesta malha)
};

void NovoFator(Fator &f, TPZGeoMesh *g, const RainData &r, const Opcoes &op)
{
    f.a.reset();
    f.m = std::make_unique<UnsatModel<TUnsatHP>>();
    f.m->SetGeoMesh(g);
    InsertBC(*f.m);
    TUnsatHP *mat = MaterialChuva(r, op);
    if (op.gim) mat->SetGravityScalesSuction(false); // GIM: a succao nao e peso, fica fixa
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
    f.a->SetVerbose(false);
    f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
    f.a->SetSRMStart(op.srmF0);
}

// FS (ou Gamma) no modelo de f com as pressoes pw (nos da malha do transiente), do estado sem tensoes, com os campos.
// chi=ponto: o material tem chi = 1 e a pressao congelada num no com p_w < 0 e chi(x) p_w, com chi(x) = tan(phi^b) /
// tan(phi'(x)) e o phi' do campo no no (ou a media): a parcela da succao fica s tan(phi^b) (reduzida por F no SRM, como
// a de tan(phi')). Com p_w P1, o produto e interpolado linearmente entre os nos
REAL CalculaFator(Fator &f, const Campo *campos, const PointLocator &loc, const std::vector<REAL> &pw, const Opcoes &op)
{
    f.m->Material()->ResetMemory(); // estado sem tensoes (apaga tambem c, phi e k dos pontos)
    if (campos) AplicaCampos(*f.m, *campos);
    f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
        if (espaco != 1) return 0.;                    // p_a = 0
        const REAL p = loc.Interpolate(x[0], x[1], pw); // p_w do transiente
        if (!op.chiPonto || p >= 0.) return p;
        const REAL phi = (campos ? (*campos)(x[0], x[1])[1] : op.res.phiDeg) * M_PI / 180.;
        const REAL phib = op.phibRazao > 0. ? op.phibRazao * phi : op.phibDeg * M_PI / 180.;
        return std::tan(phib) / std::tan(phi) * p;
    });
    return op.gim ? f.a->GIMArcLength() : f.a->SRMArcLength();
}

// Refinamento adaptativo em h do zero: malha grossa e nref niveis (marcados: sqrt(J2(eps_p)) >= 0,1 do maior);
// niveis: o fator de cada nivel; f fica com o modelo do ultimo nivel
REAL FatorAdaptativo(Fator &f, const RainData &r, const Opcoes &op, const Campo *campos, const PointLocator &loc,
                     const std::vector<REAL> &pw, std::vector<REAL> &niveis)
{
    TPZGeoMesh *g = Malha(op.malhaF);
    niveis.clear();
    for (int k = 0;; k++) {
        NovoFator(f, g, r, op);
        const REAL fs = CalculaFator(f, campos, loc, pw, op);
        niveis.push_back(fs);
        if (k == op.nref) return fs;
        const std::vector<REAL> ind = f.a->PlasticIndicator();
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        std::set<int64_t> marcados;
        for (TPZGeoEl *gel : g->ElementVec())
            if (gel && gel->Dimension() == 2 && ind[gel->Index()] >= 0.1 * vmax) marcados.insert(gel->Index());
        std::vector<int64_t> pai;
        g = f.m->RefinedGeoMesh(marcados, pai);
    }
}

// ---- gradiente do fator pela projecao modal (o GradienteModal do Projects/Cho2010, com as pressoes congeladas)
// No colapso do SRM (c* = c'/F, tan(phi*) = tan(phi')/F), pelo limite superior com Mohr-Coulomb associado e o teorema do
// envelope, com as tensoes efetivas (m_sigma: as do modelo plastico, ja com a parcela chi s da succao) e
// f = (s1 - s3) + (s1 + s3) sin(phi) - 2 c cos(phi) (tracao positiva; eps_v = 2 mu sin(phi)):
//   dD/dc = eps_v cot(phi*),   dD/dphi = -(eps_v / (2 sin(phi*))) ((s1 + s3) cos(phi*) + 2 c* sin(phi*)),
//   dF/dc'(x) = dD/dc / Den,   dF/dtan(phi')(x) = dD/dphi cos^2(phi*) / Den,
//   Den = integral de (dD/dc c* + dD/dphi sin(phi*) cos(phi*)).
// O trabalho das pressoes congeladas nao depende de c', phi' nem de F: nao entra. No GIM (F = Gamma) os parametros sao os
// reais e Den = P_mr / Gamma (P_mr = integral de c' cot(phi') eps_v). Projecao nos modos: dF/dxi_j = soma_p w_p
// (dF/dln(c') dln(c')/dxi_j + dF/dln(phi') dln(phi')/dxi_j), com ln c' = mu_c + s_c Z_c e
// ln phi' = mu_f + s_f (rho Z_c + sqrt(1 - rho^2) Z_f) lineares em xi (dZ/dxi_k = sqrt(lambda_k) f_k(x): ProjectModes).
// O k so muda as pressoes, que o fator congela: os K termos de k saem 0
std::vector<double> GradienteModal(Fator &f, const RandomField &campos, int K, REAL F, bool gim)
{
    UnsatModel<TUnsatHP> &m = *f.m;
    TUnsatHP *mat = m.Material();
    std::vector<std::array<double, 2>> xs;
    std::vector<double> ac, af; // dF/dln(c') e dF/dln(phi') de cada ponto (ainda sem 1/Den)
    REAL den = 0.;
    for (TPZCompEl *cel : m.CompMesh()->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != mat->Id()) continue;
        TPZGeoEl *gel = cel->Reference();
        auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mf) continue;
        TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
        for (int i = 0; i < mf->NMeshes(); i++)
            if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
        const int o = mat->IntegrationRuleOrder(ordens);
        std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
        TPZManVector<int, 3> ordem(gel->Dimension(), o);
        regra->SetOrder(ordem);
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        for (int ip = 0; ip < regra->NPoints(); ip++) {
            TPZManVector<REAL, 3> qsi(gel->Dimension()), x(3, 0.);
            REAL w;
            regra->Point(ip, qsi, w);
            gel->X(qsi, x);
            TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
            REAL detjac;
            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
            const TPZUnsaturatedMem &me = mat->MemItem(mem[ip]);
            const auto &pr = me.m_elastoplastic_state.fmatprop;
            if (pr.size() < 2) continue;
            const TPZTensor<REAL> &ep = me.m_elastoplastic_state.m_eps_p;
            const REAL ev = ep.XX() + ep.YY() + ep.ZZ();
            if (!(ev > 0.)) continue; // ponto elastico (ou sem dilatancia): nao dissipa
            const TPZTensor<REAL> &sg = me.m_sigma;
            const REAL cen = 0.5 * (sg.XX() + sg.YY()), r = std::hypot(0.5 * (sg.XX() - sg.YY()), sg.XY());
            const REAL s1 = std::max({cen + r, cen - r, sg.ZZ()}), s3 = std::min({cen + r, cen - r, sg.ZZ()});
            const REAL c = pr[0], phi = pr[1];
            const REAL cs = gim ? c : c / F, phs = gim ? phi : std::atan(std::tan(phi) / F); // GIM: sem reducao
            const REAL dDc = ev / std::tan(phs);
            const REAL dDp = -(ev / (2. * std::sin(phs))) * ((s1 + s3) * std::cos(phs) + 2. * cs * std::sin(phs));
            const REAL pw = w * std::abs(detjac);
            den += gim ? pw * dDc * cs / F // P_mr / Gamma
                       : pw * (dDc * cs + dDp * std::sin(phs) * std::cos(phs));
            xs.push_back({x[0], x[1]});
            ac.push_back(pw * dDc * c);                                                         // dF/dln(c') * Den
            af.push_back(pw * dDp * std::cos(phs) * std::cos(phs) / (std::cos(phi) * std::cos(phi)) * phi); // dF/dln(phi')
        }
    }
    std::vector<double> g(3 * K, 0.);
    if (xs.empty() || !(std::abs(den) > 0.)) return g;
    const double sc = campos.SigmaLn(0), sf = campos.SigmaLn(1);
    const double rho = (campos.CorrA() == 0 && campos.CorrB() == 1) ? campos.Rho() : 0.;
    const double r2 = std::sqrt(1. - rho * rho);
    std::vector<std::vector<double>> pesos(2, std::vector<double>(xs.size()));
    for (size_t i = 0; i < xs.size(); i++) {
        pesos[0][i] = (ac[i] * sc + af[i] * sf * rho) / den; // termos de c': Z_c entra em c' e em phi'
        pesos[1][i] = af[i] * sf * r2 / den;                 // termos de phi'
    }
    const std::vector<std::vector<double>> p = campos.ProjectModes(xs, pesos, K);
    for (int k = 0; k < K; k++) g[k] = p[0][k], g[K + k] = p[1][k];
    return g; // os K de k ficam 0
}

// ---- descritores do instante (as colunas do relatorio do rebaixamento, SlopeRandomField/comparacao_ceron)
// No fim do caminho do fator, nos pontos de integracao da malha do fator:
//   banda: sqrt(J2(eps_p)) >= 0,1 do maior (a superficie de ruptura; c' e phi' do Spearman da Fig. 26)
//   massa: |du| >= 0,5 do maior, du = deslocamento desde o estado geostatico no SRM (no GIM, o total): o que desliza
//          (o k do Spearman)
//   campo: as medias do solo todo
struct Mecanismo {
    REAL bArea = 0., bProf = 0., bXmin = 0., bXmax = 0., bXc = 0., bYc = 0., bYmin = 0.; // banda: m2, m
    REAL bPw = 0., bC = 0., bPhi = 0., bK = 0.; // medias na banda: p_w (kPa), c' (kPa), phi' (graus), k / k_s
    REAL mArea = 0., mXc = 0., mYc = 0., mC = 0., mPhi = 0., mK = 0.; // massa
    REAL cMed = 0., phiMed = 0., kMed = 0.;                           // solo todo
    struct Ponto {
        REAL x, y, j, du, c, phi, k;
    };
    std::vector<Ponto> pontos; // os pontos com sqrt(J2) >= 1% do maior (mecanismo.csv das figuras)
    std::vector<std::array<REAL, 2>> banda;
};

Mecanismo DescreveMecanismo(Fator &f, const PointLocator &loc, const std::vector<REAL> &pw, const Resistencia &res)
{
    UnsatModel<TUnsatHP> &m = *f.m;
    TUnsatHP *mat = m.Material();
    const int iu = mat->VariableIndex("Displacement");
    struct P {
        TPZCompEl *cel;
        TPZManVector<REAL, 3> qsi;
        REAL x, y, w, j, c, phi, k, du;
    };
    std::vector<P> ps;
    for (TPZCompEl *cel : m.CompMesh()->ElementVec()) { // a regra da montagem, como o ForEachIntegrationPoint
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != mat->Id()) continue;
        TPZGeoEl *gel = cel->Reference();
        auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
        for (int i = 0; i < mf->NMeshes(); i++)
            if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
        const int o = mat->IntegrationRuleOrder(ordens);
        std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
        TPZManVector<int, 3> ordem(gel->Dimension(), o);
        regra->SetOrder(ordem);
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        for (int ip = 0; ip < regra->NPoints(); ip++) {
            P p;
            p.cel = cel;
            p.qsi.Resize(gel->Dimension());
            REAL w;
            regra->Point(ip, p.qsi, w);
            TPZManVector<REAL, 3> x(3, 0.);
            gel->X(p.qsi, x);
            TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
            REAL detjac;
            gel->Jacobian(p.qsi, jac, axes, detjac, jacinv);
            const TPZUnsaturatedMem &me = mat->MemItem(mem[ip]);
            TPZTensor<REAL> e(me.m_elastoplastic_state.m_eps_p); // cisalhamento de engenharia -> tensor
            e.XY() *= 0.5, e.XZ() *= 0.5, e.YZ() *= 0.5;
            const auto &pr = me.m_elastoplastic_state.fmatprop;
            const bool campo = pr.size() >= 2 && (pr[0] != 0. || pr[1] != 0.);
            p.x = x[0], p.y = x[1], p.w = w * std::abs(detjac), p.j = std::sqrt(std::max<REAL>(0., e.J2()));
            p.c = campo ? pr[0] : res.c, p.phi = campo ? pr[1] * 180. / M_PI : res.phiDeg, p.k = me.fKMult, p.du = 0.;
            ps.push_back(p);
        }
    }
    // deslocamento no fim do caminho (carregado) e no estado geostatico (SRM)
    auto deslocamentos = [&](std::vector<std::array<REAL, 2>> &u) {
        u.assign(ps.size(), {0., 0.});
        for (size_t i = 0; i < ps.size(); i++) {
            TPZManVector<STATE, 9> v(2, 0.);
            ps[i].cel->Solution(ps[i].qsi, iu, v);
            u[i] = {v[0], v[1]};
        }
    };
    std::vector<std::array<REAL, 2>> uf, ug(ps.size(), {0., 0.});
    deslocamentos(uf);
    if (f.a->WeightState().Rows() == f.a->PathState().Rows()) {
        f.a->LoadSolution(f.a->WeightState());
        deslocamentos(ug);
        f.a->LoadSolution(f.a->PathState());
    }
    REAL jmax = 0., dumax = 0., area = 0.;
    Mecanismo d;
    for (size_t i = 0; i < ps.size(); i++) {
        P &p = ps[i];
        p.du = std::hypot(uf[i][0] - ug[i][0], uf[i][1] - ug[i][1]);
        jmax = std::max(jmax, p.j), dumax = std::max(dumax, p.du);
        area += p.w;
        d.cMed += p.w * p.c, d.phiMed += p.w * p.phi, d.kMed += p.w * p.k;
    }
    d.cMed /= area, d.phiMed /= area, d.kMed /= area;
    if (!(jmax > 0.)) return d;
    d.bXmin = d.bYmin = 1.e300, d.bXmax = -1.e300;
    for (const P &p : ps) {
        if (p.j >= 0.01 * jmax) d.pontos.push_back({p.x, p.y, p.j, p.du, p.c, p.phi, p.k});
        if (p.j >= 0.1 * jmax) {
            d.bArea += p.w;
            d.bXc += p.w * p.x, d.bYc += p.w * p.y;
            d.bPw += p.w * loc.Interpolate(p.x, p.y, pw);
            d.bC += p.w * p.c, d.bPhi += p.w * p.phi, d.bK += p.w * p.k;
            d.bProf = std::max(d.bProf, Superficie(p.x) - p.y);
            d.bXmin = std::min(d.bXmin, p.x), d.bXmax = std::max(d.bXmax, p.x), d.bYmin = std::min(d.bYmin, p.y);
            d.banda.push_back({p.x, p.y});
        }
        if (dumax > 0. && p.du >= 0.5 * dumax) {
            d.mArea += p.w;
            d.mXc += p.w * p.x, d.mYc += p.w * p.y;
            d.mC += p.w * p.c, d.mPhi += p.w * p.phi, d.mK += p.w * p.k;
        }
    }
    for (REAL *v : {&d.bXc, &d.bYc, &d.bPw, &d.bC, &d.bPhi, &d.bK}) *v /= d.bArea;
    if (d.mArea > 0.)
        for (REAL *v : {&d.mXc, &d.mYc, &d.mC, &d.mPhi, &d.mK}) *v /= d.mArea;
    return d;
}

// cota do nivel d'agua na vertical x: o topo da zona saturada (p_w >= 0) ligada a base, subindo de 0,1 em 0,1 m
REAL NivelAgua(const PointLocator &loc, const std::vector<REAL> &pw, REAL x)
{
    const REAL topo = Superficie(x);
    REAL y = gTalude[0][1];
    while (y + 0.1 <= topo && loc.Interpolate(x, y + 0.1, pw) >= 0.) y += 0.1;
    return y;
}

// ---- o transiente: p_w nos nos da malha do transiente em cada instante do fator (calculado, ou lido de um arquivo)
struct Pressoes {
    std::map<REAL, std::vector<REAL>> pw; // instante -> p_w por no (os instantes depois de uma parada nao estao)
    std::map<REAL, size_t> face;          // nos na face de percolacao
    std::map<REAL, std::pair<REAL, REAL>> faixa; // menor e maior p_w (kPa), nos nos com p_w (RainTransient)
    std::map<REAL, double> s;             // tempo do transiente ate o instante (s)
    REAL parou = -1.;                     // >= 0: o transiente parou nesse instante (h)
    int64_t eq = 0;                       // equacoes da malha do transiente
    std::array<double, 3> permanente = {0., 0., 0.}; // q0 > 0: passos, ultima mudanca de p_w (kPa) e tempo (s)
};

// Os dados que definem o transiente (o que tem de ser igual para reaproveitar p_w): geometria, nivel, contornos, solo,
// chuva e malha. Os campos nao entram (com pressoes/grava_pressoes o k e fixo)
std::string ChaveTransiente(const Opcoes &op, const RainData &r)
{
    std::ostringstream c;
    c << std::setprecision(12) << "geometria=";
    for (int i = 0; i < 6; i++) c << (i ? "," : "") << gTalude[i][0] << "," << gTalude[i][1];
    c << ";hw=" << r.Hw << ";na_beta=" << op.betaNA << ";succao_max=" << r.suctionMax << ";esquerda="
      << (op.esqCarga ? "carga" : "fechada") << ";direita=" << (op.dirCarga ? "carga" : "fechada") << ";carga_esq="
      << op.cargaEsq << ";carga_dir=" << op.cargaDir << ";I=" << r.intensity << ";chuva=" << r.duration << ";dt=" << r.dt
      << ";chuva_face=" << (r.rainVertical ? "vertical" : "inclinada") << ";ks=" << r.ks << ";alfa=" << op.alfa;
    if (op.vg.usa) c << ";vg=" << op.vg.ts << "," << op.vg.tr << "," << op.vg.a << "," << op.vg.n;
    else c << ";fx=" << r.aFX << "," << r.nFX << "," << r.mFX << ";power=" << r.power;
    if (op.vgTab > 0.) c << ";vg_tab=" << op.vgTab << "," << NomeFechoVG(op.vgFecho);
    c << ";n=" << r.n << ";Kw=" << r.Kw << ";transiente=" << op.malhaT << ";uT=" << op.uT;
    if (op.kFredlund > 0.) c << ";kfredlund=" << op.kFredlund << "," << op.psiR;
    if (op.q0 > 0.) c << ";q0=" << op.q0;
    if (op.dtMin > 0.) c << ";dt_min=" << op.dtMin;
    if (!op.newton.empty())
        c << ";newton=" << op.newton[0] << "," << op.newton[1] << "," << (op.newton.size() > 2 ? op.newton[2] : 4.);
    if (op.faceTol > 0.) c << ";face_tol=" << op.faceTol;
    return c.str();
}

// newton_trace= e face_tol= num transiente (so quando pedidos). O trace recebe as tolerancias postas no TPZUnsatAnalysis:
// as do newton= ou os padroes dele (1e-7, 1e-6, 40, 25)
void OpcoesTransiente(RainTransient &chuva, const Opcoes &op, const std::string &rotulo)
{
    if (!op.newtonTrace.empty())
        chuva.SetNewtonTrace(op.newtonTrace, 1.e-7, op.newton.empty() ? 1.e-6 : op.newton[0],
                             op.newton.empty() ? 40 : int(op.newton[1]), 25, rotulo);
    if (op.faceTol > 0.) chuva.SetFaceTolerance(op.faceTol);
}

// O transiente na malha de mT (com o material e os campos ja postos): ate o ultimo instante, com p_w em cada um; rotulo:
// o transiente no arquivo do newton_trace=
Pressoes CalculaTransiente(UnsatModel<TUnsatHP> &mT, const RainData &r, const Opcoes &op, const std::string &rotulo = "")
{
    Pressoes P;
    TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
    TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
    strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
    aT.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    aT.SetSolver(solver);
    if (!op.newton.empty()) aT.SetTolerances(1.e-7, op.newton[0], int(op.newton[1]), 25);
    std::set<int> laterais;
    std::map<int, REAL> cargas;
    if (op.esqCarga) laterais.insert(left);
    if (op.dirCarga) laterais.insert(right);
    if (op.esqCarga && !std::isnan(op.cargaEsq)) cargas[left] = op.cargaEsq;
    if (op.dirCarga && !std::isnan(op.cargaDir)) cargas[right] = op.cargaDir;
    RainTransient chuva(aT, r, op.res.gammaw, laterais, {topleft, face, topright}, cargas);
    OpcoesTransiente(chuva, op, rotulo);
    chuva.InitialState();
    P.eq = mT.CompMesh()->NEquations();
    if (op.q0 > 0.) { // estado inicial permanente com a infiltracao q0
        const auto tp = Clock::now();
        REAL dp = 0.;
        if (!chuva.SteadyState(op.q0 / 3.6e6, r.dt, 1.e-6, &dp)) {
            std::printf("  o permanente com q0 = %g mm/h nao estacionou (%lld passos, mudanca %.2e kPa)\n", op.q0,
                        (long long)chuva.SteadySteps(), dp);
            P.parou = 0.;
            return P;
        }
        P.permanente = {double(chuva.SteadySteps()), dp, Seg(tp)};
        std::printf("  permanente com q0 = %g mm/h: %lld passos, ultima mudanca %.1e kPa, %.1f s\n", op.q0,
                    (long long)chuva.SteadySteps(), dp, Seg(tp));
    }
    double sT = 0.;
    for (REAL t : op.tempos) {
        const auto tt = Clock::now();
        if (!chuva.AdvanceTo(t)) {
            P.parou = chuva.Time();
            std::printf("  o transiente parou em %.2f h: %s\n", P.parou, chuva.LastFailure().c_str());
            break;
        }
        sT += Seg(tt);
        P.pw[t] = chuva.NodalWaterPressure();
        P.face[t] = chuva.FaceNodes();
        P.faixa[t] = chuva.PressureRange();
        P.s[t] = sT;
        if (!op.newtonTrace.empty()) { // newton_trace=: o progresso na tela, instante a instante
            std::printf("  transiente em %g h: %lld passos, %lld iteracoes de Newton, %zu nos na face (%.0f s)\n", t,
                        (long long)chuva.Steps(), (long long)chuva.NewtonIterations(), chuva.FaceNodes(), sT);
            std::fflush(stdout);
        }
    }
    return P;
}

// arquivo do grava_pressoes / pressoes: a chave do transiente, os nos (x y) e, por instante, 't <t> <face> <s> <pmin>
// <pmax>' e p_w por no
void GravaPressoes(const std::string &arq, const std::string &chave, TPZGeoMesh *g, const Pressoes &P)
{
    const std::string tmp = arq + ".tmp";
    {
        std::ofstream f(tmp);
        f << "# SlopeRainMonteCarlo: p_w (kPa) nos nos da malha do transiente nos instantes do fator (grava_pressoes=); "
             "le-se com pressoes=\n"
          << "chave " << chave << "\n"
          << "nos " << g->NNodes() << " instantes " << P.pw.size() << " eq " << P.eq << "\n"
          << std::setprecision(17);
        for (int64_t i = 0; i < g->NNodes(); i++) f << g->NodeVec()[i].Coord(0) << " " << g->NodeVec()[i].Coord(1) << "\n";
        for (const auto &[t, pw] : P.pw) {
            f << "t " << t << " " << P.face.at(t) << " " << P.s.at(t) << " " << P.faixa.at(t).first << " "
              << P.faixa.at(t).second << "\n";
            for (REAL v : pw) f << v << "\n";
        }
        if (!f) Erro("nao consegui gravar " + tmp);
    }
    if (std::rename(tmp.c_str(), arq.c_str()) != 0) Erro("nao consegui renomear " + tmp + " para " + arq);
}

Pressoes LePressoes(const std::string &arq, const std::string &chave, TPZGeoMesh *g, const std::vector<REAL> &tempos)
{
    std::ifstream f(arq);
    if (!f) Erro("pressoes: nao abre " + arq);
    std::string linha, rot, ch;
    std::getline(f, linha); // comentario
    std::getline(f, linha);
    std::istringstream lc(linha);
    lc >> rot >> ch;
    if (rot != "chave") Erro("pressoes: " + arq + " nao e um arquivo do grava_pressoes");
    if (ch != chave)
        Erro("pressoes: o transiente de " + arq + " e outro:\n  arquivo: " + ch + "\n  pedido:  " + chave);
    int64_t nn = 0;
    size_t nt = 0;
    Pressoes P;
    f >> rot >> nn >> rot >> nt >> rot >> P.eq;
    if (!f || nn != g->NNodes()) Erro("pressoes: " + arq + " tem " + std::to_string(nn) + " nos e a malha do transiente " +
                                      std::to_string(g->NNodes()));
    for (int64_t i = 0; i < nn; i++) {
        REAL x, y;
        f >> x >> y;
        if (!f || std::abs(x - g->NodeVec()[i].Coord(0)) > 1e-9 || std::abs(y - g->NodeVec()[i].Coord(1)) > 1e-9)
            Erro("pressoes: os nos de " + arq + " nao sao os da malha do transiente (no " + std::to_string(i) + ")");
    }
    std::map<REAL, std::vector<REAL>> pw;
    std::map<REAL, size_t> face;
    std::map<REAL, std::pair<REAL, REAL>> faixa;
    std::map<REAL, double> s;
    for (size_t k = 0; k < nt; k++) {
        REAL t, pmin, pmax;
        size_t nf;
        double sk;
        f >> rot >> t >> nf >> sk >> pmin >> pmax;
        if (!f || rot != "t") Erro("pressoes: " + arq + " incompleto (instante " + std::to_string(k) + ")");
        std::vector<REAL> v(nn);
        for (REAL &x : v) f >> x;
        if (!f) Erro("pressoes: " + arq + " incompleto no instante " + std::to_string(t));
        pw[t] = std::move(v), face[t] = nf, faixa[t] = {pmin, pmax}, s[t] = sk;
    }
    for (REAL t : tempos) { // os instantes pedidos, com a tolerancia do texto
        auto it = std::find_if(pw.begin(), pw.end(), [t](const auto &e) { return std::abs(e.first - t) < 1e-9; });
        if (it == pw.end()) Erro("pressoes: " + arq + " nao tem o instante " + std::to_string(t) + " h");
        P.pw[t] = it->second, P.face[t] = face[it->first], P.faixa[t] = faixa[it->first], P.s[t] = s[it->first];
    }
    return P;
}

// uma linha do arquivo: realizacao e instante
struct Linha {
    int real = 0;
    REAL t = 0., fs = 0.;
    std::string situacao = "ok"; // ok, transiente (parou antes deste instante) ou fator (sem fator)
    bool adaptativo = false;
    std::string niveis;
    int64_t eqF = 0, eqT = 0;
    double sF = 0., sT = 0.;
    REAL pwMin = 0., pwMax = 0.;
    size_t face = 0;
    REAL nivel30 = 0., nivel35 = 0., nivel40 = 0.;
    Mecanismo mec;
};

// colunas: as do mc_resumo.csv do rebaixamento onde existem (FS_min, banda_*, massa_*, *_med, cov_*, Lx, Ly, alfa); os
// niveis sao na crista, no meio da face e no pe (x = 30, 35 e 40 m no Rahardjo e no caso A)
const char *kCabecalho =
    "real,t_h,FS_min,situacao,adaptativo,niveis,eq_fator,s_fator,eq_transiente,s_transiente,pw_min,pw_max,nos_face,"
    "nivel_x30,nivel_x35,nivel_x40,banda_area_m2,banda_prof,banda_xmin,banda_xmax,banda_xc,banda_yc,banda_ymin,banda_pw,"
    "banda_c,banda_phi_deg,banda_k_mult,massa_area_m2,massa_xc,massa_yc,massa_c,massa_phi_deg,massa_k_mult,c_med,phi_med,"
    "k_mult_med,solo,I,chuva,hw,metodo,cov_c,cov_phi,cov_k,Lx,Ly,alfa,KL_variancia,log_peso,proj_c,proj_phi,proj_k\n";

void Grava(std::ofstream &csv, const Linha &l, const Opcoes &op, REAL klVariancia, double logPeso,
           const std::vector<double> &proj)
{
    const Mecanismo &d = l.mec;
    std::ostringstream s;
    s << l.real << "," << l.t << "," << std::setprecision(7) << l.fs << "," << l.situacao << "," << int(l.adaptativo)
      << "," << l.niveis << "," << l.eqF << "," << std::setprecision(4) << l.sF << "," << l.eqT << "," << l.sT << ","
      << std::setprecision(6) << l.pwMin << "," << l.pwMax << "," << l.face << "," << l.nivel30 << "," << l.nivel35
      << "," << l.nivel40 << "," << d.bArea << "," << d.bProf << "," << d.bXmin << "," << d.bXmax << "," << d.bXc << ","
      << d.bYc << "," << d.bYmin << "," << d.bPw << "," << d.bC << "," << d.bPhi << "," << d.bK << "," << d.mArea << ","
      << d.mXc << "," << d.mYc << "," << d.mC << "," << d.mPhi << "," << d.mK << "," << d.cMed << "," << d.phiMed << ","
      << d.kMed << "," << op.solo.nome << "," << op.I << "," << op.chuva << "," << op.hw << ","
      << (op.gim ? "gim" : "srm") << "," << op.campos.covC << "," << op.campos.covPhi << "," << op.campos.covK << ","
      << op.campos.Lx << "," << op.campos.Ly << "," << op.alfa << "," << klVariancia << "," << std::setprecision(10)
      << logPeso << "," << proj[0] << "," << proj[1] << "," << proj[2] << "\n";
    csv << s.str();
    csv.flush();
}

// parametros.txt (chave=valor), como no rebaixamento
void GravaParametros(const std::string &arq, const Opcoes &op, const RandomField &campos, const RainData &r)
{
    std::ofstream p(arq);
    auto lista = [](const std::vector<REAL> &v) {
        std::string s;
        for (size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + std::to_string(v[i]).substr(0, 6);
        return s;
    };
    p << "SOLO=" << op.solo.nome << "\nA_FX=" << op.solo.a << "\nKS=" << op.solo.ks << "\nI=" << op.I
      << "\nCHUVA=" << op.chuva << "\nHW=" << op.hw << "\nMETODO=" << (op.gim ? "gim" : "srm")
      << "\nCOV_C=" << op.campos.covC << "\nCOV_PHI=" << op.campos.covPhi << "\nCOV_K=" << op.campos.covK
      << "\nLX=" << op.campos.Lx << "\nLY=" << op.campos.Ly << "\nALFA=" << op.alfa << "\nKL_M=" << campos.M()
      << "\nKL_VARIANCIA=" << campos.VarianceFraction() << "\nSEMENTE=" << op.campos.semente
      << "\nMALHA_TRANSIENTE=" << op.malhaT << "\nUT=" << op.uT << "\nMALHA_FATOR=" << op.malhaF
      << "\nNREF=" << op.nref << "\nTEMPOS=" << lista(op.tempos) << "\nADAPTA=" << lista(op.adapta)
      << "\nAI=" << (!op.deslocXi.empty() ? std::string("dado:") + op.deslocXi
                     : op.ai.piloto.empty() ? std::string("nenhuma") : op.ai.tipo + "_piloto") << "\n";
    // o problema e os campos (o transiente: a mesma chave do pressoes=)
    p << "TRANSIENTE=" << ChaveTransiente(op, r) << "\nCHI=" << (op.chiSe ? "se" : "const") << "\nC=" << op.res.c
      << "\nPHI=" << op.res.phiDeg << "\nGAMMA=" << op.res.gamma << "\nE=" << op.res.E << "\nNU=" << op.res.nu
      << "\nR_CPHI=" << op.campos.rcphi << "\nNYSTROM=" << (op.nystrom.empty() ? "nao" : op.nystrom)
      << "\nNUCLEO=" << (campos.Nystrom() ? campos.KernelDescription() : "exp analitico")
      << "\nPRESSOES=" << (op.pressoes.empty() ? "calculadas" : op.pressoes) << "\nSRM_F0=" << op.srmF0 << "\n";
    // casos B e C (so quando usados)
    if (op.chiPonto) p << "CHI=ponto\nPHIB=" << op.phibDeg << "\nPHIB_RAZAO=" << op.phibRazao << "\n";
    if (!op.fx.empty()) p << "FX=" << lista(op.fx) << "\n";
    if (op.kFredlund > 0.) p << "KFREDLUND=" << op.kFredlund << "\nPSI_R=" << op.psiR << "\n";
    if (op.vgTab > 0.) p << "VG_TAB=" << op.vgTab << "," << NomeFechoVG(op.vgFecho) << "\n";
    if (op.q0 > 0.) p << "Q0=" << op.q0 << "\nQ0_COV=" << op.q0Cov << "\n";
    if (op.chuvaAleatoria) {
        p << "CHUVA_ALEATORIA=1\nGPD=";
        for (size_t i = 0; i < op.gpd.size(); i++) p << (i ? "," : "") << op.gpd[i];
        p << "\n";
    }
    if (!op.crpc.empty()) p << "CRPC=" << lista(op.crpc) << "\nCRPC_DMAX=" << op.crpcDmax << "\nCRPC_PASSO=" << op.crpcPasso << "\n";
    if (!op.ai.piloto.empty())
        p << "PILOTO=" << op.ai.piloto << "\nPILOTO_N=" << op.ai.n << "\nPILOTO_SEMENTE=" << op.ai.semente
          << "\nELITE=" << op.ai.elite << "\nMODOS=" << op.ai.modos << "\nBANDA=" << op.banda << "\nBANDA_T=" << op.bandaT
          << "\n";
}

// mapa da faixa do mecanismo: celulas de 1 m no retangulo do talude, por instante
struct Mapa {
    int nx = 70, ny = 30;
    REAL x0 = 0., y0 = 0.;
    std::map<REAL, std::vector<int>> n, nFalha;
    int nReal = 0;
    Mapa()
    {
        const auto &p = gTalude;
        x0 = p[0][0], y0 = p[0][1];
        nx = std::max(1, int(std::ceil(p[1][0] - x0 - 1e-9))), ny = std::max(1, int(std::ceil(p[4][1] - y0 - 1e-9)));
    }
    void Soma(REAL t, const Mecanismo &d, bool falha)
    {
        std::vector<int> &v = n[t], &vf = nFalha[t];
        v.resize(nx * ny, 0), vf.resize(nx * ny, 0);
        std::set<int> celulas;
        for (const auto &p : d.banda)
            celulas.insert(std::clamp(int(p[1] - y0), 0, ny - 1) * nx + std::clamp(int(p[0] - x0), 0, nx - 1));
        for (int c : celulas) v[c]++, vf[c] += falha;
    }
    void Grava(const std::string &arq) const
    {
        std::ofstream f(arq);
        f << "t_h,x,y,n,n_falha,n_real\n";
        for (const auto &[t, v] : n)
            for (int c = 0; c < nx * ny; c++)
                if (v[c]) f << t << "," << x0 + c % nx + 0.5 << "," << y0 + c / nx + 0.5 << "," << v[c] << ","
                            << nFalha.at(t)[c] << "," << nReal << "\n";
    }
};

// os 3K numeros de um arquivo (xi= e desloc_xi=): os K de c', de phi' e de k
std::vector<std::vector<double>> LeXi(const std::string &arq, int M, const char *opcao)
{
    std::ifstream fx(arq);
    std::vector<double> v;
    for (double x; fx >> x;) v.push_back(x);
    if (v.empty() || v.size() % 3 || int(v.size() / 3) > M)
        Erro(std::string(opcao) + ": " + arq + " tem " + std::to_string(v.size()) + " numeros; precisa de 3K, com K <= M = " +
             std::to_string(M));
    const size_t K = v.size() / 3;
    std::vector<std::vector<double>> xi;
    for (int c = 0; c < 3; c++) xi.push_back(std::vector<double>(v.begin() + c * K, v.begin() + (c + 1) * K));
    return xi;
}

// ---- a chuva aleatoria e o estado inicial aleatorio (caso C): numeros da realizacao n, independentes dos campos
// (splitmix64 de (semente, n): a realizacao n tem a mesma chuva em todos os casos com a mesma semente)
struct Uniformes {
    uint64_t s;
    Uniformes(uint64_t semente, int64_t n) : s(semente * 0x9E3779B97F4A7C15ull + uint64_t(n) * 0xD1B54A32D192ED03ull + 0x632BE59BD9B4E019ull) {}
    double Proximo() // em (0, 1)
    {
        s += 0x9E3779B97F4A7C15ull;
        uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        z ^= z >> 31;
        return (double(z >> 11) + 0.5) * 0x1.0p-53;
    }
};

struct ChuvaSorteada {
    REAL I = 0., D = 0., q0 = 0., u = 0., v = 0., z = 0.; // mm/h, h, mm/h; u = F_I(I), v = F_D(D); z: o normal de q0
};

// inversa da Pareto generalizada: x = mu + sigma ((1 - u)^(-xi) - 1) / xi (xi = 0: mu - sigma ln(1 - u))
REAL InversaGPD(REAL u, REAL xi, REAL mu, REAL sigma)
{
    return std::abs(xi) < 1e-12 ? mu - sigma * std::log1p(-u) : mu + sigma * (std::pow(1. - u, -xi) - 1.) / xi;
}

// I e D da Tabela 2 de Liu e Wang (2023): u ~ U(0, 1) e v pela condicional da copula de Frank (Nelsen 2006, 4.10),
// v = -ln(1 + w (e^-omega - 1) / (w + (1 - w) e^(-omega u))) / omega com w ~ U(0, 1); I = F_I^-1(u), D = F_D^-1(v). q0:
// normal com media q0 e CoV q0_cov (Box-Muller), truncado em 0
ChuvaSorteada SorteiaChuva(const Opcoes &op, int64_t n)
{
    Uniformes U(op.campos.semente, n);
    const auto &g = op.gpd;
    ChuvaSorteada c;
    const double u = U.Proximo(), w = U.Proximo(), a = U.Proximo(), b = U.Proximo();
    const double om = g[6];
    c.u = u;
    c.v = -std::log1p(w * std::expm1(-om) / (w + (1. - w) * std::exp(-om * u))) / om;
    c.v = std::clamp(c.v, REAL(1e-16), REAL(1. - 1e-16));
    c.I = InversaGPD(c.u, g[0], g[1], g[2]);
    c.D = InversaGPD(c.v, g[3], g[4], g[5]);
    c.z = std::sqrt(-2. * std::log(a)) * std::cos(2. * M_PI * b);
    c.q0 = std::max(REAL(0.), op.q0 * (1. + op.q0Cov * c.z));
    return c;
}

// ---- a CRPC de uma realizacao (opcao crpc=; Liu e Wang 2023, Secao 3.1): para cada I (decrescente), o estado inicial
// permanente com q0, chuva I constante ate crpc_dmax e o primeiro instante com FS < 1 (D_C). Linhas no csv:
// real,q0_mmh,FS0,I,DC_h,falha,t_lo,FS_lo,t_hi,FS_hi,n_srm,s_transiente,s_fator,passos_permanente,mudanca_permanente
void CurvaCritica(int n, const Opcoes &op, const RainData &r0, const Campo *cf, std::ofstream &csv)
{
    REAL fs0 = -1.;
    const std::vector<REAL> checa = {6., 12., 24., 48., 96., 150.};
    for (REAL I : op.crpc) {
        RainData r = r0;
        r.intensity = I, r.duration = op.crpcDmax;
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        InsertBC(mT);
        mT.SetMaterial(MaterialChuva(r, op));
        mT.CreateCompMesh({{2, op.uT}, {1, 1}, {1, 1}});
        if (cf) AplicaCampos(mT, *cf);
        TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
        TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
        strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
        aT.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        aT.SetSolver(solver);
        if (!op.newton.empty()) aT.SetTolerances(1.e-7, op.newton[0], int(op.newton[1]), 25);
        std::set<int> laterais;
        std::map<int, REAL> cargas;
        if (op.esqCarga) laterais.insert(left);
        if (op.dirCarga) laterais.insert(right);
        if (op.esqCarga && !std::isnan(op.cargaEsq)) cargas[left] = op.cargaEsq;
        if (op.dirCarga && !std::isnan(op.cargaDir)) cargas[right] = op.cargaDir;
        RainTransient chuva(aT, r, op.res.gammaw, laterais, {topleft, face, topright}, cargas);
        char rotulo[64];
        std::snprintf(rotulo, sizeof rotulo, "realizacao %d, I = %g mm/h", n, I);
        OpcoesTransiente(chuva, op, rotulo);
        chuva.InitialState();
        const auto tt = Clock::now();
        REAL dp = 0.;
        if (op.q0 > 0. && !chuva.SteadyState(op.q0 / 3.6e6, r.dt, 1.e-6, &dp)) {
            std::printf("  realizacao %d: o permanente nao estacionou\n", n);
            csv << n << "," << op.q0 << ",-1," << I << ",-1,erro_permanente,0,0,0,0,0,0,0," << chuva.SteadySteps() << ","
                << dp << "\n";
            return;
        }
        const PointLocator loc(mT.GeoMesh());
        double sT = Seg(tt), sF = 0.;
        int nSrm = 0;
        std::map<REAL, std::vector<REAL>> pw; // instante (h) -> p_w nos nos
        std::map<REAL, REAL> fs;              // os FS calculados
        pw[0.] = chuva.NodalWaterPressure();
        auto fator = [&](REAL t) {
            const auto tf = Clock::now();
            Fator f;
            std::vector<REAL> niveis;
            const REAL v = FatorAdaptativo(f, r, op, cf, loc, pw.at(t), niveis);
            sF += Seg(tf), nSrm++;
            return fs[t] = v;
        };
        if (fs0 < 0.) {
            fs0 = fator(0.);
            if (!(fs0 >= 1.)) { // descartada (Liu e Wang 2023: as amostras com FS0 < 1 saem da conta)
                csv << n << "," << op.q0 << "," << fs0 << "," << I << ",0,descartada,0," << fs0 << ",0," << fs0 << ","
                    << nSrm << "," << sT << "," << sF << "," << chuva.SteadySteps() << "," << dp << "\n";
                csv.flush();
                std::printf("  realizacao %d: FS0 = %.4f < 1, descartada\n", n, fs0);
                return;
            }
        }
        fs[0.] = fs0;
        // avanca de crpc_passo em crpc_passo h, guardando p_w; o FS nos instantes de checa e em crpc_dmax
        REAL lo = 0., hi = -1., t = 0.;
        bool parou = false;
        std::vector<REAL> pontos;
        for (REAL c : checa)
            if (c < op.crpcDmax) pontos.push_back(c);
        pontos.push_back(op.crpcDmax);
        size_t ip = 0;
        while (t < op.crpcDmax - 1e-9) {
            const REAL prox = std::min(op.crpcDmax, t + op.crpcPasso);
            const auto ta = Clock::now();
            if (!chuva.AdvanceTo(prox)) {
                parou = true;
                break;
            }
            sT += Seg(ta);
            t = prox;
            pw[t] = chuva.NodalWaterPressure();
            if (ip < pontos.size() && t >= pontos[ip] - 1e-9) {
                ip++;
                if (fator(t) < 1.) {
                    hi = t;
                    break;
                }
                lo = t;
            }
        }
        if (parou) {
            // t_lo e FS_lo: o ultimo instante com FS calculado (>= 1); t_hi: onde o transiente parou
            csv << n << "," << op.q0 << "," << fs0 << "," << I << ",-1,erro_transiente," << lo << "," << fs.at(lo) << ","
                << chuva.Time() << ",-1," << nSrm << ","
                << sT << "," << sF << "," << chuva.SteadySteps() << "," << dp << "\n";
            csv.flush();
            std::printf("  realizacao %d, I = %g: o transiente parou em %.2f h: %s\n", n, I, chuva.Time(),
                        chuva.LastFailure().c_str());
            continue;
        }
        if (hi < 0.) { // nao caiu ate crpc_dmax: as I menores tambem nao (nao sao feitas)
            csv << n << "," << op.q0 << "," << fs0 << "," << I << ",-1,nao," << lo << "," << fs.at(lo) << ",-1,-1," << nSrm
                << "," << sT << "," << sF << "," << chuva.SteadySteps() << "," << dp << "\n";
            csv.flush();
            std::printf("  realizacao %d, I = %g mm/h: FS(%g h) = %.4f >= 1 (%d SRM, %.0f s + %.0f s)\n", n, I, lo,
                        fs.at(lo), nSrm, sT, sF);
            return;
        }
        // bissecao nos instantes guardados entre lo (FS >= 1) e hi (FS < 1)
        auto instantes = [&](REAL a, REAL b) {
            std::vector<REAL> v;
            for (const auto &[tk, x] : pw)
                if (tk > a + 1e-9 && tk < b - 1e-9) v.push_back(tk);
            return v;
        };
        for (std::vector<REAL> v = instantes(lo, hi); !v.empty(); v = instantes(lo, hi)) {
            const REAL m = v[v.size() / 2];
            if (fator(m) < 1.) hi = m;
            else lo = m;
        }
        const REAL flo = fs.at(lo), fhi = fs.at(hi);
        const REAL dc = lo + (flo - 1.) / (flo - fhi) * (hi - lo);
        csv << n << "," << op.q0 << "," << fs0 << "," << I << "," << dc << ",sim," << lo << "," << flo << "," << hi << ","
            << fhi << "," << nSrm << "," << sT << "," << sF << "," << chuva.SteadySteps() << "," << dp << "\n";
        csv.flush();
        std::printf("  realizacao %d, I = %g mm/h: D_C = %.2f h (FS %.4f em %g h, %.4f em %g h; %d SRM, %.0f s + %.0f s)\n",
                    n, I, dc, flo, lo, fhi, hi, nSrm, sT, sF);
        std::fflush(stdout);
    }
}

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    const Resistencia &res = op.res;
    const auto &p = gTalude;
    RainData r; // Fredlund-Xing n = m = 1, theta_s = 0,45, k_w = k_s Theta^4, phi^b = 26, succao ate 75 kPa
    r.aFX = op.solo.a, r.ks = op.solo.ks;
    r.intensity = op.I, r.duration = op.chuva, r.Hw = op.hw;
    r.xToe = p[3][0], r.yToe = p[3][1];
    if (op.betaNA != 7.) r.betaW = op.betaNA * M_PI / 180.; // o padrao do RainData e 7 graus
    r.suctionMax = op.succaoMax;
    r.rainVertical = op.chuvaVertical;
    if (op.vg.usa) r.n = op.vg.ts; // porosidade = theta_s
    if (!op.fx.empty()) r.aFX = op.fx[0], r.nFX = op.fx[1], r.mFX = op.fx[2], r.n = op.fx[3];
    if (op.phibDeg != 26.) r.phib = op.phibDeg * M_PI / 180.;
    if (op.dt != 0.1) r.dt = op.dt;
    r.dtMin = op.dtMin;
    if (op.newton.size() > 2) r.lineSearch = int(op.newton[2]);
    if (op.kFredlund > 0.) {
        gTabelaK = TabelaFredlund(r.aFX, r.nFX, r.mFX, r.n, op.psiR, op.kFredlund);
        TPZUnsaturatedHydraulics h;
        h.fA = r.aFX, h.fN = r.nFX, h.fM = r.mFX, h.fKrw = TPZUnsaturatedHydraulics::ETabulated;
        h.fTabS = gTabelaK.S, h.fTabLnK = gTabelaK.lnK;
        std::printf("  k_r de Fredlund (psi_aev = %g kPa, psi_r = %g kPa, %zu pontos):", op.kFredlund, op.psiR,
                    gTabelaK.S.size());
        for (REAL psi : {1., 2.6, 5., 10., 15., 20., 50., 100.}) {
            const auto e = h.Evaluate(-psi, 0.);
            std::printf(" k_r(%g) = %.4g;", psi, e.krw);
        }
        std::printf("\n");
    }
    if (op.vgTab > 0.) { // vg_tab=: o k_r de Mualem tabelado ate S_e = 1 - eps, conferido no material contra o EMualem
        const REAL swr = op.vg.tr / op.vg.ts;
        gTabelaVG = MualemTabelado(op.vg.n, swr, op.vg.a, op.vgTab, op.vgFecho);
        TPZUnsaturatedHydraulics h, hm;
        h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        h.fAVG = op.vg.a, h.fBVG = op.vg.n, h.fSwr = swr;
        hm = h;
        hm.fKrw = TPZUnsaturatedHydraulics::EMualem;
        UsaTabela(h, gTabelaVG);
        std::printf("  k_r de Mualem tabelado (fecho %s): eps = %.3g, psi_c = %.4g kPa, k_c = %.4f, %zu pontos, erro maximo "
                    "%.4f%% em [psi_c, 1e3 kPa];", NomeFechoVG(op.vgFecho), gTabelaVG.eps, gTabelaVG.psiC, gTabelaVG.kC,
                    gTabelaVG.S.size(), 100. * gTabelaVG.erroMax);
        for (REAL psi : {0.01, 0.1, 1., 10.})
            std::printf(" k_r(%g) = %.4f (Mualem %.4f);", psi, h.Evaluate(-psi, 0.).krw, hm.Evaluate(-psi, 0.).krw);
        std::printf("\n");
    }
    if (!op.newtonTrace.empty()) // um bloco por transiente que parar (RainTransient::SetNewtonTrace)
        std::ofstream(op.newtonTrace) << "# newton_trace do SlopeRainMonteCarlo: um bloco por transiente que parou "
                                         "(sem blocos: nenhum parou)\n";
    const Campos &cp = op.campos;
    RandomField campos(p[0][0], p[0][1], p[1][0] - p[0][0], p[4][1] - p[0][1], cp.Lx, cp.Ly, cp.M,
                       {{res.c, cp.covC}, {res.phiDeg, cp.covPhi}, {1., cp.covK}}); // c' (kPa), phi' (graus), k / k_s
    if (cp.rcphi != 0.) campos.SetCrossCorrelation(0, 1, cp.rcphi);
    if (!op.nystrom.empty()) {
        const auto tn = Clock::now();
        campos.UseNystrom(op.nystrom, cp.M, true);
        if (campos.GaussianKernel())
            std::printf("  aviso: a base de Nystrom tem o nucleo gaussiano (%s): s= nao vale\n",
                        campos.KernelDescription().c_str());
        else if (std::abs(campos.KernelLx() - cp.Lx) > 1e-9 * cp.Lx || std::abs(campos.KernelLy() - cp.Ly) > 1e-9 * cp.Ly)
            Erro("nystrom: a base " + op.nystrom + " tem o nucleo " + campos.KernelDescription() + ", e o caso pede lh = " +
                 std::to_string(cp.Lx) + " e lv = " + std::to_string(cp.Ly) + " (s=)");
        std::printf("  base de Nystrom %s: %d termos, variancia retida %.4f, nucleo %s, %s (%.1f s)\n", op.nystrom.c_str(),
                    campos.M(), campos.VarianceFraction(), campos.KernelDescription().c_str(),
                    campos.FastNystrom() ? "grade regular (avaliacao rapida)" : "soma direta", Seg(tn));
    }
    const Campo campoAleatorio = [&campos](double x, double y) { return campos.At(x, y); };
    const Campo campoFixo = [v = op.valores](double, double) { return std::vector<double>(v.begin(), v.end()); };
    std::printf("Monte Carlo com chuva: talude (%g;%g) (%g;%g) (%g;%g) (%g;%g) (%g;%g) (%g;%g); solo %s (k_s %g m/s, %s, chi "
                "%s), I = %g mm/h por %g h (%s), nivel %g m abaixo do pe a %g graus (succao ate %g kPa), laterais %s/%s | "
                "%s | realizacoes %d..%d%s\n  campos: c' %g kPa (CoV %g), phi' %g graus (CoV %g), k (CoV %g), r(c, phi) = "
                "%g; Lx = %g m, Ly = %g m; alfa = %g; KL: %d termos, %.1f%% da variancia\n  malhas: transiente %s (u P%d), "
                "fator %s + %d niveis\n",
                p[0][0], p[0][1], p[1][0], p[1][1], p[2][0], p[2][1], p[3][0], p[3][1], p[4][0], p[4][1], p[5][0], p[5][1],
                op.solo.nome.c_str(), op.solo.ks,
                op.vg.usa ? (op.vgTab > 0. ? "van Genuchten-Mualem tabelado" : "van Genuchten-Mualem") : "Fredlund-Xing",
                op.chiSe ? "S_e" : op.chiPonto ? "por ponto" : "constante", op.I, op.chuva, op.chuvaVertical ? "vertical" : "por metro inclinado",
                op.hw, op.betaNA, op.succaoMax, op.esqCarga ? "carga" : "fechada", op.dirCarga ? "carga" : "fechada",
                op.gim ? "GIM" : "SRM", op.primeira, op.primeira + op.n - 1, op.det ? " e a deterministica" : "", res.c,
                cp.covC, res.phiDeg, cp.covPhi, cp.covK, cp.rcphi, cp.Lx, cp.Ly, op.alfa, campos.M(),
                100. * campos.VarianceFraction(), op.malhaT.c_str(), op.uT, op.malhaF.c_str(), op.nref);

    // o transiente reaproveitado: so com k fixo (o LeOpcoes confere); lido do arquivo ou calculado uma vez (abaixo)
    const bool kFixo = cp.covK == 0. && (op.valores.empty() || op.valores[2] == 1.) && op.q0Cov == 0. &&
                       !op.chuvaAleatoria;
    if (op.q0 > 0.)
        std::printf("  estado inicial: permanente com q0 = %g mm/h%s\n", op.q0,
                    op.q0Cov > 0. ? (" (normal, CoV " + std::to_string(op.q0Cov) + ")").c_str() : "");
    const std::string chave = ChaveTransiente(op, r);
    std::unique_ptr<Pressoes> pFixas;
    if (!op.pressoes.empty()) {
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        pFixas = std::make_unique<Pressoes>(LePressoes(op.pressoes, chave, mT.GeoMesh(), op.tempos));
        std::printf("  pressoes do transiente lidas de %s (%zu instantes)\n", op.pressoes.c_str(), pFixas->pw.size());
    }
    if (op.soTransiente) { // so o transiente (as medias: k = 1), gravado para as outras rodadas
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        InsertBC(mT);
        mT.SetMaterial(MaterialChuva(r, op));
        mT.CreateCompMesh({{2, op.uT}, {1, 1}, {1, 1}});
        const auto t0 = Clock::now();
        const Pressoes P = CalculaTransiente(mT, r, op, "so_transiente (as medias)");
        if (P.parou >= 0.) {
            std::fprintf(stderr, "o transiente parou em %.2f h: nada gravado\n", P.parou);
            return 2;
        }
        GravaPressoes(op.gravaPressoes, chave, mT.GeoMesh(), P);
        std::printf("  transiente ate %g h em %.1f s (%lld equacoes): %s\n", op.tempos.back(), Seg(t0), (long long)P.eq,
                    op.gravaPressoes.c_str());
        return 0;
    }

    const std::string faixa = std::to_string(op.primeira) + "_" + std::to_string(op.primeira + op.n - 1);
    // amostragem por importancia: a direcao da banda (projecoes) e o deslocamento pelo piloto ou dado
    std::vector<double> dir;
    if (!op.banda.empty()) {
        dir = campos.Direction(ai::LeBanda(op.banda, op.bandaT));
        if (dir.empty()) {
            std::fprintf(stderr, "banda sem pontos: %s (instante %g)\n", op.banda.c_str(), op.bandaT);
            return 1;
        }
    }
    if (!op.ai.piloto.empty() && !ai::DeslocaPeloPiloto(campos, dir, op.ai, op.saida + "_desloc_" + faixa + ".txt"))
        return 3;
    if (!op.deslocXi.empty()) { // deslocamento dado termo a termo (ex.: o ponto de projeto do FORM)
        const std::vector<std::vector<double>> mu = LeXi(op.deslocXi, campos.M(), "desloc_xi");
        campos.SetShift(mu);
        double mu2 = 0.;
        for (const auto &m : mu)
            for (double x : m) mu2 += x * x;
        std::ofstream df(op.saida + "_desloc_" + faixa + ".txt");
        df << std::setprecision(10) << "# deslocamento da amostragem por importancia dado termo a termo\nai=dado\narquivo="
           << op.deslocXi << "\nmodos=" << mu[0].size() << "\nmu2=" << mu2 << "\ntermo,c,phi,k\n";
        for (size_t k = 0; k < mu[0].size(); k++) df << k << "," << mu[0][k] << "," << mu[1][k] << "," << mu[2][k] << "\n";
        std::printf("  amostragem por importancia: deslocamento dado em %s (%zu termos por campo, |mu|^2 = %.4f)\n",
                    op.deslocXi.c_str(), mu[0].size(), mu2);
    }
    std::vector<std::vector<double>> xiDados; // opcao xi=: os K de c', de phi' e de k
    if (!op.xi.empty()) xiDados = LeXi(op.xi, campos.M(), "xi");
    const std::string arq = op.saida + "_" + faixa + ".csv";
    std::ofstream csv(arq);
    csv << kCabecalho;
    // so o bloco que comeca na realizacao 0: no cluster, varios processos rodam blocos do mesmo caso ao mesmo tempo
    if (op.primeira == 0) GravaParametros(op.saida + "_parametros.txt", op, campos, r);
    Mapa mapa;
    std::map<REAL, std::vector<REAL>> fsInstante; // so as realizacoes aleatorias
    std::map<REAL, double> pesoFalha;             // AI: soma de peso * (FS < 1) por instante
    double pesoFalhaMin = 0.;                     // e com o menor FS dos instantes
    std::vector<REAL> fsMinimos;                  // o menor FS dos instantes de cada realizacao com todos ok
    int nPesadas = 0;                             // realizacoes com todos os instantes ok
    double somaT = 0., somaF = 0.;

    std::unique_ptr<std::ofstream> csvChuva, csvCrpc; // caso C: a chuva de cada realizacao e a CRPC
    if (op.chuvaAleatoria || op.q0Cov > 0.) {
        csvChuva = std::make_unique<std::ofstream>(op.saida + "_chuva_" + faixa + ".csv");
        *csvChuva << "real,I_mmh,D_h,q0_mmh,u,v,z_q0\n";
    }
    if (!op.crpc.empty()) {
        csvCrpc = std::make_unique<std::ofstream>(op.saida + "_crpc_" + faixa + ".csv");
        *csvCrpc << "real,q0_mmh,FS0,I,DC_h,falha,t_lo,FS_lo,t_hi,FS_hi,n_srm,s_transiente,s_fator,passos_permanente,"
                    "mudanca_permanente\n";
    }
    std::vector<int> realizacoes;
    if (op.det) realizacoes.push_back(-1);
    for (int n = op.primeira; n < op.primeira + op.n; n++) realizacoes.push_back(n);

    for (int n : realizacoes) {
        const auto t0 = Clock::now();
        const bool deterministica = n < 0;
        if (!deterministica) {
            if (xiDados.empty()) campos.Sample(cp.semente, n);
            else campos.SetXi(xiDados);
        }
        // o peso da amostragem por importancia (0 sem deslocamento) e as projecoes na direcao da banda
        const double logPeso = deterministica ? 0. : campos.LogWeight();
        const std::vector<double> proj = deterministica ? std::vector<double>(3, 0.) : campos.Projection(dir);
        // propriedades de cada ponto: o campo da realizacao; na deterministica, as medias ou os valores dados
        const Campo *cf = deterministica ? (op.valores.empty() ? nullptr : &campoFixo) : &campoAleatorio;
        // caso C: a chuva (I, D) e o q0 da realizacao, numa copia das opcoes e do RainData (daqui ate o fim da realizacao,
        // op e r sao as copias)
        Opcoes opReal = op;
        RainData rReal = r;
        if (!deterministica && (op.chuvaAleatoria || op.q0Cov > 0.)) {
            const ChuvaSorteada c = SorteiaChuva(op, n);
            if (op.q0Cov > 0.) opReal.q0 = c.q0;
            if (op.chuvaAleatoria) {
                opReal.I = rReal.intensity = c.I;
                opReal.chuva = rReal.duration = c.D;
                opReal.tempos = opReal.adapta = {0., c.D};
            }
            *csvChuva << n << "," << std::setprecision(10) << (op.chuvaAleatoria ? c.I : op.I) << ","
                      << (op.chuvaAleatoria ? c.D : op.chuva) << "," << opReal.q0 << "," << c.u << "," << c.v << "," << c.z
                      << "\n";
            csvChuva->flush();
        }
        const Opcoes &op = opReal; // NOLINT: a copia da realizacao esconde as opcoes gerais
        const RainData &r = rReal;
        if (!op.crpc.empty()) { // a curva critica no lugar do Monte Carlo
            CurvaCritica(n, op, r, cf, *csvCrpc);
            std::printf("  realizacao %d: CRPC em %.1f s\n", n, Seg(t0));
            std::fflush(stdout);
            continue;
        }

        // 1. transiente ate o ultimo instante, guardando p_w nos instantes do fator (com k fixo: o mesmo em todas as
        // realizacoes, calculado na primeira ou lido do arquivo)
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        InsertBC(mT);
        mT.SetMaterial(MaterialChuva(r, op));
        mT.CreateCompMesh({{2, op.uT}, {1, 1}, {1, 1}});
        Pressoes calc;
        const Pressoes *P = pFixas.get();
        if (!P) {
            if (cf) AplicaCampos(mT, *cf);
            calc = CalculaTransiente(mT, r, op, "realizacao " + std::to_string(n));
            P = &calc;
            if (calc.parou >= 0.) std::printf("  realizacao %d: o transiente parou em %.2f h\n", n, calc.parou);
            if (!op.gravaPressoes.empty() && calc.parou < 0.) {
                GravaPressoes(op.gravaPressoes, chave, mT.GeoMesh(), calc);
                std::printf("  pressoes do transiente gravadas em %s\n", op.gravaPressoes.c_str());
            }
            if (kFixo && calc.parou < 0.) pFixas = std::make_unique<Pressoes>(calc); // as proximas realizacoes reaproveitam
        }
        const PointLocator loc(mT.GeoMesh());

        // campos completos para figuras: o formato do animacao.py, so nos instantes do fator, e c', phi', k nos nos
        std::unique_ptr<AnimationFiles<TUnsatHP>> anim;
        const bool figura = op.figuras.count(n) > 0;
        if (figura) {
            const std::string pasta = op.saida + "_real" + std::to_string(n);
            anim = std::make_unique<AnimationFiles<TUnsatHP>>(pasta, mT);
            anim->RainSummary(op.I, op.chuva);
            std::ofstream cfs(pasta + "/campos.csv");
            cfs << "no,x,y,c,phi,k\n";
            TPZGeoMesh *g = mT.GeoMesh();
            for (int64_t i = 0; i < g->NNodes(); i++) {
                const REAL x = g->NodeVec()[i].Coord(0), y = g->NodeVec()[i].Coord(1);
                const std::vector<double> v = cf ? (*cf)(x, y) : std::vector<double>{res.c, res.phiDeg, 1.};
                cfs << i << "," << x << "," << y << "," << v[0] << "," << v[1] << "," << v[2] << "\n";
            }
        }

        std::map<REAL, Linha> linhas;
        const double sT = P == &calc ? (calc.s.empty() ? 0. : calc.s.rbegin()->second) : 0.;
        for (REAL t : op.tempos) {
            Linha &l = linhas[t];
            l.real = n, l.t = t;
            if (!P->pw.count(t)) { // o transiente parou antes: o instante fica no arquivo, sem fator
                l.situacao = "transiente", l.sT = sT;
                continue;
            }
            const std::vector<REAL> &pw = P->pw.at(t);
            l.sT = P == &calc ? P->s.at(t) : 0., l.eqT = P->eq;
            l.pwMin = P->faixa.at(t).first, l.pwMax = P->faixa.at(t).second, l.face = P->face.at(t);
            l.nivel30 = NivelAgua(loc, pw, p[4][0]), l.nivel35 = NivelAgua(loc, pw, 0.5 * (p[4][0] + p[3][0]));
            l.nivel40 = NivelAgua(loc, pw, p[3][0]);
            if (anim) anim->Snapshot(t, pw, mT.Material()->Hydraulics());
        }

        // 2. fator: primeiro os instantes adaptativos (nas realizacoes das figuras, todos), depois os outros na malha do
        // adaptativo mais proximo
        std::map<REAL, Fator> adaptados;
        std::map<REAL, REAL> malhaDe; // instante -> o instante adaptativo da malha em que foi calculado
        std::unique_ptr<std::ofstream> mecanismo; // figuras: os pontos do mecanismo em cada instante
        if (figura) {
            mecanismo = std::make_unique<std::ofstream>(op.saida + "_real" + std::to_string(n) + "/mecanismo.csv");
            *mecanismo << "t_h,x,y,epsp_raizJ2,du,c,phi_deg,k_mult\n";
        }
        auto registra = [&](REAL t, REAL fs, Fator &f, const std::vector<REAL> &niveis, bool adapt, double sF) {
            Linha &l = linhas[t];
            l.fs = fs, l.adaptativo = adapt, l.sF = sF, l.eqF = f.m->CompMesh()->NEquations();
            if (!(fs > 0.)) l.situacao = "fator";
            for (size_t i = 0; i < niveis.size(); i++) l.niveis += (i ? ";" : "") + std::to_string(niveis[i]).substr(0, 7);
            l.mec = DescreveMecanismo(f, loc, P->pw.at(t), res);
            if (anim) anim->FactorOfSafety(t, fs, *f.m, *f.a);
            if (mecanismo)
                for (const auto &q : l.mec.pontos)
                    *mecanismo << t << "," << q.x << "," << q.y << "," << q.j << "," << q.du << "," << q.c << ","
                               << q.phi << "," << q.k << "\n";
        };
        double sF = 0.;
        for (REAL t : op.tempos) {
            if (!P->pw.count(t)) continue; // o transiente parou antes
            const bool adapt = figura || std::find(op.adapta.begin(), op.adapta.end(), t) != op.adapta.end();
            if (!adapt) continue;
            const auto tf = Clock::now();
            std::vector<REAL> niveis;
            const REAL fs = FatorAdaptativo(adaptados[t], r, op, cf, loc, P->pw.at(t), niveis);
            adaptados[t].t = t, malhaDe[t] = t;
            const double s = Seg(tf);
            sF += s;
            registra(t, fs, adaptados[t], niveis, true, s);
        }
        for (REAL t : op.tempos) {
            if (!P->pw.count(t) || adaptados.count(t) || adaptados.empty()) continue;
            REAL perto = adaptados.begin()->first; // o instante adaptativo mais proximo
            for (const auto &[ta, f] : adaptados)
                if (std::abs(ta - t) < std::abs(perto - t)) perto = ta;
            const auto tf = Clock::now();
            const REAL fs = CalculaFator(adaptados[perto], cf, loc, P->pw.at(t), op);
            adaptados[perto].t = t, malhaDe[t] = perto;
            const double s = Seg(tf);
            sF += s;
            registra(t, fs, adaptados[perto], {fs}, false, s);
        }

        // 3. arquivo, mapa e resumo da realizacao
        bool todosOk = true;
        REAL fsMin = 1e30, tMin = -1.;
        for (const auto &[t, l] : linhas) {
            Grava(csv, l, op, campos.VarianceFraction(), logPeso, proj);
            if (!deterministica && l.situacao == "ok") {
                fsInstante[t].push_back(l.fs);
                mapa.Soma(t, l.mec, l.fs < 1.);
            }
            todosOk = todosOk && l.situacao == "ok";
            if (l.situacao == "ok" && l.fs < fsMin) fsMin = l.fs, tMin = t;
        }
        if (!deterministica && todosOk) {
            nPesadas++;
            for (const auto &[t, l] : linhas) pesoFalha[t] += (l.fs < 1.) * std::exp(logPeso);
            pesoFalhaMin += (fsMin < 1.) * std::exp(logPeso);
            fsMinimos.push_back(fsMin);
        }
        // o gradiente modal no instante do menor fator (o estado do colapso dele: refaz o fator se a malha ja e de outro)
        if (op.gradiente > 0 && !deterministica && tMin >= 0.) {
            const auto tg = Clock::now();
            Fator &f = adaptados[malhaDe[tMin]];
            if (f.t != tMin) CalculaFator(f, cf, loc, P->pw.at(tMin), op), f.t = tMin;
            const std::vector<double> gr = GradienteModal(f, campos, std::min(op.gradiente, campos.M()), fsMin, op.gim);
            std::ofstream gf(op.saida + "_grad_" + std::to_string(n) + ".txt");
            gf << std::setprecision(10) << "# "
               << (op.gim ? "Gamma (GIM) e o gradiente dGamma/dxi" : "FS (SRM) e o gradiente dFS/dxi")
               << " no instante do menor fator, t = " << tMin
               << " h (os K de c', de phi' e de k), projecao modal\n" << fsMin << "\n";
            for (double x : gr) gf << x << "\n";
            std::printf("  realizacao %d: gradiente modal em t = %g h (%d termos por campo) em %.1f s\n", n, tMin,
                        std::min(op.gradiente, campos.M()), Seg(tg));
        }
        if (!deterministica) mapa.nReal++, somaT += sT, somaF += sF;
        std::printf("  realizacao %4d: FS", n);
        for (const auto &[t, l] : linhas) std::printf(" %g h: %.4f%s", t, l.fs, l.adaptativo ? "" : "*");
        std::printf(" | %.1f s (transiente %.1f s)\n", Seg(t0), sT);
        std::fflush(stdout);
    }
    mapa.Grava(op.saida + "_mapa_" + faixa + ".csv");

    // ---- resumo das realizacoes aleatorias
    std::printf("\n  t (h)   FS medio   desvio   menor   P(FS < 1)   P_f lognormal\n");
    for (const auto &[t, g] : fsInstante) {
        REAL m = 0., v = 0., ml = 0., vl = 0.;
        for (REAL x : g) m += x, ml += std::log(x);
        m /= g.size(), ml /= g.size();
        for (REAL x : g) v += (x - m) * (x - m), vl += (std::log(x) - ml) * (std::log(x) - ml);
        const REAL dp = g.size() > 1 ? std::sqrt(v / (g.size() - 1)) : 0.;
        const REAL dl = g.size() > 1 ? std::sqrt(vl / (g.size() - 1)) : 0.;
        const int falhas = std::count_if(g.begin(), g.end(), [](REAL x) { return x < 1.; });
        const REAL pfLog = dl > 0. ? 0.5 * std::erfc(ml / dl / std::sqrt(2.)) : REAL(ml < 0.);
        std::printf("  %5g   %8.4f   %6.4f   %6.4f   %6.2f%%     %6.2f%%\n", t, m, dp,
                    *std::min_element(g.begin(), g.end()), 100. * falhas / g.size(), 100. * pfLog);
    }
    if (fsMinimos.size() > 1 && op.tempos.size() > 1) { // primeira passagem: o menor FS da historia
        const int falhas = std::count_if(fsMinimos.begin(), fsMinimos.end(), [](REAL x) { return x < 1.; });
        std::printf("  menor FS dos instantes (%zu realizacoes): P(FS_min < 1) = %.2f%%\n", fsMinimos.size(),
                    100. * falhas / fsMinimos.size());
    }
    if ((!op.ai.piloto.empty() || !op.deslocXi.empty()) && nPesadas) { // a P_f da amostragem por importancia, so deste bloco
        std::printf("\n  amostragem por importancia (%d realizacoes com todos os instantes): P_f ponderada", nPesadas);
        for (const auto &[t, w] : pesoFalha) std::printf(" %g h: %.3f%%", t, 100. * w / nPesadas);
        std::printf("; no menor FS: %.3f%%\n", 100. * pesoFalhaMin / nPesadas);
    }
    if (mapa.nReal)
        std::printf("\n  tempo medio por realizacao: transiente %.1f s, fator %.1f s\n", somaT / mapa.nReal,
                    somaF / mapa.nReal);
    std::printf("  arquivos: %s e %s\n", arq.c_str(), (op.saida + "_mapa_" + faixa + ".csv").c_str());
    return 0;
}
