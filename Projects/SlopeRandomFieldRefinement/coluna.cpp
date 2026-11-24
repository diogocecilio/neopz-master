// Coluna 1D de talude infinito sob chuva (Cho 2014, Eng. Geol. 171: 11-20): Richards na vertical numa faixa estreita do
// material nao saturado (RainTransient: esqueleto rigido, ar passivo, chuva no topo com a face de percolacao), um campo
// 1D lognormal de k_s (base de Nystrom do segmento, Projects/Cho2010/nystrom_base.py coluna=H) e o fator de seguranca de
// talude infinito de Lu e Godt nos nos do perfil, calculado depois de cada passo. Em cada realizacao:
//  1. a faixa [0, largura] x [0, H] (MalhaColuna): cada celula com as duas diagonais e um no no centro (tri) ou Q1
//     (quad). As duas sao simetricas no espelho x -> largura - x, e a base, as laterais e o campo tambem: p_w sai
//     exatamente 1D (a tela mostra max |p_L - p_R|). Base e laterais sem fluxo; u fixo; a chuva no topo;
//  2. o k de cada ponto de integracao vem do campo (fKMult); o estado inicial e a succao uniforme (succao0) ou o nivel
//     d'agua (nivel0); o transiente anda nos fins de passo {k dt}, os tempos e o fim da chuva, e o FS e calculado depois
//     de cada passo. A primeira travessia de FS = 1 e localizada no tempo por bissecao (localiza=1: o estado do inicio
//     do passo volta, 7 metades, dt/128, interpolacao linear no ultimo intervalo; depois o estado volta ao fim do passo,
//     na grade de dt). Se o transiente para (o Newton nao converge num sub-passo), o ultimo sub-passo aceito, em
//     t_parou, e um estado convergido: o FS dele entra no FS_passagem e a travessia em [ultimo fim de passo, t_parou]
//     e localizada do mesmo jeito (t_passagem <= t_parou). As linhas depois ficam com situacao=transiente (sem os
//     valores do instante), mas com esse t_passagem e esse FS_passagem; com para_na_falha=1 e travessia nesse
//     intervalo, ficam com falha e o estado de t_parou;
//  3. uma linha do CSV por instante de tempos= (e, se pedidos, o perfil, os xi e o campo).
// O FS (Cho 2014, Eqs. 5, 7 e 9; Lu e Godt 2008), nos nos do perfil (os cantos da lateral esquerda e os centros das
// celulas; no quad, so os cantos), de z = H - y:
//   FS(z) = [c' + (W cos^2 beta - sigma_s) tan phi'] / (W sin beta cos beta),  z > max(0, fs_zmin),
//   W = trapezios de gamma_d + theta gw desde z = 0 (theta = theta_s S_w), sigma_s = S_e p_w se p_w < 0, senao p_w.
// A referencia independente e a Python do referencias_1d/coluna_ref.py (det, mc, xi e fs <perfil.csv>).
// Uso: SlopeRainColumn [chave=valor ...], com (padrao entre parenteses: o Cho 2014):
// A corrida:
//   n (10), primeira (0)   realizacoes primeira..primeira+n-1 (o campo da realizacao n e o RandomField::Sample(semente, n)
//                          com um campo so: o mesmo do coluna_ref.py gerador=pz)
//   det (0)                1: antes, a deterministica (k_s = kmult ks em toda a coluna), como a realizacao -1
//   semente (2025)
//   saida (coluna)         prefixo dos arquivos
//   tempos (6,12,18,24)    instantes das linhas (h); o transiente vai ate o ultimo; 0 e o estado inicial
//   perfis (-1)            grava <saida>_perfil_<real>.csv das realizacoes com real < perfis (o da deterministica sempre):
//                          t_h,z,y,pw,Se,theta,W,FS nos pontos do perfil em cada instante (17 algarismos; o
//                          'coluna_ref.py fs' recalcula o FS dele)
//   grava_xi (0)           1: os M xi de cada realizacao em <saida>_xi_<real>.txt (um por linha; o 'coluna_ref.py xi' le)
//   so_campo (0)           1: so os campos, sem transiente: ln k_mult nos centros das celulas em <saida>_campo_<faixa>.csv
//                          (uma linha por realizacao) e, na tela, as medias, os desvios e a correlacao na distancia lv
//   kmult (1)              k_s / ks da deterministica
//   localiza (1)           0: t_passagem pela interpolacao linear do FS no passo da travessia, sem bissecao
//   para_na_falha (0)      1: o transiente para no passo da travessia (com t_passagem localizado); os instantes depois
//                          ficam com situacao=falha e os valores do fim desse passo (FS_min < 1). Com FS(0) < 1
//                          (t_passagem = 0) nenhum passo e dado: os instantes t > 0 ficam com falha e o estado inicial
// A malha:
//   H (2)                  altura da coluna (m)
//   ny (80)                celulas na vertical. Os V2.5 e V2.2 (~/projects/runs/coluna_estudos): o ny 40 passa no V2.5
//                          (as P_f de 18 e 24 h e o FS de 12 h do MC) mas erra o t_passagem em ate 0,5 h no V2.2
//   largura (H/ny)         largura da faixa (m)
//   elementos (quad)       quad: Q1; tri: cada celula com as duas diagonais (4 triangulos em volta do no do centro).
//                          O quad e o padrao de producao (decisao do usuario no V2.5): com ny 80 passa nos V2.5 e V2.2
//                          como o tri e custa 1/3 (7 s contra 25 s por realizacao)
//   massa (consistente     concentrada: o armazenamento nos vertices (A/3 nos cantos internos e nos centros, A/6 nos
//          no quad,        cantos do topo e da base; so triangulos); consistente: nos pontos de integracao (a unica
//          concentrada     opcao do quad)
//          no tri)
// O solo (van Genuchten-Mualem):
//   vg (0.358,3.58e-4,0.149,1.289)   theta_s,theta_r,alfa,n (alfa em 1/kPa; porosidade = theta_s; S_wr = theta_r/theta_s)
//   ks (3e-6)              k_s medio (m/s)
//   vg_tab (1e-3,reta)     eps[,fecho]: o k_r de Mualem tabelado contra S_w (MualemTabelado.h) ate o teto S_e = 1 - eps
//                          (1e-5 <= eps <= 0,1), fechado ate k_r(1) = 1 por reta, c1 ou ippisch; 0: o EMualem (com n < 2
//                          o Newton cicla perto da saturacao: so como diagnostico). O eps e parte do modelo: no Cho, a P_f
//                          de 24 h cai 1,5 ponto de eps 1e-3 para 1e-4 (V2.7, ~/projects/runs/coluna_estudos)
//   vg_tab_psi ()          psi_c[,fecho]: o mesmo com eps = 1 - S_e(psi_c) (kPa); ex.: vg_tab_psi=0.01,ippisch. Um ou outro
//   Kw (2.2e6), gw (9.81)  modulo volumetrico da agua (kPa) e peso especifico da agua (kN/m3)
// A chuva e o estado inicial:
//   I (10.8)               intensidade (mm/h), vertical, no topo (onde a agua empocaria, a face de percolacao p_w = 0)
//   chuva (24)             duracao (h); depois, fluxo nulo no topo
//   succao0 (20)           succao inicial uniforme (kPa)
//   nivel0 ()              cota y do nivel d'agua inicial (m): p_w = gw (nivel0 - y), com a succao ate succao_max (no
//                          lugar de succao0; nivel0=H: coluna cheia e hidrostatica)
//   succao_max (1e9)       kPa, so com nivel0
// A numerica:
//   dt (0.1), dt_min (0.00625)   passo (h) e o menor da subdivisao de um passo que nao converge (h; 0: para)
//   newton (1e-6,80,1)     tol_p,maxit[,busca]: |dp| < tol_p (kPa), maxit e as metades da busca linear (se nao vier: 1).
//                          O padrao vem do V2.6 (~/projects/runs/coluna_estudos): com o do TPZUnsatAnalysis (1e-6,40 e
//                          busca 4) o transiente desistia do passo em 5 a 20% das realizacoes (caso base e lv 16, tri
//                          ny 40 e 80), quando a frente chega a base (o Newton amortecido anda devagar na zona saturada);
//                          com busca 0, em 3% (CoV 1,5: a agua empocada sobre uma camada pouco permeavel faz o Newton
//                          oscilar). newton=1e-6,40,4 refaz o padrao antigo
//   solver (skyline)       skyline | pardiso: o LU do TPZUnsatAnalysis (SetUseSkyline), montagem serial. O skyline e o
//                          padrao: a matriz da faixa e de banda, o custo e o mesmo (V2.3) e nao depende das threads do MKL
//   newton_trace (<arquivo>)  diagnostico de um transiente que para (RainTransient::SetNewtonTrace), um bloco por parada
//   face_tol (0)           um no da face so volta ao fluxo se R / (dt L) > q (1 + face_tol) + 1e-10
// O FS:
//   fs (talude)            talude | nenhum (sem FS: FS_min = nan e sem travessia)
//   beta (40), c (5), phi (32), gamma_d (16)   graus, kPa, graus, kN/m3
//   fs_peso (integrado)    integrado: W = trapezios de gamma_d + theta gw; constante: W = gamma z
//   gamma (20)             kN/m3, so com fs_peso=constante
//   fs_zmin (0)            o FS so conta em z > fs_zmin (m)
// O campo (ln k_s gaussiano de media -sigma^2/2: k_mult de media 1):
//   covk (0.5)             CoV de k_s
//   lv (0.5)               distancia de autocorrelacao vertical (m), exp(-|dz|/lv); tem de ser a da base
//   M (100)                termos da KL (os M primeiros da base)
//   nystrom (<arquivo>)    a base 1D (python3 -I Projects/Cho2010/nystrom_base.py base.txt coluna=H lv=.. M=.. [dy=..]):
//                          obrigatoria com realizacoes aleatorias e covk > 0
//   xi (<arquivo>)         os xi dados (ate M numeros, um campo), sem sorteio, em todas as realizacoes
//   desloc_xi (<arquivo>)  amostragem por importancia: o deslocamento dos xi termo a termo (ate M numeros); log_peso
// Reservadas (spec B.6, depois do lancamento; hoje dao erro): chuva_tabela, chuva_cos, inclinacao, gardner, q0, base,
// carga_base, sensores, sensores_dt
// Saidas (<faixa> = primeira_ultima; so a deterministica: <saida>_0_-1.csv):
//   <saida>_<faixa>.csv    uma linha por realizacao e instante: real,t_h,FS_min,situacao,z_crit,base,FS_passagem,
//                          t_passagem,frente_m,pw_base,pw_topo,nos_face,passos,newton,s_transiente,t_parou,k_mult_med,
//                          k_mult_min,k_mult_max,log_peso,ny,dt,elementos,massa,lei_k,cov_k,lv,M,KL_variancia
//                          situacao: ok, transiente (o transiente parou antes; t_parou) ou falha (para_na_falha=1).
//                          FS_passagem: o menor FS ate o instante (todos os passos, e o ultimo sub-passo aceito de um
//                          transiente que parou); t_passagem: a primeira travessia de FS = 1 ate o instante (-1:
//                          nenhuma; numa linha transiente, ate t_parou); frente_m: o primeiro no, de cima para baixo,
//                          com p_w - p_w0 < 1 kPa (o pe da frente: no Cho, 4 a 5 cm abaixo do meio do salto, p_w =
//                          -10 kPa; o V3.1 confere o pe da frente); base: z_crit >= H - 1e-9; k_mult_*: nos pontos de integracao. Na
//                          deterministica cov_k, lv, M e KL_variancia sao 0
//   <saida>_parametros.txt chave=valor do caso (so o bloco que comeca em 0, ou a deterministica sozinha)
//   <saida>_perfil_<real>.csv, <saida>_xi_<real>.txt, <saida>_campo_<faixa>.csv   (perfis=, grava_xi=, so_campo=)
// Na tela: por realizacao, o tempo, os passos, as iteracoes de Newton, max |p_L - p_R| e |dW - infiltrado| (a agua
// armazenada contra a que entrou pelo topo); no fim, o tempo medio e a linha 'arquivos:'. Codigo de saida 2: a
// deterministica sozinha (det=1 n=0) parou antes do ultimo instante.
#include "HPMaterial.h"
#include "MualemTabelado.h"
#include "RainTransient.h"
#include "RandomField.h"
#include "UnsatModel.h"

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzquad.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
static double Seg(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }
static const REAL kNaN = std::numeric_limits<REAL>::quiet_NaN();
static const REAL kInf = std::numeric_limits<REAL>::infinity();

// ids da faixa (os do montecarlo.cpp onde existem: o topo e o id da face)
enum EIds { matSolo = 1, base = -1, direita = -2, esquerda = -5, topo = -6 };

// ---- opcoes (o cabecalho explica cada uma)
struct Opcoes {
    // corrida
    int n = 10, primeira = 0, perfis = -1;
    bool det = false, gravaXi = false, soCampo = false, localiza = true, paraNaFalha = false;
    uint64_t semente = 2025;
    std::string saida = "coluna";
    std::vector<REAL> tempos = {6., 12., 18., 24.};
    REAL kmult = 1.;
    // malha
    REAL H = 2., largura = 0.;
    int ny = 80;
    bool quad = true, concentrada = false; // o quad (V2.5)
    // solo
    REAL ts = 0.358, tr = 3.58e-4, a = 0.149, nvg = 1.289, ks = 3.e-6, Kw = 2.2e6, gw = 9.81;
    REAL vgTab = 1.e-3, vgTabPsi = -1.; // eps do teto (0: EMualem); psi_c do vg_tab_psi (< 0: nao dado)
    EFechoVG fecho = EFechoVG::Reta;
    std::string leiK = "tab_0.001_reta";
    // chuva e estado inicial
    REAL I = 10.8, chuva = 24., succao0 = 20., nivel0 = kNaN, succaoMax = 1.e9;
    // numerica
    REAL dt = 0.1, dtMin = 0.00625, faceTol = 0.;
    std::vector<REAL> newton = {1.e-6, 80., 1.}; // tol_p, maxit, busca (o V2.6)
    bool skyline = true; // o LU skyline (o V2.3 da o mesmo tempo do PARDISO)
    std::string newtonTrace;
    // FS
    bool fs = true, pesoIntegrado = true;
    REAL beta = 40., c = 5., phi = 32., gammaD = 16., gamma = 20., fsZmin = 0.;
    // campo
    REAL covk = 0.5, lv = 0.5;
    int M = 100;
    std::string nystrom, xi, deslocXi;
};

// numeros das opcoes: a cadeia inteira tem de ser o numero (std::stod/stoi aceitam '10.7' -> 10 e lancam excecao com
// 'abc'); um valor invalido e um std::invalid_argument com o texto, que o LeOpcoes transforma em Erro
REAL Numero(const std::string &s)
{
    size_t pos = 0;
    REAL v = 0.;
    try {
        v = std::stod(s, &pos);
    } catch (const std::exception &) {
        throw std::invalid_argument(s);
    }
    if (pos != s.size()) throw std::invalid_argument(s);
    return v;
}

int Inteiro(const std::string &s)
{
    size_t pos = 0;
    int v = 0;
    try {
        v = std::stoi(s, &pos);
    } catch (const std::exception &) {
        throw std::invalid_argument(s);
    }
    if (pos != s.size()) throw std::invalid_argument(s);
    return v;
}

std::vector<REAL> Lista(const std::string &s)
{
    std::vector<REAL> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(Numero(item));
    return v;
}

[[noreturn]] void Erro(const std::string &m)
{
    std::fprintf(stderr, "%s\n", m.c_str());
    std::exit(1);
}

std::string G(REAL v) // %g, para os rotulos
{
    char b[64];
    std::snprintf(b, sizeof b, "%g", v);
    return b;
}

Opcoes LeOpcoes0(int argc, char *argv[]);

// as opcoes; um numero mal formado e um Erro com a chave (e nao um abort por excecao)
Opcoes LeOpcoes(int argc, char *argv[])
{
    try {
        return LeOpcoes0(argc, argv);
    } catch (const std::invalid_argument &x) {
        Erro(std::string(x.what()) + ": valor invalido");
    }
}

Opcoes LeOpcoes0(int argc, char *argv[])
{
    const std::set<std::string> chaves = {
        "n",         "primeira",  "det",        "semente",    "saida",      "tempos",       "perfis",  "grava_xi",
        "so_campo",  "kmult",     "localiza",   "para_na_falha", "H",       "ny",           "largura", "elementos",
        "massa",     "vg",        "ks",         "vg_tab",     "vg_tab_psi", "Kw",           "gw",      "I",
        "chuva",     "succao0",   "nivel0",     "succao_max", "dt",         "dt_min",       "newton",  "solver",
        "newton_trace", "face_tol", "fs",       "beta",       "c",          "phi",          "gamma_d", "fs_peso",
        "gamma",     "fs_zmin",   "covk",       "lv",         "M",          "nystrom",      "xi",      "desloc_xi"};
    // depois do lancamento (spec B.6): reservadas
    const std::set<std::string> reservadas = {"chuva_tabela", "chuva_cos", "inclinacao", "gardner", "q0",
                                              "base",         "carga_base", "sensores", "sensores_dt"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        const std::string k = eq == std::string::npos ? a : a.substr(0, eq);
        if (eq != std::string::npos && reservadas.count(k))
            Erro(k + ": opcao reservada, ainda nao implementada (depois do lancamento, spec B.6)");
        if (eq == std::string::npos || !chaves.count(k))
            Erro("opcao desconhecida: " + a + " (as opcoes estao no cabecalho do coluna.cpp)");
        kv[k] = a.substr(eq + 1);
    }
    auto tem = [&kv](const char *k) { return kv.count(k) > 0; };
    // um numero mal formado: invalid_argument("<chave>=<valor>"), o Erro do LeOpcoes
    auto valor = [&](const char *k, auto conv) {
        try {
            return conv(kv[k]);
        } catch (const std::invalid_argument &) {
            throw std::invalid_argument(std::string(k) + "=" + kv[k]);
        }
    };
    auto real = [&](const char *k, REAL &v) {
        if (tem(k)) v = valor(k, Numero);
    };
    auto inteiro = [&](const char *k, int &v) {
        if (tem(k)) v = valor(k, Inteiro);
    };
    auto logico = [&](const char *k, bool &v) {
        if (tem(k)) v = valor(k, Inteiro) != 0;
    };
    Opcoes op;
    inteiro("n", op.n);
    inteiro("primeira", op.primeira);
    inteiro("perfis", op.perfis);
    logico("det", op.det);
    logico("grava_xi", op.gravaXi);
    logico("so_campo", op.soCampo);
    logico("localiza", op.localiza);
    logico("para_na_falha", op.paraNaFalha);
    if (tem("semente"))
        op.semente = valor("semente", [](const std::string &x) {
            size_t pos = 0;
            uint64_t v = 0;
            try {
                v = std::stoull(x, &pos);
            } catch (const std::exception &) {
                throw std::invalid_argument(x);
            }
            if (pos != x.size() || x.find('-') != std::string::npos) throw std::invalid_argument(x);
            return v;
        });
    if (tem("saida")) op.saida = kv["saida"];
    if (tem("tempos")) op.tempos = valor("tempos", Lista);
    std::sort(op.tempos.begin(), op.tempos.end());
    op.tempos.erase(std::unique(op.tempos.begin(), op.tempos.end()), op.tempos.end());
    if (op.tempos.empty() || op.tempos.front() < 0.) Erro("tempos: instantes >= 0 (h)");
    real("kmult", op.kmult);
    if (!(op.kmult > 0.)) Erro("kmult > 0");
    if (op.n < 0) Erro("n >= 0");
    if (op.primeira < 0) Erro("primeira >= 0 (a realizacao -1 e o rotulo da deterministica)");
    if (op.n == 0 && !op.det) Erro("nada a fazer: n = 0 e det = 0");
    // malha
    real("H", op.H);
    inteiro("ny", op.ny);
    if (!(op.H > 0.) || op.ny < 1) Erro("H > 0 (m) e ny >= 1");
    op.largura = op.H / op.ny;
    real("largura", op.largura);
    if (!(op.largura > 0.)) Erro("largura > 0 (m)");
    if (tem("elementos")) {
        if (kv["elementos"] != "tri" && kv["elementos"] != "quad") Erro("elementos: tri ou quad");
        op.quad = kv["elementos"] == "quad";
    }
    op.concentrada = !op.quad; // o quad so tem o armazenamento nos pontos
    if (tem("massa")) {
        if (kv["massa"] != "concentrada" && kv["massa"] != "consistente") Erro("massa: concentrada ou consistente");
        op.concentrada = kv["massa"] == "concentrada";
        if (op.quad && op.concentrada)
            Erro("massa=concentrada so com elementos=tri: o material so concentra o armazenamento nos triangulos P1 "
                 "(TPZMatUnsaturatedUPwPa.cpp:240); o quad usa massa=consistente");
    }
    // solo
    if (tem("vg")) {
        const std::vector<REAL> v = valor("vg", Lista);
        if (v.size() != 4 || !(v[0] > v[1] && v[1] >= 0. && v[0] < 1. && v[2] > 0. && v[3] > 1.))
            Erro("vg: theta_s,theta_r,alfa,n com 1 > theta_s > theta_r >= 0, alfa > 0 (1/kPa) e n > 1");
        op.ts = v[0], op.tr = v[1], op.a = v[2], op.nvg = v[3];
    }
    real("ks", op.ks);
    real("Kw", op.Kw);
    real("gw", op.gw);
    if (!(op.ks > 0. && op.Kw > 0. && op.gw > 0.)) Erro("ks, Kw e gw > 0");
    if (tem("vg_tab") && tem("vg_tab_psi")) Erro("vg_tab e vg_tab_psi: um ou outro");
    for (const std::string k : {"vg_tab", "vg_tab_psi"}) { // k_r de Mualem tabelado: eps (ou psi_c) e o fecho
        if (!tem(k.c_str())) continue;
        const size_t vir = kv[k].find(',');
        REAL x = 0.;
        try {
            x = Numero(kv[k].substr(0, vir));
        } catch (const std::invalid_argument &) {
            throw std::invalid_argument(k + "=" + kv[k]);
        }
        if (k == "vg_tab" && x == 0.) { // vg_tab=0: o EMualem
            op.vgTab = 0.;
            break;
        }
        const std::string fecho = vir == std::string::npos ? "reta" : kv[k].substr(vir + 1);
        if (fecho != "reta" && fecho != "c1" && fecho != "ippisch") Erro(k + ": o fecho e reta, c1 ou ippisch");
        op.fecho = fecho == "reta" ? EFechoVG::Reta : fecho == "c1" ? EFechoVG::C1 : EFechoVG::Ippisch;
        if (k == "vg_tab_psi" && !(x > 0.)) Erro("vg_tab_psi: psi_c > 0 (kPa)");
        op.vgTab = k == "vg_tab" ? x : EpsDePsi(op.a, op.nvg, x);
        if (k == "vg_tab_psi") op.vgTabPsi = x;
        if (!(op.vgTab >= 1.e-5 && op.vgTab <= 0.1)) {
            char b[128];
            std::snprintf(b, sizeof b, "%s: eps = %g fora de [1e-5, 0.1]%s", k.c_str(), op.vgTab,
                          k == "vg_tab" ? "" : " (eps = 1 - S_e(psi_c))");
            Erro(b);
        }
    }
    // o rotulo da lei de k (coluna lei_k; os do coluna_ref.py)
    op.leiK = op.vgTab == 0. ? "vg"
              : op.vgTabPsi > 0. ? "tabpsi_" + G(op.vgTabPsi) + "_" + NomeFechoVG(op.fecho)
                                 : "tab_" + G(op.vgTab) + "_" + NomeFechoVG(op.fecho);
    // chuva e estado inicial
    real("I", op.I);
    real("chuva", op.chuva);
    if (op.I < 0. || op.chuva < 0.) Erro("I >= 0 (mm/h) e chuva >= 0 (h)");
    if (tem("succao0") && tem("nivel0")) Erro("succao0 e nivel0: um ou outro");
    real("succao0", op.succao0);
    real("nivel0", op.nivel0);
    real("succao_max", op.succaoMax);
    if (tem("succao_max") && !tem("nivel0")) Erro("succao_max so com nivel0= (com succao0, a succao e uniforme)");
    if (op.succao0 < 0. || !(op.succaoMax > 0.)) Erro("succao0 >= 0 e succao_max > 0 (kPa)");
    // numerica
    real("dt", op.dt);
    real("dt_min", op.dtMin);
    if (!(op.dt > 0.) || op.dtMin < 0.) Erro("dt > 0 e dt_min >= 0 (h)");
    if (tem("newton")) {
        op.newton = valor("newton", Lista);
        if (op.newton.size() < 2 || op.newton.size() > 3 || !(op.newton[0] > 0. && op.newton[1] >= 1.) ||
            (op.newton.size() > 2 && op.newton[2] < 0.) || op.newton[1] != std::floor(op.newton[1]) ||
            (op.newton.size() > 2 && op.newton[2] != std::floor(op.newton[2])))
            Erro("newton: tol_p,maxit[,busca] (kPa, inteiro >= 1, metades da busca linear: inteiro >= 0)");
        if (op.newton.size() == 2) op.newton.push_back(1.); // a busca padrao
    }
    if (tem("solver")) {
        if (kv["solver"] != "pardiso" && kv["solver"] != "skyline") Erro("solver: skyline ou pardiso");
        op.skyline = kv["solver"] == "skyline";
    }
    if (tem("newton_trace")) op.newtonTrace = kv["newton_trace"];
    real("face_tol", op.faceTol);
    if (!(op.faceTol >= 0.)) Erro("face_tol >= 0");
    // FS
    if (tem("fs")) {
        if (kv["fs"] != "talude" && kv["fs"] != "nenhum") Erro("fs: talude ou nenhum");
        op.fs = kv["fs"] == "talude";
    }
    if (tem("fs_peso")) {
        if (kv["fs_peso"] != "integrado" && kv["fs_peso"] != "constante") Erro("fs_peso: integrado ou constante");
        op.pesoIntegrado = kv["fs_peso"] == "integrado";
    }
    real("beta", op.beta);
    real("c", op.c);
    real("phi", op.phi);
    real("gamma_d", op.gammaD);
    real("gamma", op.gamma);
    real("fs_zmin", op.fsZmin);
    if (!(op.beta > 0. && op.beta < 90. && op.phi >= 0. && op.phi < 90. && op.c >= 0. && op.gammaD > 0. && op.gamma > 0.))
        Erro("FS: 0 < beta < 90, 0 <= phi < 90 (graus), c >= 0 (kPa), gamma_d e gamma > 0 (kN/m3)");
    // campo
    real("covk", op.covk);
    real("lv", op.lv);
    inteiro("M", op.M);
    if (op.covk < 0. || !(op.lv > 0.) || op.M < 1) Erro("covk >= 0, lv > 0 (m) e M >= 1");
    if (tem("nystrom")) op.nystrom = kv["nystrom"];
    if (tem("xi")) op.xi = kv["xi"];
    if (tem("desloc_xi")) op.deslocXi = kv["desloc_xi"];
    if (!op.xi.empty() && !op.deslocXi.empty()) Erro("xi e desloc_xi: um ou outro");
    const bool aleatorias = op.n > 0 && (op.covk > 0. || op.soCampo);
    if (aleatorias && op.nystrom.empty())
        Erro("as realizacoes aleatorias precisam da base 1D: nystrom=<base.txt> (python3 -I Projects/Cho2010/"
             "nystrom_base.py base.txt coluna=H lv=.. M=..)");
    if (op.soCampo && op.n < 1) Erro("so_campo=1 precisa de n >= 1");
    return op;
}

// ---- a malha
// A faixa [0, largura] x [0, H] com ny celulas: nos dos cantos 2j (x = 0) e 2j + 1 (x = largura) na cota y_j = H j / ny
// e, no tri, os centros 2 (ny + 1) + j; cada celula em 4 triangulos em volta do centro (base, direita, topo, esquerda,
// todos anti-horarios) ou num Q1. Contornos: base (-1), direita (-2), topo (-6) e esquerda (-5), no sentido anti-horario.
// O LayerGMesh nao serve: o desempate dele (LayerGMesh.cpp:81) alterna as diagonais e a faixa deixa de ser simetrica
TPZGeoMesh *MalhaColuna(REAL H, REAL largura, int ny, bool quad)
{
    auto *g = new TPZGeoMesh();
    g->SetDimension(2);
    const int64_t nCantos = 2 * int64_t(ny + 1);
    g->NodeVec().Resize(nCantos + (quad ? 0 : ny));
    auto poe = [g](int64_t id, REAL x, REAL y) {
        TPZManVector<REAL, 3> c = {x, y, 0.};
        g->NodeVec()[id] = TPZGeoNode(id, c, *g);
    };
    auto Y = [H, ny](int j) { return j == ny ? H : H * j / ny; };
    for (int j = 0; j <= ny; j++) poe(2 * j, 0., Y(j)), poe(2 * j + 1, largura, Y(j));
    if (!quad)
        for (int j = 0; j < ny; j++) poe(nCantos + j, 0.5 * largura, 0.5 * (Y(j) + Y(j + 1)));
    for (int j = 0; j < ny; j++) {
        const int64_t bl = 2 * j, br = 2 * j + 1, tl = 2 * j + 2, tr = 2 * j + 3;
        if (quad) {
            TPZManVector<int64_t, 4> v = {bl, br, tr, tl};
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(v, matSolo, *g);
            continue;
        }
        const int64_t c = nCantos + j;
        for (const auto &t : {std::array<int64_t, 3>{bl, br, c}, std::array<int64_t, 3>{br, tr, c},
                              std::array<int64_t, 3>{tr, tl, c}, std::array<int64_t, 3>{tl, bl, c}}) {
            TPZManVector<int64_t, 3> v = {t[0], t[1], t[2]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(v, matSolo, *g);
        }
    }
    auto linha = [g](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> v = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(v, id, *g);
    };
    linha(0, 1, base);
    for (int j = 0; j < ny; j++) linha(2 * j + 1, 2 * j + 3, direita);
    linha(2 * int64_t(ny) + 1, 2 * int64_t(ny), topo);
    for (int j = ny; j > 0; j--) linha(2 * j, 2 * j - 2, esquerda);
    g->BuildConnectivity();
    return g;
}

// Os pontos do perfil, de cima (z = 0) para baixo (z = H): os cantos da lateral esquerda e, no tri, os centros entre eles;
// os pares (canto esquerdo, canto direito) de cada cota, para a simetria; e c' (kPa) e tan(phi') de cada ponto (hoje os
// das opcoes c= e phi=; o gancho para os artigos com a resistencia variando na profundidade, ex. o Liu 2025 por celula)
struct Perfil {
    std::vector<int64_t> nos;
    std::vector<REAL> z, y, c, tanPhi;
    std::vector<std::pair<int64_t, int64_t>> pares;
};

Perfil FazPerfil(TPZGeoMesh *g, const Opcoes &op)
{
    Perfil P;
    const int64_t nCantos = 2 * int64_t(op.ny + 1);
    const REAL tp = std::tan(op.phi * M_PI / 180.);
    auto poe = [&](int64_t no) {
        const REAL y = g->NodeVec()[no].Coord(1);
        P.nos.push_back(no), P.y.push_back(y), P.z.push_back(op.H - y), P.c.push_back(op.c), P.tanPhi.push_back(tp);
    };
    for (int j = op.ny; j >= 0; j--) {
        poe(2 * j);
        if (!op.quad && j > 0) poe(nCantos + j - 1);
    }
    for (int j = 0; j <= op.ny; j++) P.pares.push_back({2 * j, 2 * j + 1});
    return P;
}

// ---- o material (o MaterialChuva do montecarlo.cpp com van Genuchten-Mualem): esqueleto rigido e ar passivo, entao a
// resistencia, o peso e o chi nao entram no transiente
TabelaVG gTabelaVG; // k_r de Mualem tabelado ate S_e = 1 - eps (MualemTabelado.h), preenchida no main

TPZUnsaturatedHydraulics Hidraulica(const Opcoes &op)
{
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten; // S_w = S_wr + (1 - S_wr) S_e, theta = theta_s S_w
    h.fAVG = op.a, h.fBVG = op.nvg, h.fSwr = op.tr / op.ts;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    if (op.vgTab > 0.) UsaTabela(h, gTabelaVG); // vg_tab=: Mualem tabelado ate S_e = 1 - eps
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // o ar fica passivo (p_a = 0)
    return h;
}

TUnsatHP *MaterialColuna(const Opcoes &op)
{
    const REAL phi = op.phi * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(2.e4, 0.3);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, op.c, ER), ER));
    mat->SetHydraulics(Hidraulica(op));
    mat->SetBishop(TUnsat::EChiEffective); // chi = S_e, o do FS de Lu e Godt (nao entra: u fixo)
    mat->SetPorosity(op.ts);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(op.Kw);
    mat->SetConductivity(op.ks, op.gw, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(op.gamma);
    mat->SetLumpedStorage(op.concentrada);
    mat->SetHydroMechanicalCoupling(false); // esqueleto rigido na percolacao
    return mat;
}

RainData DadosChuva(const Opcoes &op)
{
    RainData r;
    r.ks = op.ks, r.n = op.ts, r.Kw = op.Kw; // so informativos: o material tem os seus
    r.intensity = op.I, r.duration = op.chuva, r.dt = op.dt, r.dtMin = op.dtMin;
    if (op.newton.size() > 2) r.lineSearch = int(op.newton[2]);
    r.rainVertical = false; // o topo e horizontal
    // o nivel do RainData e a reta yToe - Hw + (xToe - x) tan(betaW): aqui horizontal, na cota yToe
    r.xToe = 0., r.Hw = 0., r.betaW = 0.;
    if (std::isnan(op.nivel0)) r.yToe = -1.e6, r.suctionMax = op.succao0; // succao uniforme
    else r.yToe = op.nivel0, r.suctionMax = op.succaoMax;
    return r;
}

// ---- a agua na coluna (o balanco): as mesmas integrais do material
// Em cada elemento do dominio, os pontos da regra do material (ordem 2 do HPMaterial com u P1) com o peso w |det J| e as
// funcoes P1 (ou Q1) dos vertices no ponto: f(gel, peso, N)
template <class F>
void PontosDaRegra(TPZGeoMesh *g, F f)
{
    for (TPZGeoEl *gel : g->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->MaterialId() != matSolo) continue;
        std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, 2));
        TPZManVector<int, 3> ordem(2, 2);
        regra->SetOrder(ordem);
        const bool tri = gel->NCornerNodes() == 3;
        for (int ip = 0; ip < regra->NPoints(); ip++) {
            TPZManVector<REAL, 3> qsi(2), x(3);
            REAL w, detjac;
            regra->Point(ip, qsi, w);
            TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
            std::array<REAL, 4> N;
            if (tri) N = {1. - qsi[0] - qsi[1], qsi[0], qsi[1], 0.};
            else
                N = {0.25 * (1. - qsi[0]) * (1. - qsi[1]), 0.25 * (1. + qsi[0]) * (1. - qsi[1]),
                     0.25 * (1. + qsi[0]) * (1. + qsi[1]), 0.25 * (1. - qsi[0]) * (1. + qsi[1])};
            // as funcoes tem de reproduzir o mapeamento do elemento
            gel->X(qsi, x);
            for (int d = 0; d < 2; d++) {
                REAL v = 0.;
                for (int i = 0; i < gel->NCornerNodes(); i++) v += N[i] * gel->NodePtr(i)->Coord(d);
                if (std::abs(v - x[d]) > 1.e-12 * (1. + std::abs(x[d]))) DebugStop();
            }
            f(gel, w * std::abs(detjac), N);
        }
    }
}

// n S_w integrado como no armazenamento do material: concentrado nos vertices dos triangulos ((int N_j) n S(p_j)) ou nos
// pontos (n S(p) no ponto); m3 por metro de espessura
REAL AguaArmazenada(TPZGeoMesh *g, const std::vector<REAL> &pw, const TPZUnsaturatedHydraulics &h, REAL n, bool concentrada)
{
    REAL soma = 0.;
    PontosDaRegra(g, [&](TPZGeoEl *gel, REAL peso, const std::array<REAL, 4> &N) {
        const int nv = gel->NCornerNodes();
        if (concentrada && nv == 3) {
            for (int i = 0; i < 3; i++) {
                REAL S, Gm;
                h.Retention(-pw[gel->NodeIndex(i)], S, Gm);
                soma += peso * N[i] * n * S;
            }
            return;
        }
        REAL p = 0.;
        for (int i = 0; i < nv; i++) p += N[i] * pw[gel->NodeIndex(i)];
        REAL S, Gm;
        h.Retention(-p, S, Gm);
        soma += peso * n * S;
    });
    return soma;
}

// O termo da compressibilidade da agua num passo, nos pontos (n S(p) / K_w (p - p_ant)): o resto da agua armazenada
REAL Compressibilidade(TPZGeoMesh *g, const std::vector<REAL> &pw, const std::vector<REAL> &pwAnt,
                       const TPZUnsaturatedHydraulics &h, REAL n, REAL Kw)
{
    REAL soma = 0.;
    PontosDaRegra(g, [&](TPZGeoEl *gel, REAL peso, const std::array<REAL, 4> &N) {
        REAL p = 0., pa = 0.;
        for (int i = 0; i < gel->NCornerNodes(); i++) p += N[i] * pw[gel->NodeIndex(i)], pa += N[i] * pwAnt[gel->NodeIndex(i)];
        REAL S, Gm;
        h.Retention(-p, S, Gm);
        soma += peso * n * S / Kw * (p - pa);
    });
    return soma;
}

// ---- o FS de Lu e Godt nos pontos do perfil
struct Estado {
    REAL fs = kNaN, zCrit = kNaN, frente = kNaN, pwBase = kNaN, pwTopo = kNaN, simetria = 0.;
    int base = 0;
    size_t face = 0;
    std::vector<REAL> pw, Se, theta, W, FS; // nos pontos do perfil, de cima para baixo
};

Estado Avalia(const std::vector<REAL> &pwNos, const Perfil &P, const TPZUnsaturatedHydraulics &h, const Opcoes &op,
              const std::vector<REAL> &pw0)
{
    Estado e;
    const size_t np = P.nos.size();
    e.pw.resize(np), e.Se.resize(np), e.theta.resize(np), e.W.assign(np, 0.), e.FS.assign(np, kInf);
    for (size_t i = 0; i < np; i++) {
        const REAL p = pwNos[P.nos[i]];
        REAL S, Gm, Se, dSe;
        h.Retention(-p, S, Gm);
        h.Effective(S, Se, dSe);
        e.pw[i] = p, e.Se[i] = Se, e.theta[i] = op.ts * S;
    }
    for (size_t i = 1; i < np; i++) { // W desde z = 0: trapezios de gamma_t = gamma_d + theta gw, ou gamma z
        const REAL gi = op.gammaD + e.theta[i] * op.gw, ga = op.gammaD + e.theta[i - 1] * op.gw;
        e.W[i] = op.pesoIntegrado ? e.W[i - 1] + 0.5 * (gi + ga) * (P.z[i] - P.z[i - 1]) : op.gamma * P.z[i];
    }
    if (op.fs) {
        const REAL b = op.beta * M_PI / 180.;
        const REAL c2 = std::cos(b) * std::cos(b), sc = std::sin(b) * std::cos(b), zmin = std::max(REAL(0.), op.fsZmin);
        size_t k = np;
        for (size_t i = 0; i < np; i++) {
            if (!(P.z[i] > zmin + 1.e-12)) continue;
            const REAL sig = e.pw[i] < 0. ? e.Se[i] * e.pw[i] : e.pw[i];
            e.FS[i] = (P.c[i] + (e.W[i] * c2 - sig) * P.tanPhi[i]) / (e.W[i] * sc);
            if (k == np || e.FS[i] < e.FS[k]) k = i;
        }
        if (k < np) e.fs = e.FS[k], e.zCrit = P.z[k], e.base = P.z[k] >= op.H - 1.e-9;
    }
    // a frente: o primeiro ponto, de cima para baixo, que ainda nao molhou (p_w - p_w0 < 1 kPa)
    e.frente = P.z.back();
    for (size_t i = 0; i < np; i++)
        if (e.pw[i] - pw0[i] < 1.) {
            e.frente = P.z[i];
            break;
        }
    e.pwBase = e.pw.back(), e.pwTopo = e.pw.front();
    for (const auto &[l, r] : P.pares) e.simetria = std::max(e.simetria, std::abs(pwNos[l] - pwNos[r]));
    return e;
}

// ---- saidas
const char *kCabecalho =
    "real,t_h,FS_min,situacao,z_crit,base,FS_passagem,t_passagem,frente_m,pw_base,pw_topo,nos_face,passos,newton,"
    "s_transiente,t_parou,k_mult_med,k_mult_min,k_mult_max,log_peso,ny,dt,elementos,massa,lei_k,cov_k,lv,M,KL_variancia\n";

std::string F(REAL v) // %.10g (nan e inf por extenso)
{
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0. ? "inf" : "-inf";
    char b[64];
    std::snprintf(b, sizeof b, "%.10g", v);
    return b;
}

// O que vai em todas as linhas de uma realizacao
struct Fixos {
    int real = -1;
    REAL kMed = 1., kMin = 1., kMax = 1., logPeso = 0., covk = 0., lv = 0., klVar = 0.;
    int M = 0;
};

struct Linha {
    REAL t = 0., fsPass = kNaN, tPass = -1., tParou = -1.;
    std::string situacao = "ok";
    Estado e;
    int64_t passos = 0, newton = 0;
    double s = 0.;
};

void Grava(std::ostream &csv, const Linha &l, const Fixos &f, const Opcoes &op)
{
    const Estado &e = l.e;
    csv << f.real << "," << F(l.t) << "," << F(e.fs) << "," << l.situacao << "," << F(e.zCrit) << "," << e.base << ","
        << F(l.fsPass) << "," << F(l.tPass) << "," << F(e.frente) << "," << F(e.pwBase) << "," << F(e.pwTopo) << ","
        << e.face << "," << l.passos << "," << l.newton << "," << F(l.s) << "," << F(l.tParou) << "," << F(f.kMed) << ","
        << F(f.kMin) << "," << F(f.kMax) << "," << F(f.logPeso) << "," << op.ny << "," << F(op.dt) << ","
        << (op.quad ? "quad" : "tri") << "," << (op.concentrada ? "concentrada" : "consistente") << "," << op.leiK
        << "," << F(f.covk) << "," << F(f.lv) << "," << f.M << "," << F(f.klVar) << "\n";
}

void GravaParametros(const std::string &arq, const Opcoes &op, const RandomField *campos, int argc, char *argv[])
{
    std::ofstream p(arq);
    p << std::setprecision(12) << "# SlopeRainColumn";
    for (int i = 1; i < argc; i++) p << " " << argv[i];
    auto lista = [](const std::vector<REAL> &v) {
        std::string s;
        for (size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + G(v[i]);
        return s;
    };
    p << "\nn=" << op.n << "\nprimeira=" << op.primeira << "\ndet=" << op.det << "\nsemente=" << op.semente
      << "\nsaida=" << op.saida << "\ntempos=" << lista(op.tempos) << "\nperfis=" << op.perfis
      << "\ngrava_xi=" << op.gravaXi << "\nso_campo=" << op.soCampo << "\nkmult=" << op.kmult
      << "\nlocaliza=" << op.localiza << "\npara_na_falha=" << op.paraNaFalha << "\nH=" << op.H << "\nny=" << op.ny
      << "\nlargura=" << op.largura << "\nelementos=" << (op.quad ? "quad" : "tri")
      << "\nmassa=" << (op.concentrada ? "concentrada" : "consistente") << "\nvg=" << op.ts << "," << op.tr << ","
      << op.a << "," << op.nvg << "\nks=" << op.ks << "\nlei_k=" << op.leiK << "\nKw=" << op.Kw << "\ngw=" << op.gw
      << "\nI=" << op.I << "\nchuva=" << op.chuva;
    if (std::isnan(op.nivel0)) p << "\nsuccao0=" << op.succao0;
    else p << "\nnivel0=" << op.nivel0 << "\nsuccao_max=" << op.succaoMax;
    p << "\ndt=" << op.dt << "\ndt_min=" << op.dtMin << "\nnewton="
      << (op.newton.empty() ? std::string("padrao") : lista(op.newton)) << "\nsolver="
      << (TPZUnsatAnalysis::UsesPardiso() ? "pardiso" : "skyline") << "\nface_tol=" << op.faceTol
      << "\nfs=" << (op.fs ? "talude" : "nenhum") << "\nbeta=" << op.beta << "\nc=" << op.c << "\nphi=" << op.phi
      << "\ngamma_d=" << op.gammaD << "\nfs_peso=" << (op.pesoIntegrado ? "integrado" : "constante")
      << "\ngamma=" << op.gamma << "\nfs_zmin=" << op.fsZmin << "\ncovk=" << op.covk << "\nlv=" << op.lv
      << "\nM=" << (campos ? campos->M() : op.M) << "\nnystrom=" << (op.nystrom.empty() ? "nao" : op.nystrom)
      << "\nxi=" << (op.xi.empty() ? "nao" : op.xi) << "\ndesloc_xi=" << (op.deslocXi.empty() ? "nao" : op.deslocXi);
    if (campos) p << "\nnucleo=" << campos->KernelDescription() << "\nKL_variancia=" << campos->VarianceFraction();
    if (op.vgTab > 0.)
        p << "\ntabela=eps " << gTabelaVG.eps << " psi_c " << gTabelaVG.psiC << " kPa k_c " << gTabelaVG.kC << " pontos "
          << gTabelaVG.S.size() << " erro_max " << gTabelaVG.erroMax;
    p << "\n";
}

// Os numeros de um arquivo (xi=, desloc_xi=): ate M, um campo
std::vector<double> LeXi(const std::string &arq, int M, const char *opcao)
{
    std::ifstream f(arq);
    if (!f) Erro(std::string(opcao) + ": nao abre " + arq);
    std::vector<double> v;
    for (double x; f >> x;) v.push_back(x);
    if (v.empty() || int(v.size()) > M)
        Erro(std::string(opcao) + ": " + arq + " tem " + std::to_string(v.size()) + " numeros; precisa de 1 a M = " +
             std::to_string(M));
    return v;
}

// O cabecalho da base de Nystrom: n, M e as chaves (nucleo, lh, lv, dx, dy)
std::map<std::string, std::string> CabecalhoBase(const std::string &arq)
{
    std::ifstream f(arq);
    std::string l;
    if (!f || !std::getline(f, l)) Erro("nystrom: nao abre " + arq);
    std::istringstream s(l);
    std::map<std::string, std::string> cab;
    std::string n, M;
    s >> n >> M;
    cab["n"] = n, cab["M"] = M;
    for (std::string t; s >> t;) {
        const size_t eq = t.find('=');
        if (eq != std::string::npos) cab[t.substr(0, eq)] = t.substr(eq + 1);
    }
    return cab;
}

// ---- uma realizacao
struct Resumo {
    REAL tParou = -1., tPass = -1., simetria = 0., balanco = 0., infiltrado = 0., armazenado = 0.;
    int64_t passos = 0, newton = 0;
    double s = 0.;
    std::vector<std::pair<REAL, REAL>> fs; // (t, FS) das linhas
};

// Os fins de passo (h): k dt, os tempos e o fim da chuva (o mesmo grade() do coluna_ref.py)
std::vector<REAL> Grade(const Opcoes &op)
{
    const REAL T = op.tempos.back();
    auto arred = [](REAL x) { return std::round(x * 1.e9) / 1.e9; };
    std::set<REAL> ts;
    const int nk = int(std::ceil(T / op.dt - 1.e-9));
    for (int k = 1; k <= nk; k++) ts.insert(arred(k * op.dt));
    for (REAL t : op.tempos)
        if (t > 0.) ts.insert(arred(t));
    if (op.chuva > 0. && op.chuva < T) ts.insert(arred(op.chuva));
    std::vector<REAL> v;
    for (REAL t : ts)
        if (t <= T + 1.e-9) v.push_back(t);
    return v;
}

Resumo Realizacao(const Fixos &fixos0, const Opcoes &op, const RandomField *campos, std::ostream &csv)
{
    const auto t0 = Clock::now();
    Fixos fx = fixos0;
    const int real = fx.real;
    UnsatModel<TUnsatHP> mT;
    mT.SetGeoMesh(MalhaColuna(op.H, op.largura, op.ny, op.quad));
    // sem condicoes de contorno no modelo: o RainTransient fixa todo u (e p_a) pelas equacoes, e um elemento de contorno
    // de u so custaria montagem (o InsertBC do montecarlo.cpp levava 20% do tempo da coluna)
    mT.SetMaterial(MaterialColuna(op));
    mT.CreateCompMesh({{2, 1}, {1, 1}, {1, 1}}); // u P1 (fixo), p_w e p_a P1
    // k de cada ponto de integracao: o campo ou o kmult da deterministica
    double soma = 0.;
    int64_t npts = 0;
    fx.kMin = 1.e300, fx.kMax = -1.e300;
    mT.ForEachIntegrationPoint([&](const TPZVec<REAL> &x, TPZUnsaturatedMem &mem) {
        const REAL k = campos ? campos->At(x[0], x[1])[0] : op.kmult;
        mem.fKMult = k;
        soma += k, npts++;
        fx.kMin = std::min(fx.kMin, k), fx.kMax = std::max(fx.kMax, k);
    });
    fx.kMed = soma / npts; // os pontos tem todos o mesmo peso (malha uniforme, regras de pesos iguais)
    TPZUnsatAnalysis aT(mT.CompMesh(), {});
    if (!op.newton.empty()) aT.SetTolerances(1.e-7, op.newton[0], int(op.newton[1]), 25);
    RainTransient chuva(aT, DadosChuva(op), op.gw, {}, {topo});
    const std::string rotulo = real < 0 ? "deterministica" : "realizacao " + std::to_string(real);
    auto liga = [&](bool sim) { // newton_trace: desligado na bissecao (um passo dela que falha nao e uma parada)
        if (op.newtonTrace.empty()) return;
        chuva.SetNewtonTrace(sim ? op.newtonTrace : "", 1.e-7, op.newton.empty() ? 1.e-6 : op.newton[0],
                             op.newton.empty() ? 40 : int(op.newton[1]), 25, rotulo);
    };
    liga(true);
    if (op.faceTol > 0.) chuva.SetFaceTolerance(op.faceTol);
    chuva.InitialState();

    TPZGeoMesh *g = mT.GeoMesh();
    const Perfil P = FazPerfil(g, op);
    const TPZUnsaturatedHydraulics &h = mT.Material()->Hydraulics();
    std::vector<REAL> pwNos = chuva.NodalWaterPressure();
    std::vector<REAL> pw0(P.nos.size());
    for (size_t i = 0; i < P.nos.size(); i++) pw0[i] = pwNos[P.nos[i]];
    const REAL agua0 = AguaArmazenada(g, pwNos, h, op.ts, op.concentrada);
    REAL compress = 0.;

    // o perfil (a deterministica sempre; as realizacoes com real < perfis)
    std::unique_ptr<std::ofstream> perfil;
    if (real < 0 || real < op.perfis) {
        perfil = std::make_unique<std::ofstream>(op.saida + "_perfil_" + std::to_string(real) + ".csv");
        *perfil << "t_h,z,y,pw,Se,theta,W,FS\n" << std::setprecision(17);
    }
    Resumo res;
    REAL fsAnt = kNaN, fsHist = kInf, tPass = -1.;
    auto avalia = [&]() {
        pwNos = chuva.NodalWaterPressure();
        Estado e = Avalia(pwNos, P, h, op, pw0);
        e.face = chuva.FaceNodes();
        res.simetria = std::max(res.simetria, e.simetria);
        return e;
    };
    auto linha = [&](REAL t, const Estado &e, const std::string &situacao) {
        Linha l;
        l.t = t, l.e = e, l.situacao = situacao;
        l.fsPass = op.fs ? fsHist : kNaN;
        l.tPass = tPass >= 0. && tPass <= t + 1.e-9 ? tPass : -1.;
        l.passos = chuva.Steps(), l.newton = chuva.NewtonIterations(), l.s = Seg(t0), l.tParou = res.tParou;
        Grava(csv, l, fx, op);
        res.fs.push_back({t, e.fs});
        if (perfil && situacao == "ok")
            for (size_t i = 0; i < P.nos.size(); i++)
                *perfil << G(t) << "," << P.z[i] << "," << P.y[i] << "," << e.pw[i] << "," << e.Se[i] << "," << e.theta[i]
                        << "," << e.W[i] << "," << e.FS[i] << "\n";
    };
    auto ehTempo = [&](REAL t) {
        return std::any_of(op.tempos.begin(), op.tempos.end(), [t](REAL x) { return std::abs(x - t) < 1.e-9; });
    };

    // A primeira travessia de FS = 1 em [tL, tR] (FS(tL) = fL >= 1 > fR = FS(tR)), com a malha no estado de tR: bissecao
    // a partir do estado eL de tL (localiza=1; um passo da bissecao que falha deixa o intervalo que ja se tem) e a
    // interpolacao linear no intervalo final. A malha volta ao estado de tR.
    auto travessia = [&](RainTransient::Snapshot eL, REAL tL, REAL fL, REAL tR, REAL fR) {
        if (op.localiza) {
            const RainTransient::Snapshot fim = chuva.Save();
            liga(false);
            for (int k = 0; k < 7; k++) {
                const REAL tm = 0.5 * (tL + tR);
                chuva.Load(eL);
                if (!chuva.AdvanceTo(tm)) break; // fica o intervalo que ja se tem
                const REAL fm = Avalia(chuva.NodalWaterPressure(), P, h, op, pw0).fs;
                if (fm < 1.) tR = tm, fR = fm;
                else tL = tm, fL = fm, eL = chuva.Save();
            }
            liga(true);
            chuva.Load(fim);
        }
        return fL > fR ? tL + (fL - 1.) / (fL - fR) * (tR - tL) : tR;
    };

    Estado e = avalia();
    if (op.fs) {
        fsAnt = e.fs, fsHist = e.fs;
        if (e.fs < 1.) tPass = 0.;
    }
    std::set<REAL> feitos; // os tempos ja gravados
    if (ehTempo(0.)) linha(0., e, "ok"), feitos.insert(0.);
    REAL tcur = 0.;
    // para_na_falha com FS(0) < 1: a falha ja esta no estado inicial; nenhum passo, os instantes t > 0 ficam com falha
    bool parouNaFalha = op.fs && op.paraNaFalha && tPass == 0.;
    for (REAL tn : Grade(op)) {
        if (parouNaFalha) break;
        const RainTransient::Snapshot e0 = chuva.Save();
        const std::vector<REAL> pwAnt = pwNos;
        if (!chuva.AdvanceTo(tn)) {
            // o transiente parou: a malha guarda o ultimo sub-passo aceito, em Time() (tcur <= Time() < tn), que e um
            // estado convergido. O FS dele entra no FS_passagem e na regra da travessia, como num passo normal.
            res.tParou = chuva.Time();
            if (res.tParou > tcur + 1.e-12) {
                e = avalia();
                compress += Compressibilidade(g, pwNos, pwAnt, h, op.ts, op.Kw);
                if (op.fs) {
                    if (tPass < 0. && e.fs < 1. && fsAnt >= 1.) {
                        tPass = travessia(e0, tcur, fsAnt, res.tParou, e.fs);
                        parouNaFalha = op.paraNaFalha;
                    }
                    fsHist = std::min(fsHist, e.fs);
                    fsAnt = e.fs;
                }
            }
            break;
        }
        e = avalia();
        compress += Compressibilidade(g, pwNos, pwAnt, h, op.ts, op.Kw);
        if (op.fs) {
            if (tPass < 0. && e.fs < 1. && fsAnt >= 1.) { // a primeira travessia de FS = 1, neste passo
                tPass = travessia(e0, tcur, fsAnt, tn, e.fs);
                parouNaFalha = op.paraNaFalha;
            }
            fsHist = std::min(fsHist, e.fs);
            fsAnt = e.fs;
        }
        tcur = tn;
        if (ehTempo(tn)) linha(tn, e, "ok"), feitos.insert(tn);
        if (parouNaFalha) break;
    }
    // os instantes que faltam: o transiente parou (situacao transiente) ou parou na falha (para_na_falha)
    for (REAL t : op.tempos) {
        if (std::any_of(feitos.begin(), feitos.end(), [t](REAL x) { return std::abs(x - t) < 1.e-9; })) continue;
        if (parouNaFalha) {
            linha(t, e, "falha");
            continue;
        }
        Estado vazio;
        vazio.face = chuva.FaceNodes();
        linha(t, vazio, "transiente");
    }
    res.tPass = tPass;
    res.passos = chuva.Steps(), res.newton = chuva.NewtonIterations();
    res.infiltrado = chuva.Infiltrated();
    res.armazenado = AguaArmazenada(g, chuva.NodalWaterPressure(), h, op.ts, op.concentrada) - agua0 + compress;
    res.balanco = std::abs(res.armazenado - res.infiltrado);
    res.s = Seg(t0);
    return res;
}

// ---- so_campo=1: ln k_mult nos centros das celulas
void SoCampo(const Opcoes &op, RandomField &campos, const std::string &faixa, std::vector<std::string> &arquivos)
{
    const std::string arq = op.saida + "_campo_" + faixa + ".csv";
    std::ofstream f(arq);
    std::vector<REAL> y(op.ny);
    for (int j = 0; j < op.ny; j++) y[j] = op.H * (op.ny - 1 - j + 0.5) / op.ny; // de cima para baixo
    f << "# SlopeRainColumn so_campo=1: ln k_mult nos centros das celulas (x = largura/2), de cima para baixo; z (m):";
    for (REAL v : y) f << " " << G(op.H - v);
    f << "\nreal";
    for (int j = 0; j < op.ny; j++) f << ",lnk_" << j;
    f << "\n" << std::setprecision(10);
    // estatisticas: media e desvio em 3 profundidades; a correlacao de todos os pares na distancia lv (se for um numero
    // inteiro de celulas)
    const int L = int(std::lround(op.lv * op.ny / op.H));
    const bool temL = L >= 1 && L < op.ny && std::abs(L * op.H / op.ny - op.lv) < 1.e-9 * op.lv;
    const std::array<int, 3> js = {op.ny / 4, op.ny / 2, (3 * op.ny) / 4};
    std::vector<double> s1(op.ny, 0.), s2(op.ny, 0.);
    double sxy = 0.;
    int64_t npar = 0;
    for (int n = op.primeira; n < op.primeira + op.n; n++) {
        if (op.xi.empty()) campos.Sample(op.semente, n);
        std::vector<double> z(op.ny);
        for (int j = 0; j < op.ny; j++) z[j] = std::log(campos.At(0.5 * op.largura, y[j])[0]);
        f << n;
        for (double v : z) f << "," << v;
        f << "\n";
        for (int j = 0; j < op.ny; j++) s1[j] += z[j], s2[j] += z[j] * z[j];
        if (temL)
            for (int j = 0; j + L < op.ny; j++) sxy += z[j] * z[j + L], npar++;
    }
    const double sl = std::sqrt(std::log(1. + op.covk * op.covk)), mu = -0.5 * sl * sl, N = op.n;
    std::printf("  campo de %d realizacoes (ln k_mult; esperado: media %.5f, desvio %.5f):\n", op.n, mu, sl);
    double m1 = 0., m2 = 0.;
    for (int j = 0; j < op.ny; j++) m1 += s1[j] / N, m2 += s2[j] / N;
    m1 /= op.ny, m2 /= op.ny;
    for (int j : js) {
        const double m = s1[j] / N, sd = std::sqrt(std::max(0., s2[j] / N - m * m) * N / std::max(1., N - 1.));
        std::printf("    z = %6.4f m: media %+.5f (%+.2f erros-padrao), desvio %.5f (%+.2f%%)\n", op.H - y[j], m,
                    (m - mu) / (sl / std::sqrt(N)), sd, 100. * (sd / sl - 1.));
    }
    if (temL) {
        const double var = m2 - m1 * m1;
        std::printf("    correlacao na distancia lv = %g m (%d celulas, %lld pares): %.4f (exp(-1) = %.4f)\n", op.lv, L,
                    (long long)npar, (sxy / npar - m1 * m1) / var, std::exp(-1.));
    }
    arquivos.push_back(arq);
}

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    TPZUnsatAnalysis::SetUseSkyline(op.skyline); // antes de qualquer analise
    if (op.vgTab > 0.) { // o k_r de Mualem tabelado ate S_e = 1 - eps, conferido no material contra o EMualem
        const REAL swr = op.tr / op.ts;
        gTabelaVG = MualemTabelado(op.nvg, swr, op.a, op.vgTab, op.fecho);
        TPZUnsaturatedHydraulics h = Hidraulica(op), hm = h;
        hm.fKrw = TPZUnsaturatedHydraulics::EMualem;
        std::printf("  k_r de Mualem tabelado (fecho %s): eps = %.3g, psi_c = %.4g kPa, k_c = %.4f, %zu pontos, erro maximo "
                    "%.4f%% em [psi_c, 1e3 kPa];", NomeFechoVG(op.fecho), gTabelaVG.eps, gTabelaVG.psiC, gTabelaVG.kC,
                    gTabelaVG.S.size(), 100. * gTabelaVG.erroMax);
        for (REAL psi : {0.01, 0.1, 1., 10.})
            std::printf(" k_r(%g) = %.4f (Mualem %.4f);", psi, h.Evaluate(-psi, 0.).krw, hm.Evaluate(-psi, 0.).krw);
        std::printf("\n");
    }
    if (!op.newtonTrace.empty()) // um bloco por transiente que parar (RainTransient::SetNewtonTrace)
        std::ofstream(op.newtonTrace) << "# newton_trace do SlopeRainColumn: um bloco por transiente que parou "
                                         "(sem blocos: nenhum parou)\n";
    // o campo: a base 1D de Nystrom (lh = 1e30: o nucleo so depende de y; a grade rapida com nx = 1)
    std::unique_ptr<RandomField> campos;
    if (!op.nystrom.empty()) {
        const auto tn = Clock::now();
        const std::map<std::string, std::string> cab = CabecalhoBase(op.nystrom);
        if (cab.count("dy")) { // a base cobre [0, H]
            const REAL Hb = std::stod(cab.at("n")) * std::stod(cab.at("dy"));
            if (std::abs(Hb - op.H) > 1.e-9 * op.H)
                Erro("nystrom: a base " + op.nystrom + " cobre " + G(Hb) + " m e a coluna tem H = " + G(op.H) + " m");
        }
        campos = std::make_unique<RandomField>(0., 0., op.largura, op.H, op.lv, op.lv, op.M,
                                               std::vector<RandomField::Stats>{{1., op.covk}});
        campos->UseNystrom(op.nystrom, op.M, true);
        if (!(campos->KernelLx() > 1.e20))
            Erro("nystrom: a base " + op.nystrom + " nao e 1D (" + campos->KernelDescription() +
                 "): use nystrom_base.py coluna=H");
        if (std::abs(campos->KernelLy() - op.lv) > 1.e-9 * op.lv)
            Erro("nystrom: a base " + op.nystrom + " tem o nucleo " + campos->KernelDescription() + ", e o caso pede lv = " +
                 G(op.lv) + " (lv=)");
        if (campos->M() < op.M)
            std::printf("  aviso: a base tem %d termos (M = %d pedido)\n", campos->M(), op.M);
        if (!op.deslocXi.empty()) campos->SetShift(std::vector<std::vector<double>>{LeXi(op.deslocXi, campos->M(), "desloc_xi")});
        if (!op.xi.empty()) campos->SetXi({LeXi(op.xi, campos->M(), "xi")});
        std::printf("  base de Nystrom %s: %d termos, variancia retida %.4f, nucleo %s, %s (%.2f s)\n", op.nystrom.c_str(),
                    campos->M(), campos->VarianceFraction(), campos->KernelDescription().c_str(),
                    campos->FastNystrom() ? "grade regular (avaliacao rapida)" : "soma direta", Seg(tn));
    } else if (!op.xi.empty() || !op.deslocXi.empty())
        Erro("xi= e desloc_xi= precisam da base: nystrom=<base.txt>");
    std::printf("Coluna de chuva: H = %g m, faixa de %g m com %d celulas (%s, massa %s), solver %s | solo VG theta_s %g, "
                "theta_r %g, alfa %g 1/kPa, n %g, k_s %g m/s, k_r %s | chuva %g mm/h por %g h, %s | dt %g h (dt_min %g h)"
                " | FS %s (beta %g, c' %g kPa, phi' %g, gamma_d %g) | campo CoV %g, lv %g m, M %d | realizacoes %d..%d%s\n",
                op.H, op.largura, op.ny, op.quad ? "Q1" : "4 triangulos por celula", op.concentrada ? "concentrada" : "consistente",
                TPZUnsatAnalysis::UsesPardiso() ? "PARDISO" : "skyline", op.ts, op.tr, op.a, op.nvg, op.ks, op.leiK.c_str(),
                op.I, op.chuva,
                std::isnan(op.nivel0) ? ("succao inicial " + G(op.succao0) + " kPa").c_str()
                                      : ("nivel inicial em y = " + G(op.nivel0) + " m").c_str(),
                op.dt, op.dtMin, op.fs ? "talude" : "nenhum", op.beta, op.c, op.phi, op.gammaD, op.covk, op.lv, op.M,
                op.primeira, op.primeira + op.n - 1, op.det ? " e a deterministica" : "");
    if (op.skyline == TPZUnsatAnalysis::UsesPardiso())
        std::printf("  aviso: solver=pardiso pedido, mas o neopz foi compilado sem MKL: skyline\n");

    const std::string faixa = std::to_string(op.primeira) + "_" + std::to_string(op.primeira + op.n - 1);
    std::vector<std::string> arquivos;
    if (op.soCampo) {
        SoCampo(op, *campos, faixa, arquivos);
        if (op.primeira == 0) {
            GravaParametros(op.saida + "_parametros.txt", op, campos.get(), argc, argv);
            arquivos.push_back(op.saida + "_parametros.txt");
        }
        std::printf("  arquivos:");
        for (const std::string &a : arquivos) std::printf(" %s", a.c_str());
        std::printf("\n");
        return 0;
    }
    const std::string arq = op.saida + "_" + faixa + ".csv";
    std::ofstream csv(arq);
    csv << kCabecalho;
    arquivos.push_back(arq);
    // so o bloco que comeca na realizacao 0 (no cluster varios blocos do mesmo caso rodam juntos)
    if (op.primeira == 0) {
        GravaParametros(op.saida + "_parametros.txt", op, campos.get(), argc, argv);
        arquivos.push_back(op.saida + "_parametros.txt");
    }
    std::vector<int> reais;
    if (op.det) reais.push_back(-1);
    for (int n = op.primeira; n < op.primeira + op.n; n++) reais.push_back(n);
    double somaS = 0.;
    int nAleat = 0;
    bool detParou = false;
    for (int n : reais) {
        Fixos fx;
        fx.real = n;
        const RandomField *cf = nullptr;
        if (n >= 0 && campos) {
            if (op.xi.empty()) campos->Sample(op.semente, n);
            fx.logPeso = campos->LogWeight(), fx.covk = op.covk, fx.lv = op.lv, fx.M = campos->M();
            fx.klVar = campos->VarianceFraction();
            if (op.gravaXi) {
                std::ofstream fxi(op.saida + "_xi_" + std::to_string(n) + ".txt");
                arquivos.push_back(op.saida + "_xi_" + std::to_string(n) + ".txt");
                fxi << std::setprecision(17);
                for (double x : campos->Xi()[0]) fxi << x << "\n";
            }
            cf = campos.get();
        }
        const Resumo r = Realizacao(fx, op, cf, csv);
        csv.flush();
        if (n < 0 && r.tParou >= 0.) detParou = true;
        if (n >= 0) somaS += r.s, nAleat++;
        if (n < 0 || n < op.perfis) arquivos.push_back(op.saida + "_perfil_" + std::to_string(n) + ".csv");
        std::printf("  %s: %.2f s, %lld passos, %lld iteracoes de Newton, max |p_L - p_R| = %.1e kPa, |dW - infiltrado| = "
                    "%.2e m3/m (%.1e do infiltrado %.4g m3/m)%s; FS",
                    n < 0 ? "deterministica" : ("realizacao " + std::to_string(n)).c_str(), r.s, (long long)r.passos,
                    (long long)r.newton, r.simetria, r.balanco, r.infiltrado > 0. ? r.balanco / r.infiltrado : 0.,
                    r.infiltrado, r.tPass >= 0. ? (", t_passagem " + G(r.tPass) + " h").c_str() : "");
        for (const auto &[t, f] : r.fs) std::printf(" %g h: %.4f", t, f);
        if (r.tParou >= 0.) std::printf(" | o transiente parou em %.4f h", r.tParou);
        std::printf("\n");
        std::fflush(stdout);
    }
    if (!op.newtonTrace.empty()) arquivos.push_back(op.newtonTrace);
    if (nAleat) std::printf("  tempo medio por realizacao: %.2f s (%d realizacoes)\n", somaS / nAleat, nAleat);
    std::printf("  arquivos:");
    for (const std::string &a : arquivos) std::printf(" %s", a.c_str());
    std::printf("\n");
    return op.det && op.n == 0 && detParou ? 2 : 0;
}
