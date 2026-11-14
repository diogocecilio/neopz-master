// Monte Carlo do talude de Rahardjo et al. (2007) sob chuva, com campos aleatorios lognormais de c', phi' e k_s
// (RandomField: Karhunen-Loeve analitico; padrao: Tabela 2 de Vargas Ceron et al. 2025). Talude da Fig. 2 do Rahardjo
// (H = 10 m a 45 graus) na TriGMesh do main.cpp. Em cada realizacao:
//  1. transiente da chuva (RainTransient: esqueleto rigido, ar passivo) na malha do transiente, com o k de cada ponto do
//     campo, ate o ultimo instante pedido; guarda p_w em cada instante do fator;
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
//                             chuva (nivel d'agua em x = 30, 35 e 40 m, p_w na banda, nos na face, ...); niveis: o FS
//                             de cada nivel do refinamento
//   <saida>_parametros.txt    chave=valor do caso
//   <saida>_mapa_<faixa>.csv  por instante e celula de 1 m: quantas realizacoes tem a banda do mecanismo ali (n) e quantas
//                             dessas tem FS < 1 (n_falha); n_real: realizacoes somadas (somar os arquivos das partes)
//   <saida>_real<n>/          campos de uma realizacao: o formato do LearningMohrCoulomb2/animacao.py, campos.csv (c', phi',
//                             k nos nos) e mecanismo.csv (pontos com sqrt(J2) >= 1% do maior, como no rebaixamento)
#include "AnimationFiles.h"
#include "HPMaterial.h"
#include "ImportanceSampling.h"
#include "LayerGMesh.h"
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

// talude da Fig. 2 de Rahardjo et al. (2007), H = 10 m a 45 graus: cantos no sentido anti-horario (base esquerda, base
// direita, terreno do pe, pe, crista, crista a esquerda), como no main.cpp
const std::vector<std::vector<REAL>> kTalude = {{0, 0}, {70, 0}, {70, 20}, {40, 20}, {30, 30}, {0, 30}};

// ---- dados do problema
struct Solo { // Tabela 1 do Rahardjo: o que muda de um solo para outro
    std::string nome;
    REAL a, ks; // parametro a de Fredlund-Xing (kPa) e k_s (m/s)
};
const std::vector<Solo> kSolos = {{"f10", 10., 1.e-4}, {"f50", 50., 1.e-5}, {"f100", 100., 1.e-6}};
struct Resistencia { // igual nos tres solos
    REAL c = 10., phiDeg = 26.; // medias de c' (kPa) e phi' (graus)
    REAL gamma = 20., gammaw = 9.81, E = 2.e4, nu = 0.3;
};
struct Campos { // Tabela 2 do Ceron: CoV e distancias de autocorrelacao (m); k: multiplicador de k_s com media 1
    REAL covC = 0.3, covPhi = 0.1, covK = 0.6, Lx = 20., Ly = 2.;
    int M = 2000;
    uint64_t semente = 2025; // a realizacao n tem os mesmos numeros aleatorios em todos os casos
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
    std::string banda;         // amostragem por importancia: o mecanismo.csv da banda e o instante dela
    REAL bandaT = -1.;
    ai::Opcoes ai;             // o deslocamento pelo piloto
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

Opcoes LeOpcoes(int argc, char *argv[])
{
    const std::set<std::string> chaves = {"n",      "primeira", "det",  "solo", "I",          "chuva", "hw",
                                          "tempos", "adapta",   "metodo", "covc", "covphi",   "covk",  "s",
                                          "alfa",   "transiente", "uT", "fator", "nref",      "figuras", "saida", "valores",
                                          "semente", "banda", "banda_t", "piloto", "piloto_n", "piloto_semente",
                                          "ai", "elite", "modos"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        if (eq == std::string::npos || !chaves.count(a.substr(0, eq))) {
            std::fprintf(stderr, "opcao desconhecida: %s (as opcoes estao no cabecalho do montecarlo.cpp)\n", a.c_str());
            std::exit(1);
        }
        kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    auto tem = [&kv](const char *k) { return kv.count(k) > 0; };
    Opcoes op;
    if (tem("n")) op.n = std::stoi(kv["n"]);
    if (tem("primeira")) op.primeira = std::stoi(kv["primeira"]);
    if (tem("det")) op.det = std::stoi(kv["det"]) != 0;
    if (tem("solo")) {
        auto it = std::find_if(kSolos.begin(), kSolos.end(), [&](const Solo &s) { return s.nome == kv["solo"]; });
        if (it == kSolos.end()) {
            std::fprintf(stderr, "solo desconhecido: %s (use f10, f50 ou f100)\n", kv["solo"].c_str());
            std::exit(1);
        }
        op.solo = *it;
    }
    op.I = tem("I") ? std::stod(kv["I"]) : op.solo.ks * 3.6e6; // padrao: 1 k_s (m/s -> mm/h)
    if (tem("chuva")) op.chuva = std::stod(kv["chuva"]);
    if (tem("hw")) op.hw = std::stod(kv["hw"]);
    if (tem("tempos")) op.tempos = Lista(kv["tempos"]);
    std::sort(op.tempos.begin(), op.tempos.end());
    op.adapta = tem("adapta") ? Lista(kv["adapta"]) : op.tempos;
    for (REAL t : op.adapta)
        if (std::find(op.tempos.begin(), op.tempos.end(), t) == op.tempos.end()) {
            std::fprintf(stderr, "adapta: %g nao esta em tempos\n", t);
            std::exit(1);
        }
    if (op.adapta.empty()) {
        std::fprintf(stderr, "adapta precisa de pelo menos um instante\n");
        std::exit(1);
    }
    if (tem("metodo")) {
        if (kv["metodo"] != "srm" && kv["metodo"] != "gim") {
            std::fprintf(stderr, "metodo: srm ou gim\n");
            std::exit(1);
        }
        op.gim = kv["metodo"] == "gim";
    }
    if (tem("covc")) op.campos.covC = std::stod(kv["covc"]);
    if (tem("covphi")) op.campos.covPhi = std::stod(kv["covphi"]);
    if (tem("covk")) op.campos.covK = std::stod(kv["covk"]);
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
    op.ai.semente = tem("piloto_semente") ? std::stoull(kv["piloto_semente"]) : Campos{}.semente;
    if (tem("piloto")) op.ai.piloto = kv["piloto"];
    if (tem("piloto_n")) op.ai.n = std::stoi(kv["piloto_n"]);
    if (tem("ai")) op.ai.tipo = kv["ai"];
    if (tem("elite")) op.ai.elite = std::stod(kv["elite"]);
    if (tem("modos")) op.ai.modos = std::stoi(kv["modos"]);
    if (op.ai.tipo != "banda" && op.ai.tipo != "modos") {
        std::fprintf(stderr, "ai=%s: banda ou modos\n", op.ai.tipo.c_str());
        std::exit(1);
    }
    if (!op.ai.piloto.empty() && op.ai.tipo == "banda" && op.banda.empty()) {
        std::fprintf(stderr, "piloto com ai=banda precisa de banda=<mecanismo.csv>\n");
        std::exit(1);
    }
    if (tem("valores")) {
        op.valores = Lista(kv["valores"]);
        if (op.valores.size() != 3) {
            std::fprintf(stderr, "valores: c,phi,k (c' em kPa, phi' em graus, k / k_s)\n");
            std::exit(1);
        }
    }
    if (op.malhaF.rfind("camq", 0) == 0) {
        std::fprintf(stderr, "a malha do fator precisa de triangulos (a bisseccao do refinamento): use triK ou camH\n");
        std::exit(1);
    }
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
    const auto &p = kTalude;
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
    if (nome.rfind("tri", 0) == 0) return MalhaPlana(TriGMesh(kTalude, std::stoi(nome.substr(3))));
    std::fprintf(stderr, "malha desconhecida: %s (use triK, camH[pD] ou camqH[pD])\n", nome.c_str());
    std::exit(1);
}

// cota da superficie do talude em x
REAL Superficie(REAL x)
{
    const auto &p = kTalude;
    if (x <= p[4][0]) return p[4][1];
    if (x >= p[3][0]) return p[3][1];
    return p[4][1] + (p[3][1] - p[4][1]) * (x - p[4][0]) / (p[3][0] - p[4][0]);
}

// ---- material do solo (o MaterialRahardjo do rahardjo.cpp), com a anisotropia k_h / k_v
TUnsatHP *MaterialChuva(const RainData &r, REAL alfa)
{
    const Resistencia s;
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EFredlundXing; // S_w(p_c) com C(psi) = 1
    h.fA = r.aFX, h.fN = r.nFX, h.fM = r.mFX;
    h.fKrw = TPZUnsaturatedHydraulics::EPowerLaw; // k_w = k_s Theta^p
    h.fPower = r.power;
    h.fKra = TPZUnsaturatedHydraulics::EBrooksCorey; // o ar fica passivo (p_a = 0)
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, std::tan(r.phib) / std::tan(phi)); // (u_a - u_w) tan(phi^b)
    mat->SetPorosity(r.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(r.Kw);
    mat->SetConductivity(r.ks, s.gammaw, 1.e-3 / 1.8e-5);
    if (alfa != 1.) mat->SetAnisotropicConductivity(alfa * r.ks, r.ks); // k_x = alfa k_s, k_y = k_s
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
};

void NovoFator(Fator &f, TPZGeoMesh *g, const RainData &r, const Opcoes &op)
{
    f.a.reset();
    f.m = std::make_unique<UnsatModel<TUnsatHP>>();
    f.m->SetGeoMesh(g);
    InsertBC(*f.m);
    TUnsatHP *mat = MaterialChuva(r, op.alfa);
    if (op.gim) mat->SetGravityScalesSuction(false); // GIM: a succao nao e peso, fica fixa
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
    f.a->SetVerbose(false);
    f.a->SetLoadFactorFunction([mat](REAL l) { mat->SetGravityFactor(l); }, false);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
}

// FS (ou Gamma) no modelo de f com as pressoes pw (nos da malha do transiente), do estado sem tensoes, com os campos
REAL CalculaFator(Fator &f, const Campo *campos, const PointLocator &loc, const std::vector<REAL> &pw, bool gim)
{
    f.m->Material()->ResetMemory(); // estado sem tensoes (apaga tambem c, phi e k dos pontos)
    if (campos) AplicaCampos(*f.m, *campos);
    f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
        return espaco == 1 ? loc.Interpolate(x[0], x[1], pw) : 0.; // p_w do transiente; p_a = 0
    });
    return gim ? f.a->GIMArcLength() : f.a->SRMArcLength();
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
        const REAL fs = CalculaFator(f, campos, loc, pw, op.gim);
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

Mecanismo DescreveMecanismo(Fator &f, const PointLocator &loc, const std::vector<REAL> &pw)
{
    const Resistencia res;
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
    REAL y = kTalude[0][1];
    while (y + 0.1 <= topo && loc.Interpolate(x, y + 0.1, pw) >= 0.) y += 0.1;
    return y;
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

// colunas: as do mc_resumo.csv do rebaixamento onde existem (FS_min, banda_*, massa_*, *_med, cov_*, Lx, Ly, alfa)
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
void GravaParametros(const std::string &arq, const Opcoes &op, const RandomField &campos)
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
      << "\nAI=" << (op.ai.piloto.empty() ? "nenhuma" : op.ai.tipo + "_piloto") << "\n";
    if (!op.ai.piloto.empty())
        p << "PILOTO=" << op.ai.piloto << "\nPILOTO_N=" << op.ai.n << "\nPILOTO_SEMENTE=" << op.ai.semente
          << "\nELITE=" << op.ai.elite << "\nMODOS=" << op.ai.modos << "\nBANDA=" << op.banda << "\nBANDA_T=" << op.bandaT
          << "\n";
}

// mapa da faixa do mecanismo: celulas de 1 m no retangulo do talude, por instante
struct Mapa {
    static constexpr int nx = 70, ny = 30;
    std::map<REAL, std::vector<int>> n, nFalha;
    int nReal = 0;
    void Soma(REAL t, const Mecanismo &d, bool falha)
    {
        std::vector<int> &v = n[t], &vf = nFalha[t];
        v.resize(nx * ny, 0), vf.resize(nx * ny, 0);
        std::set<int> celulas;
        for (const auto &p : d.banda)
            celulas.insert(std::clamp(int(p[1]), 0, ny - 1) * nx + std::clamp(int(p[0]), 0, nx - 1));
        for (int c : celulas) v[c]++, vf[c] += falha;
    }
    void Grava(const std::string &arq) const
    {
        std::ofstream f(arq);
        f << "t_h,x,y,n,n_falha,n_real\n";
        for (const auto &[t, v] : n)
            for (int c = 0; c < nx * ny; c++)
                if (v[c]) f << t << "," << c % nx + 0.5 << "," << c / nx + 0.5 << "," << v[c] << "," << nFalha.at(t)[c] << ","
                            << nReal << "\n";
    }
};

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    const Resistencia res;
    RainData r; // Fredlund-Xing n = m = 1, theta_s = 0,45, k_w = k_s Theta^4, phi^b = 26, succao ate 75 kPa
    r.aFX = op.solo.a, r.ks = op.solo.ks;
    r.intensity = op.I, r.duration = op.chuva, r.Hw = op.hw;
    r.xToe = kTalude[3][0], r.yToe = kTalude[3][1];
    const Campos &cp = op.campos;
    RandomField campos(0., 0., kTalude[1][0], kTalude[4][1], cp.Lx, cp.Ly, cp.M,
                       {{res.c, cp.covC}, {res.phiDeg, cp.covPhi}, {1., cp.covK}}); // c' (kPa), phi' (graus), k / k_s
    const Campo campoAleatorio = [&campos](double x, double y) { return campos.At(x, y); };
    const Campo campoFixo = [v = op.valores](double, double) { return std::vector<double>(v.begin(), v.end()); };
    std::printf("Monte Carlo com chuva (Rahardjo, Fig. 2): solo %s, I = %g mm/h por %g h, nivel %g m abaixo do pe | %s | "
                "realizacoes %d..%d%s\n  campos: CoV c' %g, phi' %g, k %g; Lx = %g m, Ly = %g m; alfa = %g; KL: %d termos, "
                "%.1f%% da variancia\n  malhas: transiente %s (u P%d), fator %s + %d niveis\n",
                op.solo.nome.c_str(), op.I, op.chuva, op.hw, op.gim ? "GIM" : "SRM", op.primeira, op.primeira + op.n - 1,
                op.det ? " e a deterministica" : "", cp.covC, cp.covPhi, cp.covK, cp.Lx, cp.Ly, op.alfa, campos.M(),
                100. * campos.VarianceFraction(), op.malhaT.c_str(), op.uT, op.malhaF.c_str(), op.nref);

    const std::string faixa = std::to_string(op.primeira) + "_" + std::to_string(op.primeira + op.n - 1);
    // amostragem por importancia: a direcao da banda (projecoes) e o deslocamento pelo piloto
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
    const std::string arq = op.saida + "_" + faixa + ".csv";
    std::ofstream csv(arq);
    csv << kCabecalho;
    // so o bloco que comeca na realizacao 0: no cluster, varios processos rodam blocos do mesmo caso ao mesmo tempo
    if (op.primeira == 0) GravaParametros(op.saida + "_parametros.txt", op, campos);
    Mapa mapa;
    std::map<REAL, std::vector<REAL>> fsInstante; // so as realizacoes aleatorias
    std::map<REAL, double> pesoFalha;             // AI: soma de peso * (FS < 1) por instante
    double pesoFalhaMin = 0.;                     // e com o menor FS dos instantes
    int nPesadas = 0;                             // realizacoes com todos os instantes ok
    double somaT = 0., somaF = 0.;

    std::vector<int> realizacoes;
    if (op.det) realizacoes.push_back(-1);
    for (int n = op.primeira; n < op.primeira + op.n; n++) realizacoes.push_back(n);

    for (int n : realizacoes) {
        const auto t0 = Clock::now();
        const bool deterministica = n < 0;
        if (!deterministica) campos.Sample(cp.semente, n);
        // o peso da amostragem por importancia (0 sem deslocamento) e as projecoes na direcao da banda
        const double logPeso = deterministica ? 0. : campos.LogWeight();
        const std::vector<double> proj = deterministica ? std::vector<double>(3, 0.) : campos.Projection(dir);
        // propriedades de cada ponto: o campo da realizacao; na deterministica, as medias ou os valores dados
        const Campo *cf = deterministica ? (op.valores.empty() ? nullptr : &campoFixo) : &campoAleatorio;

        // 1. transiente ate o ultimo instante, guardando p_w nos instantes do fator
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        InsertBC(mT);
        mT.SetMaterial(MaterialChuva(r, op.alfa));
        mT.CreateCompMesh({{2, op.uT}, {1, 1}, {1, 1}});
        if (cf) AplicaCampos(mT, *cf);
        TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
        TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
        strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
        aT.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        aT.SetSolver(solver);
        RainTransient chuva(aT, r, res.gammaw, {left, right}, {topleft, face, topright});
        chuva.InitialState();
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
        std::map<REAL, std::vector<REAL>> pws;
        double sT = 0.;
        bool parou = false;
        for (REAL t : op.tempos) {
            const auto tt = Clock::now();
            if (!parou && !chuva.AdvanceTo(t)) {
                std::printf("  realizacao %d: o transiente parou em %.2f h\n", n, chuva.Time());
                parou = true;
            }
            sT += Seg(tt);
            if (parou) { // os instantes que faltam ficam no arquivo, sem fator
                Linha &l = linhas[t];
                l.real = n, l.t = t, l.situacao = "transiente", l.sT = sT;
                continue;
            }
            Linha &l = linhas[t];
            l.real = n, l.t = t, l.sT = sT, l.eqT = mT.CompMesh()->NEquations();
            pws[t] = chuva.NodalWaterPressure();
            const auto [pmin, pmax] = chuva.PressureRange();
            l.pwMin = pmin, l.pwMax = pmax, l.face = chuva.FaceNodes();
            l.nivel30 = NivelAgua(loc, pws[t], 30.), l.nivel35 = NivelAgua(loc, pws[t], 35.);
            l.nivel40 = NivelAgua(loc, pws[t], 40.);
            if (anim) anim->Snapshot(t, pws[t], mT.Material()->Hydraulics());
        }

        // 2. fator: primeiro os instantes adaptativos (nas realizacoes das figuras, todos), depois os outros na malha do
        // adaptativo mais proximo
        std::map<REAL, Fator> adaptados;
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
            l.mec = DescreveMecanismo(f, loc, pws[t]);
            if (anim) anim->FactorOfSafety(t, fs, *f.m, *f.a);
            if (mecanismo)
                for (const auto &p : l.mec.pontos)
                    *mecanismo << t << "," << p.x << "," << p.y << "," << p.j << "," << p.du << "," << p.c << ","
                               << p.phi << "," << p.k << "\n";
        };
        double sF = 0.;
        for (REAL t : op.tempos) {
            if (!pws.count(t)) continue; // o transiente parou antes
            const bool adapt = figura || std::find(op.adapta.begin(), op.adapta.end(), t) != op.adapta.end();
            if (!adapt) continue;
            const auto tf = Clock::now();
            std::vector<REAL> niveis;
            const REAL fs = FatorAdaptativo(adaptados[t], r, op, cf, loc, pws[t], niveis);
            const double s = Seg(tf);
            sF += s;
            registra(t, fs, adaptados[t], niveis, true, s);
        }
        for (REAL t : op.tempos) {
            if (!pws.count(t) || adaptados.count(t) || adaptados.empty()) continue;
            REAL perto = adaptados.begin()->first; // o instante adaptativo mais proximo
            for (const auto &[ta, f] : adaptados)
                if (std::abs(ta - t) < std::abs(perto - t)) perto = ta;
            const auto tf = Clock::now();
            const REAL fs = CalculaFator(adaptados[perto], cf, loc, pws[t], op.gim);
            const double s = Seg(tf);
            sF += s;
            registra(t, fs, adaptados[perto], {fs}, false, s);
        }

        // 3. arquivo, mapa e resumo da realizacao
        bool todosOk = true;
        REAL fsMin = 1e30;
        for (const auto &[t, l] : linhas) {
            Grava(csv, l, op, campos.VarianceFraction(), logPeso, proj);
            if (!deterministica && l.situacao == "ok") {
                fsInstante[t].push_back(l.fs);
                mapa.Soma(t, l.mec, l.fs < 1.);
            }
            todosOk = todosOk && l.situacao == "ok", fsMin = std::min(fsMin, l.fs);
        }
        if (!deterministica && todosOk) {
            nPesadas++;
            for (const auto &[t, l] : linhas) pesoFalha[t] += (l.fs < 1.) * std::exp(logPeso);
            pesoFalhaMin += (fsMin < 1.) * std::exp(logPeso);
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
    if (!op.ai.piloto.empty() && nPesadas) { // a P_f da amostragem por importancia, so deste bloco
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
