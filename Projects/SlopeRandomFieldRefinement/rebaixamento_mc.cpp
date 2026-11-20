// Monte Carlo do rebaixamento rapido do talude de Vargas Ceron et al. (2025), Secoes 6.1 a 6.5, com o transiente: o
// reservatorio desce da crista ate o nivel final em 0,1 h (k_v = 1e-6 m/s), e Gamma (GIM, como no artigo) e calculado no
// fim do rebaixamento, com a pressao do reservatorio sobre a face e o terreno do pe crescendo com lambda como o peso. E o
// TransientDrawDown do main.cpp (Gamma = 1,3339 no talude ceron2 com 3 niveis) com campos
// aleatorios lognormais de c', phi' e k_v (RandomField: Karhunen-Loeve analitico, Tabela 2 do artigo). Em cada
// realizacao:
//  1. transiente na malha do transiente (esqueleto rigido, ar passivo; o k de cada ponto do campo) ate o fim do
//     rebaixamento;
//  2. Gamma com as pressoes levadas para a malha do fator (PointLocator + SetFrozenPressures), com o refinamento
//     adaptativo em h feito do zero (malha grossa e nref niveis);
//  3. descritores do mecanismo: banda (sqrt(J2(eps_p)) >= 0,1 do maior), massa que desliza (|du| do ultimo incremento
//     do caminho >= 0,5 do maior, como na campanha do rebaixamento de antes) e as medias dos campos.
// Uso: SlopeDrawdownMonteCarlo [chave=valor ...], com (padrao entre parenteses):
//   n (10), primeira (0)   realizacoes primeira..primeira+n-1
//   det (0)                1: antes, o caso deterministico (as medias), como a realizacao -1
//   valores ()             c',phi',k: a deterministica com esses valores em todo o solo (c' kPa, phi' graus, k / k_v)
//   hw (1)                 h_w / H: o rebaixamento (1: total, ate o pe; Secao 6.5 do artigo: 0,5 a 0,9)
//   trebaix (0.1)          duracao do rebaixamento (h)
//   covc (0.3), covphi (0.1), covk (0.6)   CoV de c', phi' e k_v (Secao 6.2)
//   s (1)                  escala das distancias de autocorrelacao: 20 s e 2 s m (Secao 6.3)
//   alfa (1)               k_h / k_v (Secao 6.4)
//   transiente (tri3)      malha do transiente (TriGMesh com K divisoes); fator (tri1), nref (3): a do fator
//   figuras ()             realizacoes com os campos gravados (pasta <saida>_real<n>), ex.: -1,0
// Amostragem por importancia (o estimador ponderado de P_f e sem vies; posproc/amostragem_importancia.py):
//   banda (<arquivo>)      o mecanismo.csv do deterministico: os pontos com sqrt(J2) >= 0,1 do maior dao a direcao
//                          d dos xi (o gradiente da media do campo gaussiano na banda)
//   desloc (s_c,s_phi,s_k) desloca a media dos xi de cada campo para s d; o arquivo ganha log_peso
//   semente (20250416)     semente dos campos (outra para uma amostra independente de outra rodada)
//   media_falhas (<arquivo>), modos (30)   so a media dos modos 1..K de cada campo nas realizacoes listadas no
//                          arquivo (um numero por linha; as que falharam no piloto): o deslocamento da entropia
//                          cruzada regularizado nos K primeiros modos, gravado em <saida>_desloc_modos.csv
//   desloc_modos (<arquivo>) desloca a media dos xi pelo <saida>_desloc_modos.csv (termo,c,phi,k)
//   so_projecoes (0)       1: so sorteia os campos e grava as projecoes xi . d e o log_peso (<saida>_proj_<faixa>.csv)
//                          e a direcao d (<saida>_direcao.csv), sem elementos finitos (para conferir os pesos e estimar
//                          o deslocamento pelas falhas de um piloto)
//   piloto (<pasta>)       o deslocamento sai de um Monte Carlo direto ja feito (o piloto: as realizacoes
//                          0..piloto_n-1 dos arquivos <prefixo>_<a>_<b>.csv da pasta), pela entropia cruzada: a media
//                          dos xi do conjunto elite (as que falharam ou, se forem menos que a fracao elite das validas,
//                          as de menor Gamma), sorteados de novo com a semente do piloto. Todos os blocos que leem o
//                          mesmo piloto completo calculam o mesmo deslocamento (gravado em <saida>_desloc_<faixa>.txt).
//                          Sai com o codigo 3 se o piloto tiver menos de 90% das realizacoes
//   piloto_n (1000), piloto_semente (20250416), elite (0.1)
//   ai (banda)             com piloto: banda (s_c = media de xi_c . d, a direcao da banda; precisa de banda=) ou modos
//                          (a media dos modos 1..K de cada campo, K = modos)
//   saida (mc)             prefixo dos arquivos
//   xi (<arquivo>)         os xi dados, sem sorteio (um ponto do FORM, Projects/Cho2010/form_kl.py prog=rebaixamento):
//                          3K numeros, os K de c', os K de phi' e os K de k (K <= M; os outros termos ficam em 0);
//                          a realizacao primeira usa esses xi
//   gradiente (0)          K > 0: depois de Gamma, o gradiente dGamma/dxi nos K primeiros termos de cada campo pela
//                          projecao modal da sensibilidade de limite superior (GIM, Mohr-Coulomb associado, tensoes
//                          efetivas do colapso com as pressoes do transiente congeladas): <saida>_grad_<real>.txt (Gamma
//                          e os 3K numeros, c', phi' e k). O k so muda as pressoes, que a analise limite congela: os K
//                          termos de k saem 0
// Saidas (<faixa> = primeira_ultima): <saida>_<faixa>.csv (uma linha por realizacao, com as colunas do mc_resumo.csv da
// campanha do rebaixamento: FS_min = Gamma, banda_*, massa_*, *_med, parametros), <saida>_mapa_<faixa>.csv (frequencia
// da banda em celulas de 1 m), <saida>_parametros.txt e <saida>_real<n>/ (animacao, campos.csv e mecanismo.csv).
// posproc/junta_rebaixamento.py monta com elas as pastas que o SlopeRandomField/comparacao_ceron/relatorio_sat.py le.
//
// [Atencao, 2026-11-20 (Projects/CeronCinematico): (1) o 1,336 do artigo e de percolacao PERMANENTE saturada logo apos
//  um rebaixamento instantaneo (Eqs. 19-21), nao deste transiente nao saturado de 0,1 h: os dois nao se comparam. (2) o
//  GIM com 3 niveis nao esta convergido: na exportacao det, NREF=3 da 1,3253, NREF=5 da 1,2851 e o extrapolado e ~1,276
//  (o cinematico da 1,279); o acordo antigo 1,3339 x 1,336 foi coincidencia de malha grossa. (3) kv abaixo e usado
//  como condutividade em m/s; a Tabela 2 do artigo da k_v = 1e-6 m^4/(kN s), isto e k = 9,81e-6 m/s.]
#include "AnimationFiles.h"
#include "DrawdownTransient.h"
#include "HPMaterial.h"
#include "PointLocator.h"
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
#include "pzmultiphysicselement.h"
#include "pzstepsolver.h"
#include "tpzgeoelrefpattern.h"

#include <dirent.h>

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
// propriedades num ponto (x, y): {c' (kPa), phi' (graus), k / k_v}
using Campo = std::function<std::vector<double>(double, double)>;
static double Seg(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// ids da TriGMesh do main
enum EIds { matSolo = 1, base = -1, right = -2, topright = -3, topleft = -4, left = -5, face = -6 };

// o talude do Ceron (H = 5 m a 45 graus) no dominio do main.cpp (ceron2): base esquerda, base direita, terreno do pe, pe,
// crista e crista a esquerda
const std::vector<std::vector<REAL>> kTalude = {{0, 0}, {70, 0}, {70, 30}, {40, 30}, {35, 35}, {0, 35}};

// ---- dados (o UnsatParams do main.cpp: Tabela 2 do artigo e o material nao saturado padrao)
struct Solo {
    REAL c = 10., phiDeg = 30.; // medias de c' (kPa) e phi' (graus)
    REAL gamma = 20., gammaw = 9.81, E = 2.e4, nu = 0.3;
    REAL kv = 1.e-6, n = 0.45, Kw = 2.2e6; // k_v (m/s), porosidade, modulo da agua (kPa)
};
struct Campos { // Tabela 2 do Ceron
    REAL covC = 0.3, covPhi = 0.1, covK = 0.6, Lx = 20., Ly = 2.;
    int M = 2000;
    uint64_t semente = 20250416; // a da campanha do rebaixamento de antes
};

struct Opcoes {
    int n = 10, primeira = 0;
    bool det = false;
    REAL hw = 1., trebaix = 0.1, alfa = 1., s = 1.;
    Campos campos;
    std::string malhaT = "tri3", malhaF = "tri1", saida = "mc";
    int nref = 3;
    std::set<int> figuras;
    std::vector<REAL> valores;
    std::string banda;          // amostragem por importancia: o mecanismo.csv da banda
    std::vector<REAL> desloc;   // s_c, s_phi, s_k
    std::string mediaFalhas, deslocModos; // deslocamento nos K primeiros modos
    int modos = 30;
    bool soProjecoes = false;
    std::string piloto, ai = "banda"; // deslocamento pelo piloto: a pasta do Monte Carlo direto e o tipo
    int pilotoN = 1000;
    uint64_t pilotoSemente = Campos{}.semente;
    double elite = 0.1;
    std::string xi;    // xi=: os xi dados (FORM)
    int gradiente = 0; // gradiente=K: o gradiente modal de Gamma nos K primeiros termos
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
    const std::set<std::string> chaves = {"n",    "primeira", "det",        "valores", "hw",   "trebaix", "covc",
                                          "covphi", "covk",   "s",          "alfa",    "transiente", "fator", "nref",
                                          "figuras", "saida", "banda", "desloc", "so_projecoes", "semente", "media_falhas", "modos", "desloc_modos",
                                          "piloto", "piloto_n", "piloto_semente", "ai", "elite", "xi", "gradiente"};
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        const size_t eq = a.find('=');
        if (eq == std::string::npos || !chaves.count(a.substr(0, eq))) {
            std::fprintf(stderr, "opcao desconhecida: %s (as opcoes estao no cabecalho do rebaixamento_mc.cpp)\n", a.c_str());
            std::exit(1);
        }
        kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    auto tem = [&kv](const char *k) { return kv.count(k) > 0; };
    Opcoes op;
    if (tem("n")) op.n = std::stoi(kv["n"]);
    if (tem("primeira")) op.primeira = std::stoi(kv["primeira"]);
    if (tem("det")) op.det = std::stoi(kv["det"]) != 0;
    if (tem("hw")) op.hw = std::stod(kv["hw"]);
    if (tem("trebaix")) op.trebaix = std::stod(kv["trebaix"]);
    if (tem("covc")) op.campos.covC = std::stod(kv["covc"]);
    if (tem("covphi")) op.campos.covPhi = std::stod(kv["covphi"]);
    if (tem("covk")) op.campos.covK = std::stod(kv["covk"]);
    if (tem("s")) op.s = std::stod(kv["s"]);
    op.campos.Lx = 20. * op.s, op.campos.Ly = 2. * op.s;
    if (tem("alfa")) op.alfa = std::stod(kv["alfa"]);
    if (tem("transiente")) op.malhaT = kv["transiente"];
    if (tem("fator")) op.malhaF = kv["fator"];
    if (tem("nref")) op.nref = std::stoi(kv["nref"]);
    if (tem("figuras"))
        for (REAL v : Lista(kv["figuras"])) op.figuras.insert(int(v));
    if (tem("saida")) op.saida = kv["saida"];
    if (tem("banda")) op.banda = kv["banda"];
    if (tem("semente")) op.campos.semente = std::stoull(kv["semente"]);
    if (tem("media_falhas")) op.mediaFalhas = kv["media_falhas"];
    if (tem("modos")) op.modos = std::stoi(kv["modos"]);
    if (tem("desloc_modos")) op.deslocModos = kv["desloc_modos"];
    if (tem("desloc")) op.desloc = Lista(kv["desloc"]);
    if (tem("so_projecoes")) op.soProjecoes = std::stoi(kv["so_projecoes"]) != 0;
    if (tem("piloto")) op.piloto = kv["piloto"];
    if (tem("piloto_n")) op.pilotoN = std::stoi(kv["piloto_n"]);
    if (tem("piloto_semente")) op.pilotoSemente = std::stoull(kv["piloto_semente"]);
    if (tem("ai")) op.ai = kv["ai"];
    if (tem("elite")) op.elite = std::stod(kv["elite"]);
    if (tem("xi")) op.xi = kv["xi"];
    if (tem("gradiente")) op.gradiente = std::stoi(kv["gradiente"]);
    if (int(!op.desloc.empty()) + int(!op.deslocModos.empty()) + int(!op.piloto.empty()) > 1) {
        std::fprintf(stderr, "desloc, desloc_modos e piloto: so um deles\n");
        std::exit(1);
    }
    if (op.ai != "banda" && op.ai != "modos") {
        std::fprintf(stderr, "ai=%s: banda ou modos\n", op.ai.c_str());
        std::exit(1);
    }
    if ((!op.desloc.empty() || op.soProjecoes || (!op.piloto.empty() && op.ai == "banda")) && op.banda.empty()) {
        std::fprintf(stderr, "desloc, so_projecoes e piloto com ai=banda precisam de banda=<mecanismo.csv>\n");
        std::exit(1);
    }
    if (tem("valores")) {
        op.valores = Lista(kv["valores"]);
        if (op.valores.size() != 3) {
            std::fprintf(stderr, "valores: c,phi,k (c' em kPa, phi' em graus, k / k_v)\n");
            std::exit(1);
        }
    }
    for (const std::string *m : {&op.malhaT, &op.malhaF})
        if (m->rfind("tri", 0) != 0) {
            std::fprintf(stderr, "malha desconhecida: %s (use triK)\n", m->c_str());
            std::exit(1);
        }
    return op;
}

// ---- malha: a TriGMesh do main.cpp (sem gravar o VTK) e a copia so com as folhas
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

// So as folhas: o TPZUnsatAnalysis acha os nos do reservatorio pelos elementos de contorno da malha geometrica
TPZGeoMesh *Malha(const std::string &nome)
{
    UnsatModel<TUnsatHP> tmp;
    tmp.SetGeoMesh(TriGMesh(kTalude, std::stoi(nome.substr(3)))); // tmp fica dono e apaga
    std::vector<int64_t> pai;
    return tmp.RefinedGeoMesh({}, pai);
}

REAL Superficie(REAL x)
{
    const auto &p = kTalude;
    if (x <= p[4][0]) return p[4][1];
    if (x >= p[3][0]) return p[3][1];
    return p[4][1] + (p[3][1] - p[4][1]) * (x - p[4][0]) / (p[3][0] - p[4][0]);
}

// ---- material: o UnsatMaterial do main.cpp (Mohr-Coulomb associado, van Genuchten-Mualem, chi = S_w), com k_h / k_v
TUnsatHP *MaterialCeron(REAL alfa, bool acoplado)
{
    const Solo s;
    const REAL phi = s.phiDeg * M_PI / 180.;
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    auto *mat = new TUnsatHP(matSolo);
    mat->SetPlasticModel(TMC(TPZYCMohrCoulombPV2(phi, phi, s.c, ER), ER));
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
    h.fKrw = TPZUnsaturatedHydraulics::EMualem;
    h.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    mat->SetHydraulics(h);
    mat->SetPorosity(s.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(s.Kw);
    mat->SetConductivity(s.kv, s.gammaw, 1.e-3 / 1.8e-5);
    if (alfa != 1.) mat->SetAnisotropicConductivity(alfa * s.kv, s.kv); // k_x = alfa k_v, k_y = k_v
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(s.gamma);
    mat->SetGravityDirection(TPZManVector<REAL, 3>({0., -1., 0.}));
    mat->SetBishop(TUnsat::EChiSaturation, 1.);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(acoplado);
    return mat;
}

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

// ---- fator: GIM com a agua do reservatorio sobre o talude submerso (pressao gamma_w (nivel - y) normal a face e ao
// terreno do pe), que cresce com lambda como o peso; estado = {lambda, nivel}, lido pelas funcoes das tracoes
struct Fator {
    std::shared_ptr<std::array<REAL, 2>> estado;
    std::unique_ptr<UnsatModel<TUnsatHP>> m;
    std::unique_ptr<UnsatAnalysis<TUnsatHP>> a; // declarada depois: destruida antes do modelo
};

void NovoFator(Fator &f, TPZGeoMesh *g, const Opcoes &op, REAL nivel)
{
    f.a.reset();
    f.estado = std::make_shared<std::array<REAL, 2>>(std::array<REAL, 2>{1., nivel});
    f.m = std::make_unique<UnsatModel<TUnsatHP>>();
    f.m->SetGeoMesh(g);
    InsertBC(*f.m);
    const Solo s;
    auto pressao = [est = f.estado, gw = s.gammaw](REAL nx, REAL ny) { // normal para fora do solo
        return [est, gw, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
            const REAL pr = (*est)[0] * gw * std::max(REAL(0), (*est)[1] - x[1]);
            v2[0] = -pr * nx, v2[1] = -pr * ny;
        };
    };
    const auto &p = kTalude;
    const REAL dx = p[3][0] - p[4][0], dy = p[4][1] - p[3][1], L = std::hypot(dx, dy);
    f.m->AddNaturalBC(face, TUnsat::ENeumannU, {0., 0.}, pressao(dy / L, dx / L));
    f.m->AddNaturalBC(topright, TUnsat::ENeumannU, {0., 0.}, pressao(0., 1.));
    TUnsatHP *mat = MaterialCeron(op.alfa, true);
    f.m->SetMaterial(mat);
    f.m->CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
    f.a = std::make_unique<UnsatAnalysis<TUnsatHP>>(*f.m);
    f.a->SetVerbose(false);
    auto est = f.estado;
    f.a->SetLoadFactorFunction([mat, est](REAL l) { mat->SetGravityFactor(l), (*est)[0] = l; }, false);
    f.a->SetStrengthReductionFunction([mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); });
}

REAL CalculaFator(Fator &f, const Campo *campos, const PointLocator &loc, const std::vector<REAL> &pw)
{
    f.m->Material()->ResetMemory(); // estado sem tensoes (apaga tambem c, phi e k dos pontos)
    if (campos) AplicaCampos(*f.m, *campos);
    f.a->SetFrozenPressures([&](int espaco, const TPZVec<REAL> &x) {
        return espaco == 1 ? loc.Interpolate(x[0], x[1], pw) : 0.; // p_w do transiente; p_a = 0
    });
    return f.a->GIMArcLength();
}

REAL FatorAdaptativo(Fator &f, const Opcoes &op, REAL nivel, const Campo *campos, const PointLocator &loc,
                     const std::vector<REAL> &pw, std::vector<REAL> &niveis)
{
    TPZGeoMesh *g = Malha(op.malhaF);
    niveis.clear();
    for (int k = 0;; k++) {
        NovoFator(f, g, op, nivel);
        const REAL fs = CalculaFator(f, campos, loc, pw);
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

// ---- descritores (as colunas do mc_resumo.csv da campanha do rebaixamento)
struct Mecanismo {
    REAL bArea = 0., bProf = 0., bXmin = 0., bXmax = 0., bXc = 0., bYc = 0., bYmin = 0., bPw = 0., bC = 0., bPhi = 0.,
         bK = 0.;
    REAL mArea = 0., mXc = 0., mYc = 0., mC = 0., mPhi = 0., mK = 0.;
    REAL cMed = 0., phiMed = 0., kMed = 0.;
    struct Ponto {
        REAL x, y, j, du, c, phi, k;
    };
    std::vector<Ponto> pontos; // sqrt(J2) >= 1% do maior (mecanismo.csv)
    std::vector<std::array<REAL, 2>> banda;
};

Mecanismo DescreveMecanismo(Fator &f, const PointLocator &loc, const std::vector<REAL> &pw)
{
    const Solo res;
    UnsatModel<TUnsatHP> &m = *f.m;
    TUnsatHP *mat = m.Material();
    const int iu = mat->VariableIndex("Displacement");
    struct P {
        TPZCompEl *cel;
        TPZManVector<REAL, 3> qsi;
        REAL x, y, w, j, c, phi, k, du;
    };
    std::vector<P> ps;
    for (TPZCompEl *cel : m.CompMesh()->ElementVec()) {
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
    // du: o deslocamento do ultimo incremento do caminho (a forma do mecanismo); sem ele, o deslocamento final
    const bool incremento = f.a->PathIncrement().Rows() == f.a->PathState().Rows();
    if (incremento) f.a->LoadSolution(f.a->PathIncrement());
    for (P &p : ps) {
        TPZManVector<STATE, 9> v(2, 0.);
        p.cel->Solution(p.qsi, iu, v);
        p.du = std::hypot(v[0], v[1]);
    }
    if (incremento) f.a->LoadSolution(f.a->PathState());
    REAL jmax = 0., dumax = 0., area = 0.;
    Mecanismo d;
    for (const P &p : ps) {
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

// ---- gradiente de Gamma pela projecao modal (o GradienteModal do Projects/Cho2010, ramo GIM)
// No colapso do GIM, com Mohr-Coulomb associado e o teorema do envelope (as cargas externas, o peso e a agua do
// reservatorio, crescem com Gamma): dGamma/dc'(x) = dD/dc / Den e dGamma/dtan(phi')(x) = dD/dphi cos^2(phi') / Den, com
//   dD/dc = eps_v cot(phi'),  dD/dphi = -(eps_v / (2 sin(phi'))) ((s1 + s3) cos(phi') + 2 c' sin(phi')),
//   Den = P_mr / Gamma = integral de c' cot(phi') eps_v / Gamma,
// s1, s3 as tensoes efetivas principais do colapso (m_sigma: as do modelo plastico, ja com as pressoes congeladas) e
// eps_v o traco da deformacao plastica. Projecao nos modos: d/dxi_j pela regra da cadeia com ln c' e ln phi' lineares em
// xi (dZ/dxi = sqrt(lambda) f(x), ModeValues). O k nao entra na dissipacao: so nas pressoes, que aqui ficam congeladas
std::vector<double> GradienteModal(Fator &f, RandomField &campos, int K, REAL Gamma)
{
    UnsatModel<TUnsatHP> &m = *f.m;
    TUnsatHP *mat = m.Material();
    struct P {
        REAL x, y, w, gc, gt, c, phi;
    };
    std::vector<P> ps;
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
            if (!(ev > 0.)) continue; // ponto elastico: nao dissipa
            const TPZTensor<REAL> &sg = me.m_sigma;
            const REAL cen = 0.5 * (sg.XX() + sg.YY()), r = std::hypot(0.5 * (sg.XX() - sg.YY()), sg.XY());
            const REAL s1 = std::max({cen + r, cen - r, sg.ZZ()}), s3 = std::min({cen + r, cen - r, sg.ZZ()});
            const REAL c = pr[0], phi = pr[1];
            const REAL dDc = ev / std::tan(phi);
            const REAL dDp = -(ev / (2. * std::sin(phi))) * ((s1 + s3) * std::cos(phi) + 2. * c * std::sin(phi));
            P q;
            q.x = x[0], q.y = x[1], q.w = w * std::abs(detjac), q.c = c, q.phi = phi;
            q.gc = dDc, q.gt = dDp * std::cos(phi) * std::cos(phi);
            den += q.w * dDc * c / Gamma; // P_mr / Gamma
            ps.push_back(q);
        }
    }
    std::vector<double> g(3 * K, 0.);
    if (ps.empty() || !(std::abs(den) > 0.)) return g;
    const double sc = campos.SigmaLn(0), sf = campos.SigmaLn(1); // c' e phi' independentes (sem correlacao cruzada)
    for (const P &q : ps) {
        const std::vector<double> mv = campos.ModeValues(q.x, q.y, K);
        const double sec2 = 1. / (std::cos(q.phi) * std::cos(q.phi));
        const double ac = q.w * q.gc / den * q.c;          // dGamma/dln(c') do ponto
        const double af = q.w * q.gt / den * sec2 * q.phi; // dGamma/dln(phi') do ponto
        for (int k = 0; k < K; k++) g[k] += ac * sc * mv[k], g[K + k] += af * sf * mv[k];
    }
    return g; // os K de k ficam 0
}

struct Linha {
    int real = 0;
    REAL fs = 0., nivel = 0.;
    std::string situacao = "ok", niveis;
    int64_t eqF = 0, eqT = 0;
    double sF = 0., sT = 0.;
    REAL pwMin = 0., pwMax = 0.;
    Mecanismo mec;
    REAL logPeso = 0.;
    std::vector<double> proj = {0., 0., 0.}; // xi . d de c', phi' e k
};

const char *kCabecalho =
    "real,t_h,FS_min,situacao,niveis,eq_fator,s_fator,eq_transiente,s_transiente,pw_min,pw_max,nivel_reservatorio,"
    "banda_area_m2,banda_prof,banda_xmin,banda_xmax,banda_xc,banda_yc,banda_ymin,banda_pw,banda_c,banda_phi_deg,"
    "banda_k_mult,massa_area_m2,massa_xc,massa_yc,massa_c,massa_phi_deg,massa_k_mult,c_med,phi_med,k_mult_med,hw_H,"
    "nivel_final,trebaix,cov_c,cov_phi,cov_k,Lx,Ly,alfa,KL_variancia,log_peso,proj_c,proj_phi,proj_k\n";

void Grava(std::ofstream &csv, const Linha &l, const Opcoes &op, REAL nivelFinal, REAL klVariancia)
{
    const Mecanismo &d = l.mec;
    std::ostringstream s;
    s << l.real << "," << op.trebaix << "," << std::setprecision(7) << l.fs << "," << l.situacao << "," << l.niveis
      << "," << l.eqF << "," << std::setprecision(4) << l.sF << "," << l.eqT << "," << l.sT << "," << std::setprecision(6)
      << l.pwMin << "," << l.pwMax << "," << l.nivel << "," << d.bArea << "," << d.bProf << "," << d.bXmin << ","
      << d.bXmax << "," << d.bXc << "," << d.bYc << "," << d.bYmin << "," << d.bPw << "," << d.bC << "," << d.bPhi << ","
      << d.bK << "," << d.mArea << "," << d.mXc << "," << d.mYc << "," << d.mC << "," << d.mPhi << "," << d.mK << ","
      << d.cMed << "," << d.phiMed << "," << d.kMed << "," << op.hw << "," << nivelFinal << "," << op.trebaix << ","
      << op.campos.covC << "," << op.campos.covPhi << "," << op.campos.covK << "," << op.campos.Lx << ","
      << op.campos.Ly << "," << op.alfa << "," << klVariancia << "," << std::setprecision(10) << l.logPeso << ","
      << l.proj[0] << "," << l.proj[1] << "," << l.proj[2] << "\n";
    csv << s.str();
    csv.flush();
}

// mapa da banda do mecanismo: celulas de 1 m no retangulo do talude
struct Mapa {
    static constexpr int nx = 70, ny = 35;
    std::vector<int> n = std::vector<int>(nx * ny, 0), nFalha = std::vector<int>(nx * ny, 0);
    int nReal = 0;
    void Soma(const Mecanismo &d, bool falha)
    {
        std::set<int> celulas;
        for (const auto &p : d.banda)
            celulas.insert(std::clamp(int(p[1]), 0, ny - 1) * nx + std::clamp(int(p[0]), 0, nx - 1));
        for (int c : celulas) n[c]++, nFalha[c] += falha;
    }
    void Grava(const std::string &arq, REAL t) const
    {
        std::ofstream f(arq);
        f << "t_h,x,y,n,n_falha,n_real\n";
        for (int c = 0; c < nx * ny; c++)
            if (n[c]) f << t << "," << c % nx + 0.5 << "," << c / nx + 0.5 << "," << n[c] << "," << nFalha[c] << "," << nReal << "\n";
    }
};

// pontos da banda do mecanismo.csv (colunas x, y, epsp_raizJ2): os com sqrt(J2) >= 0,1 do maior, com peso 1
std::vector<std::array<double, 3>> LeBanda(const std::string &arq)
{
    std::ifstream f(arq);
    std::string linha;
    std::getline(f, linha);
    std::vector<std::string> cab;
    std::stringstream sc(linha);
    for (std::string c; std::getline(sc, c, ',');) cab.push_back(c);
    auto col = [&cab](const char *nome) { return int(std::find(cab.begin(), cab.end(), nome) - cab.begin()); };
    const int ix = col("x"), iy = col("y"), ij = col("epsp_raizJ2");
    std::vector<std::array<double, 3>> p; // x, y, sqrt(J2)
    while (std::getline(f, linha)) {
        std::vector<double> v;
        std::stringstream sl(linha);
        for (std::string c; std::getline(sl, c, ',');) v.push_back(std::stod(c));
        if (int(v.size()) > std::max({ix, iy, ij})) p.push_back({v[ix], v[iy], v[ij]});
    }
    double jmax = 0.;
    for (const auto &q : p) jmax = std::max(jmax, q[2]);
    std::vector<std::array<double, 3>> banda;
    for (const auto &q : p)
        if (q[2] >= 0.1 * jmax) banda.push_back({q[0], q[1], 1.});
    return banda;
}

// os campos de uma linha de CSV
std::vector<std::string> Divide(const std::string &linha)
{
    std::vector<std::string> v;
    std::stringstream sl(linha);
    for (std::string c; std::getline(sl, c, ',');) v.push_back(c);
    return v;
}

// O piloto da amostragem por importancia: as realizacoes 0..np-1 do Monte Carlo direto da pasta (os <prefixo>_<a>_<b>.csv
// da fila do cluster, com as colunas real, FS_min e situacao) e o conjunto elite da entropia cruzada (Rubinstein e Kroese,
// 2004): as que falharam (Gamma < 1) ou, se forem menos que a fracao rho das validas, as rho com menor Gamma (um nivel
// intermediario: Gamma <= nivel). Ordenado por (Gamma, realizacao): o mesmo em todos os blocos
struct Piloto {
    int presentes = 0, validas = 0, falhas = 0;
    double nivel = 1.;
    std::vector<int64_t> elite;
};

Piloto LePiloto(const std::string &pasta, int np, double rho)
{
    std::map<int64_t, double> g; // realizacao -> Gamma (-1: a analise nao terminou)
    DIR *d = opendir(pasta.c_str());
    if (!d) {
        std::fprintf(stderr, "piloto: nao abri a pasta %s\n", pasta.c_str());
        std::exit(1);
    }
    for (dirent *e; (e = readdir(d));) {
        const std::string nome = e->d_name;
        // <prefixo>_<a>_<b>.csv com o prefixo sem '_' (nao pega mc_mapa_*, mc_proj_*)
        if (nome.size() < 5 || nome.compare(nome.size() - 4, 4, ".csv")) continue;
        const std::string base = nome.substr(0, nome.size() - 4);
        const size_t u = base.find('_');
        if (u == std::string::npos || std::count(base.begin(), base.end(), '_') != 2 ||
            base.find_first_not_of("0123456789_", u) != std::string::npos)
            continue;
        std::ifstream f(pasta + "/" + nome);
        std::string linha;
        if (!std::getline(f, linha)) continue;
        const std::vector<std::string> cab = Divide(linha);
        auto col = [&cab](const char *c) { return int(std::find(cab.begin(), cab.end(), c) - cab.begin()); };
        const int ir = col("real"), ig = col("FS_min"), is = col("situacao");
        if (std::max({ir, ig, is}) >= int(cab.size())) continue;
        while (std::getline(f, linha)) {
            const std::vector<std::string> v = Divide(linha);
            if (int(v.size()) <= std::max({ir, ig, is})) continue; // linha incompleta
            try {
                const int64_t r = std::stoll(v[ir]);
                if (r < 0 || r >= np) continue;
                const double gamma = std::stod(v[ig]);
                g[r] = v[is] == "ok" && gamma > 0. ? gamma : -1.;
            } catch (const std::exception &) {
            }
        }
    }
    closedir(d);
    Piloto p;
    p.presentes = g.size();
    std::vector<std::pair<double, int64_t>> v;
    for (const auto &x : g)
        if (x.second > 0.) v.push_back({x.second, x.first});
    std::sort(v.begin(), v.end());
    p.validas = v.size();
    p.falhas = std::count_if(v.begin(), v.end(), [](const std::pair<double, int64_t> &x) { return x.first < 1.; });
    const int m = std::min(p.validas, std::max(1, int(std::ceil(rho * p.validas))));
    const int ne = p.falhas >= m ? p.falhas : m;
    for (int i = 0; i < ne; i++) p.elite.push_back(v[i].second);
    if (p.falhas < m && m > 0) p.nivel = v[m - 1].first;
    return p;
}

int main(int argc, char *argv[])
{
    const Opcoes op = LeOpcoes(argc, argv);
    const Solo res;
    const Campos &cp = op.campos;
    const auto &p = kTalude;
    const REAL H = p[4][1] - p[3][1], nivelFinal = p[4][1] - op.hw * H;
    RandomField campos(0., 0., p[1][0], p[4][1], cp.Lx, cp.Ly, cp.M,
                       {{res.c, cp.covC}, {res.phiDeg, cp.covPhi}, {1., cp.covK}}); // c' (kPa), phi' (graus), k / k_v
    const Campo campoAleatorio = [&campos](double x, double y) { return campos.At(x, y); };
    const Campo campoFixo = [v = op.valores](double, double) { return std::vector<double>(v.begin(), v.end()); };
    // amostragem por importancia: a direcao da banda e o deslocamento
    std::vector<double> dir;
    if (!op.banda.empty()) {
        dir = campos.Direction(LeBanda(op.banda));
        if (dir.empty()) {
            std::fprintf(stderr, "banda sem pontos: %s\n", op.banda.c_str());
            return 1;
        }
        if (!op.desloc.empty()) campos.SetShift(dir, std::vector<double>(op.desloc.begin(), op.desloc.end()));
        std::printf("  amostragem por importancia: banda de %s, deslocamento %s\n", op.banda.c_str(),
                    op.desloc.empty() ? "nenhum" : (std::to_string(op.desloc[0]) + ", " + std::to_string(op.desloc[1]) +
                                                   ", " + std::to_string(op.desloc[2])).c_str());
    }
    std::printf("Monte Carlo do rebaixamento do Ceron com o transiente: h_w/H = %g (nivel %g -> %g m em %g h), Gamma pelo GIM "
                "| realizacoes %d..%d%s\n  campos: CoV c' %g, phi' %g, k %g; Lx = %g m, Ly = %g m; alfa = %g; KL: %d "
                "termos, %.1f%% da variancia\n  malhas: transiente %s, fator %s + %d niveis\n",
                op.hw, p[4][1], nivelFinal, op.trebaix, op.primeira, op.primeira + op.n - 1,
                op.det ? " e a deterministica" : "", cp.covC, cp.covPhi, cp.covK, cp.Lx, cp.Ly, op.alfa, campos.M(),
                100. * campos.VarianceFraction(), op.malhaT.c_str(), op.malhaF.c_str(), op.nref);

    const std::string faixa = std::to_string(op.primeira) + "_" + std::to_string(op.primeira + op.n - 1);
    const std::string arq = op.saida + "_" + faixa + ".csv";
    std::ofstream csv(arq);
    csv << kCabecalho;
    if (op.primeira == 0) { // so o bloco que comeca na realizacao 0 (no cluster, varios blocos rodam ao mesmo tempo)
        std::ofstream pa(op.saida + "_parametros.txt");
        pa << "HW_H=" << op.hw << "\nNIVEL_FINAL=" << nivelFinal << "\nTREBAIX=" << op.trebaix << "\nKV=" << res.kv
           << "\nCOV_C=" << cp.covC << "\nCOV_PHI=" << cp.covPhi << "\nCOV_K=" << cp.covK << "\nLX=" << cp.Lx
           << "\nLY=" << cp.Ly << "\nALFA=" << op.alfa << "\nKL_M=" << campos.M()
           << "\nKL_VARIANCIA=" << campos.VarianceFraction() << "\nSEMENTE=" << cp.semente
           << "\nMALHA_TRANSIENTE=" << op.malhaT << "\nMALHA_FATOR=" << op.malhaF << "\nNREF=" << op.nref
           << "\nFATOR=gim_arco\nAI="
           << (!op.piloto.empty()      ? op.ai + "_piloto"
               : !op.desloc.empty()    ? std::string("banda_fixo")
               : !op.deslocModos.empty() ? std::string("modos_fixo")
                                         : std::string("nenhuma"))
           << "\n";
        if (!op.piloto.empty())
            pa << "PILOTO=" << op.piloto << "\nPILOTO_N=" << op.pilotoN << "\nPILOTO_SEMENTE=" << op.pilotoSemente
               << "\nELITE=" << op.elite << "\nMODOS=" << op.modos << "\n";
    }
    Mapa mapa;
    std::vector<REAL> gammas;
    double somaT = 0., somaF = 0.;

    if (!op.mediaFalhas.empty()) { // so a media dos K primeiros modos de cada campo nas realizacoes listadas
        std::ifstream lf(op.mediaFalhas);
        std::vector<std::vector<double>> soma(3, std::vector<double>(op.modos, 0.));
        int nf = 0;
        for (int64_t r; lf >> r; nf++) {
            campos.Sample(cp.semente, r);
            for (int c = 0; c < 3; c++)
                for (int k = 0; k < op.modos && k < campos.M(); k++) soma[c][k] += campos.Xi()[c][k];
        }
        if (!nf) {
            std::fprintf(stderr, "nenhuma realizacao em %s\n", op.mediaFalhas.c_str());
            return 1;
        }
        std::ofstream df(op.saida + "_desloc_modos.csv");
        df << "termo,c,phi,k\n" << std::setprecision(10);
        double n2 = 0.;
        for (int k = 0; k < op.modos && k < campos.M(); k++) {
            df << k;
            for (int c = 0; c < 3; c++) df << "," << soma[c][k] / nf, n2 += std::pow(soma[c][k] / nf, 2);
            df << "\n";
        }
        std::printf("  media de %d realizacoes nos %d primeiros modos: |mu|^2 = %.3f; %s\n", nf, op.modos, n2,
                    (op.saida + "_desloc_modos.csv").c_str());
        return 0;
    }
    if (!op.deslocModos.empty()) { // deslocamento por modo, do arquivo
        std::ifstream df(op.deslocModos);
        std::string linha;
        std::getline(df, linha);
        std::vector<std::vector<double>> mu(3);
        while (std::getline(df, linha)) {
            std::stringstream sl(linha);
            std::string v;
            std::getline(sl, v, ','); // termo
            for (int c = 0; c < 3 && std::getline(sl, v, ','); c++) mu[c].push_back(std::stod(v));
        }
        campos.SetShift(mu);
        double n2 = 0.;
        for (const auto &m : mu)
            for (double v : m) n2 += v * v;
        std::printf("  amostragem por importancia: deslocamento por modo de %s (%zu termos, |mu|^2 = %.3f)\n",
                    op.deslocModos.c_str(), mu[0].size(), n2);
    }
    if (!op.piloto.empty()) { // deslocamento da entropia cruzada no conjunto elite do piloto
        const Piloto pl = LePiloto(op.piloto, op.pilotoN, op.elite);
        if (pl.presentes < 0.9 * op.pilotoN || pl.elite.empty()) {
            std::fprintf(stderr, "piloto incompleto em %s: %d das realizacoes 0..%d, %d validas\n", op.piloto.c_str(),
                         pl.presentes, op.pilotoN - 1, pl.validas);
            return 3;
        }
        const double ne = pl.elite.size();
        double mu2 = 0.;
        std::ostringstream reg; // o registro do deslocamento
        reg << std::setprecision(10);
        if (op.ai == "banda") { // s_c = media de xi_c . d
            std::vector<double> s(3, 0.);
            for (int64_t r : pl.elite) {
                campos.Sample(op.pilotoSemente, r); // sem deslocamento ainda: os xi do piloto
                const std::vector<double> pr = campos.Projection(dir);
                for (int c = 0; c < 3; c++) s[c] += pr[c] / ne;
            }
            campos.SetShift(dir, s);
            for (double x : s) mu2 += x * x;
            reg << "s=" << s[0] << "," << s[1] << "," << s[2] << "\nmu2=" << mu2 << "\n";
        } else { // a media dos K primeiros modos de cada campo
            const int K = std::min(op.modos, campos.M());
            std::vector<std::vector<double>> mu(3, std::vector<double>(K, 0.));
            for (int64_t r : pl.elite) {
                campos.Sample(op.pilotoSemente, r);
                for (int c = 0; c < 3; c++)
                    for (int k = 0; k < K; k++) mu[c][k] += campos.Xi()[c][k] / ne;
            }
            campos.SetShift(mu);
            for (const auto &m : mu)
                for (double x : m) mu2 += x * x;
            reg << "modos=" << K << "\nmu2=" << mu2 << "\ntermo,c,phi,k\n";
            for (int k = 0; k < K; k++) reg << k << "," << mu[0][k] << "," << mu[1][k] << "," << mu[2][k] << "\n";
        }
        std::ofstream df(op.saida + "_desloc_" + faixa + ".txt");
        df << "# deslocamento da amostragem por importancia: a media dos xi do conjunto elite do piloto\nai=" << op.ai
           << "\npiloto=" << op.piloto << "\npiloto_n=" << op.pilotoN << "\npiloto_semente=" << op.pilotoSemente
           << "\npresentes=" << pl.presentes << "\nvalidas=" << pl.validas << "\nfalhas=" << pl.falhas
           << "\nelite=" << pl.elite.size() << "\nnivel=" << pl.nivel << "\n"
           << reg.str();
        std::printf("  amostragem por importancia (%s) pelo piloto %s: %d validas de %d, %d falhas, elite de %zu "
                    "(Gamma < %.4f); |mu|^2 = %.4f\n",
                    op.ai.c_str(), op.piloto.c_str(), pl.validas, op.pilotoN, pl.falhas, pl.elite.size(), pl.nivel, mu2);
    }
    if (op.soProjecoes) { // so os campos: as projecoes e o peso de cada realizacao
        std::ofstream pj(op.saida + "_proj_" + faixa + ".csv");
        pj << "real,log_peso,proj_c,proj_phi,proj_k\n" << std::setprecision(10);
        for (int n = op.primeira; n < op.primeira + op.n; n++) {
            campos.Sample(cp.semente, n);
            const std::vector<double> pr = campos.Projection(dir);
            pj << n << "," << campos.LogWeight() << "," << pr[0] << "," << pr[1] << "," << pr[2] << "\n";
        }
        std::printf("  projecoes em %s\n", (op.saida + "_proj_" + faixa + ".csv").c_str());
        std::ofstream dd(op.saida + "_direcao.csv"); // a direcao unitaria d da banda, por termo da KL
        dd << "termo,d\n" << std::setprecision(10);
        for (size_t k = 0; k < dir.size(); k++) dd << k << "," << dir[k] << "\n";
        return 0;
    }
    std::vector<std::vector<double>> xiDados; // opcao xi=: os K de c', de phi' e de k
    if (!op.xi.empty()) {
        std::ifstream fx(op.xi);
        std::vector<double> v;
        for (double x; fx >> x;) v.push_back(x);
        if (v.empty() || v.size() % 3 || int(v.size() / 3) > campos.M()) {
            std::fprintf(stderr, "xi: %s tem %zu numeros; precisa de 3K, com K <= M = %d\n", op.xi.c_str(), v.size(), campos.M());
            return 1;
        }
        const size_t K = v.size() / 3;
        for (int c = 0; c < 3; c++) xiDados.push_back(std::vector<double>(v.begin() + c * K, v.begin() + (c + 1) * K));
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
        const Campo *cf = deterministica ? (op.valores.empty() ? nullptr : &campoFixo) : &campoAleatorio;
        Linha l;
        l.real = n;
        if (!deterministica) l.logPeso = campos.LogWeight(), l.proj = campos.Projection(dir);

        // 1. transiente do rebaixamento: so a agua, com o k de cada ponto
        UnsatModel<TUnsatHP> mT;
        mT.SetGeoMesh(Malha(op.malhaT));
        InsertBC(mT);
        mT.SetMaterial(MaterialCeron(op.alfa, false));
        mT.CreateCompMesh({{2, 2}, {1, 1}, {1, 1}});
        if (cf) AplicaCampos(mT, *cf);
        TPZUnsatAnalysis aT(mT.CompMesh(), {{base, {true, true}}, {right, {true, false}}, {left, {true, false}}});
        TPZSpStructMatrix<STATE> strmat(mT.CompMesh());
        strmat.SetNumThreads(std::getenv("PZ_THREADS") ? std::atoi(std::getenv("PZ_THREADS")) : 8);
        aT.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        aT.SetSolver(solver);
        DrawdownData d;
        d.gammaw = res.gammaw;
        d.nivelInicial = p[4][1]; // reservatorio cheio: na crista
        d.nivelFinal = nivelFinal;
        d.nivelMacico = p[4][1];  // NA do macico (lateral esquerda), fixo
        d.tRebaixamento = op.trebaix;
        DrawdownTransient rebaixamento(aT, d, left, {right, topright, face, topleft});
        rebaixamento.InitialState();
        const PointLocator loc(mT.GeoMesh());
        const auto tt = Clock::now();
        const bool ok = rebaixamento.AdvanceTo(op.trebaix);
        l.sT = Seg(tt), l.eqT = mT.CompMesh()->NEquations(), l.nivel = rebaixamento.Level();
        const std::vector<REAL> pw = rebaixamento.NodalWaterPressure();
        l.pwMin = *std::min_element(pw.begin(), pw.end()), l.pwMax = *std::max_element(pw.begin(), pw.end());

        // 2. Gamma no fim do rebaixamento
        if (!ok) {
            l.situacao = "transiente";
            std::printf("  realizacao %d: o transiente parou em %.3f h\n", n, rebaixamento.Time());
        } else {
            const bool figura = op.figuras.count(n) > 0;
            std::unique_ptr<AnimationFiles<TUnsatHP>> anim;
            if (figura) {
                const std::string pasta = op.saida + "_real" + std::to_string(n);
                anim = std::make_unique<AnimationFiles<TUnsatHP>>(pasta, mT);
                anim->Summary(op.trebaix, d.nivelInicial, d.nivelFinal);
                anim->Snapshot(op.trebaix, pw, l.nivel, mT.Material()->Hydraulics());
                std::ofstream cfs(pasta + "/campos.csv");
                cfs << "no,x,y,c,phi,k\n";
                TPZGeoMesh *g = mT.GeoMesh();
                for (int64_t i = 0; i < g->NNodes(); i++) {
                    const REAL x = g->NodeVec()[i].Coord(0), y = g->NodeVec()[i].Coord(1);
                    const std::vector<double> v = cf ? (*cf)(x, y) : std::vector<double>{res.c, res.phiDeg, 1.};
                    cfs << i << "," << x << "," << y << "," << v[0] << "," << v[1] << "," << v[2] << "\n";
                }
            }
            const auto tf = Clock::now();
            Fator f;
            std::vector<REAL> niveis;
            l.fs = FatorAdaptativo(f, op, l.nivel, cf, loc, pw, niveis);
            l.sF = Seg(tf), l.eqF = f.m->CompMesh()->NEquations();
            if (!(l.fs > 0.)) l.situacao = "fator";
            for (size_t i = 0; i < niveis.size(); i++) l.niveis += (i ? ";" : "") + std::to_string(niveis[i]).substr(0, 7);
            l.mec = DescreveMecanismo(f, loc, pw);
            if (op.gradiente > 0 && !deterministica && l.fs > 0.) {
                const std::vector<double> gr = GradienteModal(f, campos, std::min(op.gradiente, campos.M()), l.fs);
                std::ofstream gf(op.saida + "_grad_" + std::to_string(n) + ".txt");
                gf << std::setprecision(10)
                   << "# Gamma (GIM) e o gradiente dGamma/dxi (os K de c', de phi' e de k), projecao modal\n" << l.fs << "\n";
                for (double x : gr) gf << x << "\n";
            }
            if (anim) {
                anim->FactorOfSafety(op.trebaix, l.fs, *f.m, *f.a);
                std::ofstream mf(op.saida + "_real" + std::to_string(n) + "/mecanismo.csv");
                mf << "t_h,x,y,epsp_raizJ2,du,c,phi_deg,k_mult\n";
                for (const auto &q : l.mec.pontos)
                    mf << op.trebaix << "," << q.x << "," << q.y << "," << q.j << "," << q.du << "," << q.c << ","
                       << q.phi << "," << q.k << "\n";
            }
        }
        Grava(csv, l, op, nivelFinal, campos.VarianceFraction());
        if (!deterministica && l.situacao == "ok") {
            gammas.push_back(l.fs);
            mapa.Soma(l.mec, l.fs < 1.);
            mapa.nReal++, somaT += l.sT, somaF += l.sF;
        }
        std::printf("  realizacao %4d: Gamma = %.4f (%s) | %.1f s (transiente %.1f s)\n", n, l.fs, l.situacao.c_str(),
                    Seg(t0), l.sT);
        std::fflush(stdout);
    }
    mapa.Grava(op.saida + "_mapa_" + faixa + ".csv", op.trebaix);

    if (!gammas.empty()) {
        REAL m = 0., v = 0.;
        for (REAL x : gammas) m += x;
        m /= gammas.size();
        for (REAL x : gammas) v += (x - m) * (x - m);
        const int falhas = std::count_if(gammas.begin(), gammas.end(), [](REAL x) { return x < 1.; });
        std::printf("\n  %zu realizacoes: Gamma medio %.4f, desvio %.4f, menor %.4f, P(Gamma < 1) = %.2f%%\n",
                    gammas.size(), m, gammas.size() > 1 ? std::sqrt(v / (gammas.size() - 1)) : 0.,
                    *std::min_element(gammas.begin(), gammas.end()), 100. * falhas / gammas.size());
        std::printf("  tempo medio por realizacao: transiente %.1f s, fator %.1f s\n", somaT / mapa.nReal,
                    somaF / mapa.nReal);
    }
    std::printf("  arquivos: %s e %s\n", arq.c_str(), (op.saida + "_mapa_" + faixa + ".csv").c_str());
    return 0;
}
