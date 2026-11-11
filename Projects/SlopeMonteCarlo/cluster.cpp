// cluster.cpp (alvo SlopeMonteCarloCluster): as Secoes 6.1 a 6.5 de Vargas Ceron et al. (2025) no cluster (fermi), com o
// mesmo modelo do main.cpp deste projeto: malha ceron_h5_a45_g15, refinamento adaptativo em h (nunif = 0, nref = 3,
// pmax = 2), percolacao permanente do rebaixamento e Gamma pelo GIM com arc-length no material nao saturado, campos
// lognormais de c', phi' e k_v por Karhunen-Loeve (M = 2000). Ruptura: Gamma < 1.
//
// Os casos sao os das Tabelas 3 a 6 do artigo, e as pastas e arquivos sao os do exemplo do fermi (campanha ceron_sat do
// SlopeRandomField), para que SlopeRandomField/comparacao_ceron/dados_sat.py e relatorio_sat.py os leiam sem mudancas:
//   <raiz>/<prefixo>_fase1/f50_I0_iii                  referencia: Fig. 25 e Spearman (6.1)
//   <raiz>/<prefixo>_fase2/f50_I0_iii_covk0, ...        CoV de k_v, c' e phi' (6.2), s (6.3) e alfa = k_h/k_v (6.4)
//   <raiz>/<prefixo>_fase3/f50_I0_iii_nivelfinal7.5,... h_w/H de 0,5 a 0,9 (6.5)
// O nome f50_I0_iii vem do programa antigo (solo f50, sem chuva, variante iii) e aqui so identifica as pastas. Em cada
// caso: mc_resumo.csv (uma linha por realizacao: a juncao de partes/p*.csv) e parametros.txt; nos casos com
// deterministico proprio (referencia, alfa e h_w): <caso>_det/r00000/ com resumo.csv, mecanismo.csv e resultado.txt.
//
// Uso: SlopeMonteCarloCluster casos [casos = todos]
//      SlopeMonteCarloCluster rodar <raiz> <primeira> <ultima> <k> <K> [casos = todos]
//      SlopeMonteCarloCluster junta <raiz>
//      SlopeMonteCarloCluster resumo <raiz>
// rodar: o processo k de K faz as tarefas i com i % K == k de uma fila com os deterministicos primeiro e depois as
// realizacoes primeira..ultima dos casos, intercaladas (a 0 de todos, depois a 1, ...). Pula as que ja estao em algum
// partes/*.csv do caso: rodar de novo continua, e com <ultima> maior acrescenta realizacoes. Com o arquivo <raiz>/PARAR,
// nenhuma tarefa nova comeca. casos: ids separados por virgulas (ref, covk=0.75, covc=0.5, covphi=0.2, s=2, alfa=3,
// hw=0.5), fases (fase1, fase2, fase3) ou todos.
// junta: o mc_resumo.csv de cada caso. resumo: junta e compara cada caso com o artigo (o Spearman da 6.1 e as figuras
// ficam com o dados_sat.py e o relatorio_sat.py).
// Ambiente: SMC_PREFIXO (padrao ceron_mc) e SMC_VERSAO (vai para a coluna versao).
#include "SlopeAnalysis.h"
#include "SlopeFields.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <unistd.h>

using namespace slope;
namespace fs = std::filesystem;

// ---- o modelo da ultima simulacao (main.cpp)
const char *kMalha = "ceron_h5_a45_g15.msh";
const int kNunif = 0, kNref = 3, kPmax = 2;
const uint64_t kSemente = 20250416; // a realizacao n tem os mesmos numeros aleatorios em todos os casos
const std::string kBase = "f50_I0_iii";

// ---- um caso: o que muda em relacao a referencia (Tabela 2) e os valores do artigo
struct Caso {
    std::string id;     // como no dados_sat.py: ref, covk=0.75, s=2, ...
    int fase;
    std::string sufixo; // pasta kBase + sufixo
    REAL covK = 0.6, covC = 0.3, covPhi = 0.1, s = 1., alfa = 1., hwH = 1.;
    bool detProprio = false; // o deterministico muda com o caso
    REAL mu, sigma, cov, pf, covPf, det; // artigo: lognormal de Gamma, P_f (%) e Gamma deterministico (0: o da referencia)
};

std::vector<Caso> Casos()
{
    std::vector<Caso> v;
    auto add = [&v](Caso c) { v.push_back(c); };
    //    id            fase sufixo             covK  covC  covPhi s     alfa hwH  det    mu     sigma  cov   pf     covPf  det
    add({"ref",          1, "",                 0.6,  0.3,  0.1,   1.,   1.,  1.,  true,  1.353, 0.318, 23.5, 11.50, 2.77, 1.336});
    add({"covk=0",       2, "_covk0",           0.0,  0.3,  0.1,   1.,   1.,  1.,  false, 1.329, 0.288, 21.7, 11.21, 2.81, 0.});
    add({"covk=0.75",    2, "_covk0.75",        0.75, 0.3,  0.1,   1.,   1.,  1.,  false, 1.375, 0.341, 24.8, 11.41, 2.79, 0.});
    add({"covk=0.9",     2, "_covk0.9",         0.9,  0.3,  0.1,   1.,   1.,  1.,  false, 1.378, 0.354, 25.7, 11.63, 2.76, 0.});
    add({"covk=1",       2, "_covk1",           1.0,  0.3,  0.1,   1.,   1.,  1.,  false, 1.391, 0.360, 25.9, 11.44, 2.78, 0.});
    add({"covc=0.1",     2, "_covc0.1",         0.6,  0.1,  0.1,   1.,   1.,  1.,  false, 1.367, 0.180, 13.2, 0.45,  4.97, 0.});
    add({"covc=0.5",     2, "_covc0.5",         0.6,  0.5,  0.1,   1.,   1.,  1.,  false, 1.324, 0.472, 35.7, 25.78, 1.70, 0.});
    add({"covc=0.7",     2, "_covc0.7",         0.6,  0.7,  0.1,   1.,   1.,  1.,  false, 1.296, 0.620, 47.9, 36.34, 1.32, 0.});
    add({"covphi=0.05",  2, "_covphi0.05",      0.6,  0.3,  0.05,  1.,   1.,  1.,  false, 1.351, 0.300, 22.2, 10.20, 2.97, 0.});
    add({"covphi=0.15",  2, "_covphi0.15",      0.6,  0.3,  0.15,  1.,   1.,  1.,  false, 1.357, 0.346, 25.5, 12.71, 2.62, 0.});
    add({"covphi=0.2",   2, "_covphi0.2",       0.6,  0.3,  0.2,   1.,   1.,  1.,  false, 1.361, 0.384, 28.2, 15.18, 2.36, 0.});
    add({"s=1.5",        2, "_lx30_ly3",        0.6,  0.3,  0.1,   1.5,  1.,  1.,  false, 1.356, 0.352, 25.9, 13.94, 2.48, 0.});
    add({"s=2",          2, "_lx40_ly4",        0.6,  0.3,  0.1,   2.,   1.,  1.,  false, 1.361, 0.372, 27.3, 14.99, 2.38, 0.});
    add({"s=5",          2, "_lx100_ly10",      0.6,  0.3,  0.1,   5.,   1.,  1.,  false, 1.350, 0.404, 29.9, 18.93, 2.07, 0.});
    add({"s=10",         2, "_lx200_ly20",      0.6,  0.3,  0.1,   10.,  1.,  1.,  false, 1.348, 0.419, 31.1, 19.85, 2.01, 0.});
    add({"s=20",         2, "_lx400_ly40",      0.6,  0.3,  0.1,   20.,  1.,  1.,  false, 1.355, 0.430, 31.7, 20.23, 1.99, 0.});
    add({"s=400",        2, "_lx8000_ly800",    0.6,  0.3,  0.1,   400., 1.,  1.,  false, 1.348, 0.432, 32.1, 21.93, 1.89, 0.});
    add({"alfa=2",       2, "_alfa2",           0.6,  0.3,  0.1,   1.,   2.,  1.,  true,  1.562, 0.390, 24.9, 3.90,  3.51, 1.533});
    add({"alfa=3",       2, "_alfa3",           0.6,  0.3,  0.1,   1.,   3.,  1.,  true,  1.704, 0.432, 25.3, 1.83,  4.23, 1.674});
    add({"alfa=4",       2, "_alfa4",           0.6,  0.3,  0.1,   1.,   4.,  1.,  true,  1.821, 0.474, 26.0, 1.00,  4.98, 1.783});
    add({"alfa=5",       2, "_alfa5",           0.6,  0.3,  0.1,   1.,   5.,  1.,  true,  1.910, 0.504, 26.4, 0.63,  4.74, 1.872});
    add({"hw=0.5",       3, "_nivelfinal7.5",   0.6,  0.3,  0.1,   1.,   1.,  0.5, true,  1.670, 0.378, 22.6, 1.38,  4.88, 1.671});
    add({"hw=0.6",       3, "_nivelfinal7",     0.6,  0.3,  0.1,   1.,   1.,  0.6, true,  1.502, 0.345, 22.9, 4.20,  4.78, 1.494});
    add({"hw=0.7",       3, "_nivelfinal6.5",   0.6,  0.3,  0.1,   1.,   1.,  0.7, true,  1.393, 0.314, 22.5, 8.54,  3.27, 1.383});
    add({"hw=0.8",       3, "_nivelfinal6",     0.6,  0.3,  0.1,   1.,   1.,  0.8, true,  1.338, 0.309, 23.1, 11.99, 2.71, 1.322});
    add({"hw=0.9",       3, "_nivelfinal5.5",   0.6,  0.3,  0.1,   1.,   1.,  0.9, true,  1.325, 0.306, 23.1, 12.58, 2.64, 1.307});
    return v;
}

std::string Prefixo()
{
    const char *p = std::getenv("SMC_PREFIXO");
    return p && *p ? p : "ceron_mc";
}
std::string Pasta(const std::string &raiz, const Caso &c)
{
    return raiz + "/" + Prefixo() + "_fase" + std::to_string(c.fase) + "/" + kBase + c.sufixo;
}
std::string PastaDet(const std::string &raiz, const Caso &c) { return Pasta(raiz, c) + "_det/r00000"; }

// Dados do talude do caso (Data padrao = o talude do artigo) e a estatistica dos campos
Data DadosCaso(const Caso &c)
{
    Data d;
    d.alpha = c.alfa;
    d.hw = c.hwH * d.H;
    return d;
}
FieldStats CamposCaso(const Caso &c)
{
    FieldStats s;
    s.covC = c.covC, s.covPhi = c.covPhi, s.covK = c.covK;
    s.Lx = 20. * c.s, s.Ly = 2. * c.s;
    return s;
}

std::string Num(REAL v)
{
    if (!std::isfinite(v)) return "";
    std::ostringstream o;
    o << std::setprecision(10) << v;
    return o.str();
}

// ---- grandezas de uma realizacao a partir dos pontos do ultimo nivel (como o programa antigo): estatistica dos campos
// ponderada pela area; banda de ruptura (sqrt(J2(eps_p)) >= 0,1 do maximo) e massa que desliza (|du| >= 0,5 do maximo,
// du o ultimo incremento do caminho do GIM) com area, centroide e propriedades medias
struct Grupo {
    REAL area = 0., xc = 0., yc = 0., ymin = 1.e300, c = 0., phi = 0., k = 0.;
};
Grupo Media(const std::vector<PointState> &pts, const std::function<bool(const PointState &)> &dentro)
{
    Grupo g;
    for (const auto &p : pts) {
        if (!dentro(p)) continue;
        g.area += p.area, g.xc += p.area * p.x, g.yc += p.area * p.y, g.ymin = std::min(g.ymin, p.y);
        g.c += p.area * p.c, g.phi += p.area * p.phi * 180. / M_PI, g.k += p.area * p.kmult;
    }
    if (g.area > 0.) g.xc /= g.area, g.yc /= g.area, g.c /= g.area, g.phi /= g.area, g.k /= g.area;
    else g = Grupo{0., NAN, NAN, NAN, NAN, NAN, NAN};
    return g;
}

using Linha = std::vector<std::pair<std::string, std::string>>;

Linha Resultado(const Caso &c, int64_t real, bool campos, const RandomFields *rf, const std::vector<LevelSummary> &niveis,
                double seg)
{
    const LevelSummary &fim = niveis.back();
    const auto &pts = fim.points;
    Linha r = {{"real", std::to_string(real)}, {"campos", campos ? "1" : "0"}, {"semente", std::to_string(kSemente)},
               {"KL_M", campos ? std::to_string(rf->M()) : "0"}, {"KL_variancia", campos ? Num(rf->VarianceFraction()) : ""},
               {"cov_c", Num(c.covC)}, {"cov_phi", Num(c.covPhi)}, {"cov_k", Num(c.covK)}, {"Lx", Num(20. * c.s)},
               {"Ly", Num(2. * c.s)}, {"alfa", Num(c.alfa)}, {"hw_H", Num(c.hwH)},
               {"nivel_final", Num(DadosCaso(c).Nivel())}};
    // campos no solo: media e desvio ponderados pela area, minimo e maximo
    REAL area = 0.;
    for (const auto &p : pts) area += p.area;
    r.push_back({"area_solo_m2", Num(area)});
    auto estat = [&](const std::string &nome, const std::function<REAL(const PointState &)> &f) {
        REAL m = 0., v = 0., lo = 1.e300, hi = -1.e300;
        for (const auto &p : pts) m += p.area * f(p) / area;
        for (const auto &p : pts) v += p.area * (f(p) - m) * (f(p) - m) / area, lo = std::min(lo, f(p)), hi = std::max(hi, f(p));
        r.push_back({nome + "_med", Num(m)});
        r.push_back({nome + "_dp", Num(std::sqrt(v))});
        r.push_back({nome + "_min", Num(lo)});
        r.push_back({nome + "_max", Num(hi)});
    };
    estat("c", [](const PointState &p) { return p.c; });
    estat("phi", [](const PointState &p) { return p.phi * 180. / M_PI; });
    estat("k_mult", [](const PointState &p) { return p.kmult; });
    r.push_back({"FS_min", Num(fim.value)}); // o Gamma da realizacao (o nome do programa antigo, lido pelo dados_sat.py)
    for (const auto &l : niveis) r.push_back({"gamma_nivel" + std::to_string(l.level), Num(l.value)});
    REAL epsMax = 0., duMax = 0.;
    for (const auto &p : pts) epsMax = std::max(epsMax, p.sqrtJ2), duMax = std::max(duMax, p.du);
    const Grupo b = Media(pts, [&](const PointState &p) { return epsMax > 0. && p.sqrtJ2 >= 0.1 * epsMax; });
    const Grupo m = Media(pts, [&](const PointState &p) { return duMax > 0. && p.du >= 0.5 * duMax; });
    for (const auto &[nome, g] : {std::pair<std::string, Grupo>{"banda", b}, {"massa", m}}) {
        r.push_back({nome + "_area_m2", Num(g.area)});
        r.push_back({nome + "_xc", Num(g.xc)});
        r.push_back({nome + "_yc", Num(g.yc)});
        if (nome == "banda") r.push_back({"banda_ymin", Num(g.ymin)});
        r.push_back({nome + "_c", Num(g.c)});
        r.push_back({nome + "_phi_deg", Num(g.phi)});
        r.push_back({nome + "_k_mult", Num(g.k)});
    }
    std::string versao = std::getenv("SMC_VERSAO") ? std::getenv("SMC_VERSAO") : "-";
    std::replace(versao.begin(), versao.end(), ',', ';');
    std::replace(versao.begin(), versao.end(), ' ', '_');
    r.insert(r.end(), {{"triangulos", std::to_string(fim.ntri)}, {"equacoes", std::to_string(fim.neq)},
                       {"passos", std::to_string(fim.steps)}, {"fim", fim.end},
                       {"situacao", fim.value > 0. ? "ok" : "sem_equilibrio"}, {"seg", Num(seg)}, {"versao", versao}});
    return r;
}

std::string Cabecalho(const Linha &l)
{
    std::string s;
    for (size_t i = 0; i < l.size(); i++) s += (i ? "," : "") + l[i].first;
    return s;
}
std::string Valores(const Linha &l)
{
    std::string s;
    for (size_t i = 0; i < l.size(); i++) s += (i ? "," : "") + l[i].second;
    return s;
}
void GravaAtomico(const std::string &arq, const std::string &conteudo) // tmp + rename: nunca fica pela metade
{
    char host[256] = "";
    ::gethostname(host, sizeof host - 1);
    const std::string tmp = arq + ".tmp." + host + "." + std::to_string(::getpid()); // nos diferentes podem ter o mesmo pid
    std::ofstream(tmp) << conteudo;
    fs::rename(tmp, arq);
}

void GravaParametros(const std::string &pasta, const Caso &c)
{
    fs::create_directories(pasta);
    std::ostringstream o;
    o << "MESH=" << kMalha << "\nNUNIF=" << kNunif << "\nNREF=" << kNref << "\nPMAX=" << kPmax
      << "\nFATOR=gim_arco\nAGUA=percolacao_permanente\nMATERIAL=nao_saturado\nCOES=10\nPHI_DEG=30\nGAMMA=20\nKV=1e-6\n"
      << "H=5\nBETA_DEG=45\nXCRISTA=10\nYCRISTA=10\nXTOE=15\nYTOE=5\nCOV_C=" << c.covC << "\nCOV_PHI=" << c.covPhi
      << "\nCOV_K=" << c.covK << "\nLX=" << 20. * c.s << "\nLY=" << 2. * c.s << "\nALFA=" << c.alfa << "\nHW_H=" << c.hwH
      << "\nNIVEL_FINAL=" << DadosCaso(c).Nivel() << "\nKL_M=" << FieldStats().M << "\nSEMENTE=" << kSemente << "\n";
    GravaAtomico(pasta + "/parametros.txt", o.str());
}

// As linhas completas de partes/*.csv de um caso, por realizacao (a primeira de cada uma), e o cabecalho. Uma linha com
// outro numero de colunas (cortada por um processo morto no meio da escrita) nao conta
std::map<int64_t, std::string> LinhasPartes(const std::string &pasta, std::string &cab)
{
    std::map<int64_t, std::string> linhas;
    if (!fs::is_directory(pasta + "/partes")) return linhas;
    for (const auto &e : fs::directory_iterator(pasta + "/partes")) {
        if (e.path().extension() != ".csv") continue;
        std::ifstream in(e.path());
        std::string c, l;
        if (!std::getline(in, c)) continue;
        if (cab.empty()) cab = c;
        const auto virgulas = std::count(c.begin(), c.end(), ',');
        while (std::getline(in, l))
            if (!l.empty() && std::isdigit((unsigned char)l[0]) && std::count(l.begin(), l.end(), ',') == virgulas)
                linhas.emplace(std::stoll(l.substr(0, l.find(','))), l);
    }
    return linhas;
}

// Realizacoes ja feitas de um caso
std::set<int64_t> Feitas(const std::string &pasta)
{
    std::string cab;
    std::set<int64_t> f;
    for (const auto &[r, l] : LinhasPartes(pasta, cab)) f.insert(r);
    return f;
}

// Uma tarefa: o deterministico (real < 0) ou uma realizacao do caso; devolve a linha do resumo e os pontos do ultimo nivel
struct Tarefa {
    Linha linha;
    std::vector<PointState> pontos;
};
Tarefa RodaTarefa(const Caso &c, int64_t real, const FlatMesh &m0, const std::string &raiz)
{
    const auto t0 = Clock::now();
    const Data d = DadosCaso(c);
    SlopeAnalysis an(d, MaterialKind::Unsat, WaterModel::Drawdown, FactorMethod::GIMArcLength);
    an.SetRefinement({kNunif, kNref, kPmax, 0.1});
    an.SetRecordPoints(true);
    std::unique_ptr<RandomFields> rf;
    if (real >= 0) {
        rf = std::make_unique<RandomFields>(CamposCaso(c), m0); // modos de KL do caso (Lx e Ly mudam com s)
        rf->Sample(kSemente, real);
        RandomFields *p = rf.get();
        an.SetPointProperties([p](const TPZVec<REAL> &x, TPZElastoPlasticMem &mem) {
            const FieldValues v = p->At(x[0], x[1]);
            mem.m_elastoplastic_state.fmatprop = {v.c, v.phi, v.phi}; // {c, phi, psi}: fluxo associado
        });
        an.SetPointPermeability([p](const TPZVec<REAL> &x) { return p->At(x[0], x[1]).kmult; });
    } else {
        an.SetTarget(c.det);
        an.SetOutputPrefix(PastaDet(raiz, c) + "/gim"); // VTK e CSV de cada nivel do deterministico
    }
    std::vector<LevelSummary> niveis;
    {
        QuietCout q;
        niveis = an.Run(m0);
    }
    return {Resultado(c, real < 0 ? 0 : real, real >= 0, rf.get(), niveis, Seconds(t0)), niveis.back().points};
}

void GravaDeterministico(const std::string &pasta, const Linha &l, const std::vector<PointState> &pts)
{
    GravaAtomico(pasta + "/resumo.csv", Cabecalho(l) + "\n" + Valores(l) + "\n");
    std::ostringstream m; // mecanismo.csv: os pontos com sqrt(J2(eps_p)) >= 0,01 do maximo
    m << std::setprecision(8) << "x,y,epsp_raizJ2,du,c,phi_deg,k_mult\n";
    REAL epsMax = 0.;
    for (const auto &p : pts) epsMax = std::max(epsMax, p.sqrtJ2);
    for (const auto &p : pts)
        if (epsMax > 0. && p.sqrtJ2 >= 0.01 * epsMax)
            m << p.x << "," << p.y << "," << p.sqrtJ2 << "," << p.du << "," << p.c << "," << p.phi * 180. / M_PI << ","
              << p.kmult << "\n";
    GravaAtomico(pasta + "/mecanismo.csv", m.str());
    std::string gamma, sit;
    for (const auto &[k, v] : l) {
        if (k == "FS_min") gamma = v;
        if (k == "situacao") sit = v;
    }
    GravaAtomico(pasta + "/resultado.txt", "real=0 FS_min=" + gamma + " situacao=" + sit + "\n"); // marca: concluido
}

std::vector<Caso> Seleciona(const std::string &lista)
{
    const auto todos = Casos();
    if (lista.empty() || lista == "todos") return todos;
    std::vector<Caso> v;
    std::stringstream ss(lista);
    std::string t;
    while (std::getline(ss, t, ',')) {
        bool achou = false;
        for (const auto &c : todos)
            if (c.id == t || (t.size() == 5 && t.rfind("fase", 0) == 0 && c.fase == t[4] - '0')) {
                if (std::none_of(v.begin(), v.end(), [&](const Caso &x) { return x.id == c.id; })) v.push_back(c);
                achou = true;
            }
        if (!achou) throw std::runtime_error("caso desconhecido: " + t + " (veja SlopeMonteCarloCluster casos)");
    }
    return v;
}

// ---- rodar: a fila de tarefas deste processo
int Rodar(const std::string &raiz, int64_t primeira, int64_t ultima, int k, int K, const std::vector<Caso> &casos)
{
    FlatMesh m0;
    {
        QuietCout q; // o leitor do gmsh imprime cada grupo fisico
        TPZGeoMesh *gmesh = ReadGMesh(std::string(MESH_DIR) + "/" + kMalha);
        m0 = FromGeoMesh(gmesh);
        delete gmesh;
    }
    // a fila (a mesma em todos os processos): os deterministicos e depois as realizacoes intercaladas
    std::vector<std::pair<size_t, int64_t>> fila;
    for (size_t i = 0; i < casos.size(); i++)
        if (casos[i].detProprio) fila.push_back({i, -1});
    for (int64_t r = primeira; r <= ultima; r++)
        for (size_t i = 0; i < casos.size(); i++) fila.push_back({i, r});
    std::vector<std::set<int64_t>> feitas;
    for (const auto &c : casos) {
        GravaParametros(Pasta(raiz, c), c);
        if (c.detProprio) GravaParametros(Pasta(raiz, c) + "_det", c);
        feitas.push_back(Feitas(Pasta(raiz, c)));
    }
    Print("processo %d de %d | %zu casos | realizacoes %lld..%lld | %zu tarefas na fila, %zu deste processo", k, K,
          casos.size(), (long long)primeira, (long long)ultima, fila.size(),
          fila.size() > (size_t)k ? (fila.size() - k + K - 1) / K : 0);
    int feitasAgora = 0;
    const auto t0 = Clock::now();
    for (size_t i = k; i < fila.size(); i += K) {
        if (fs::exists(raiz + "/PARAR")) {
            Print("%s/PARAR existe: nenhuma tarefa nova", raiz.c_str());
            break;
        }
        const Caso &c = casos[fila[i].first];
        const int64_t r = fila[i].second;
        const std::string pasta = Pasta(raiz, c);
        if (r < 0 ? fs::exists(PastaDet(raiz, c) + "/resultado.txt") : feitas[fila[i].first].count(r) > 0) continue;
        if (r < 0) fs::create_directories(PastaDet(raiz, c));
        const auto tr = Clock::now();
        const Tarefa t = RodaTarefa(c, r, m0, raiz);
        REAL gamma = 0.;
        for (const auto &[n, v] : t.linha)
            if (n == "FS_min") gamma = std::atof(v.c_str());
        if (r < 0) {
            GravaDeterministico(PastaDet(raiz, c), t.linha, t.pontos);
        } else { // uma linha por realizacao no arquivo deste processo
            fs::create_directories(pasta + "/partes");
            const std::string arq = pasta + "/partes/p" + std::to_string(k) + "de" + std::to_string(K) + ".csv";
            const bool novo = !fs::exists(arq) || fs::file_size(arq) == 0;
            std::ofstream out(arq, std::ios::app);
            if (novo) out << Cabecalho(t.linha) << "\n";
            out << Valores(t.linha) << std::endl;
        }
        feitasAgora++;
        Print("%-12s %s %6lld | Gamma = %.4f%s | %.1f s | %d tarefas em %.0f s", c.id.c_str(), r < 0 ? "det " : "real",
              (long long)(r < 0 ? 0 : r), gamma, gamma < 1. ? " RUPTURA" : "", Seconds(tr), feitasAgora, Seconds(t0));
    }
    return 0;
}

// ---- junta: mc_resumo.csv de cada caso (as linhas de partes/*.csv sem repeticao, em ordem de realizacao)
void Junta(const std::string &raiz)
{
    std::map<int, std::vector<std::string>> porFase;
    for (const auto &c : Casos()) {
        const std::string pasta = Pasta(raiz, c);
        if (!fs::is_directory(pasta + "/partes")) continue;
        std::string cab;
        const auto linhas = LinhasPartes(pasta, cab);
        std::string txt = cab + "\n";
        for (const auto &[r, l] : linhas) txt += l + "\n";
        GravaAtomico(pasta + "/mc_resumo.csv", txt);
        porFase[c.fase].push_back(c.id + " " + std::to_string(linhas.size()));
    }
    for (const auto &[f, v] : porFase) { // campanha.txt: os casos da fase e as realizacoes de cada um
        std::string txt;
        for (const auto &l : v) txt += l + "\n";
        GravaAtomico(raiz + "/" + Prefixo() + "_fase" + std::to_string(f) + "/campanha.txt", txt);
    }
}

// Colunas de um CSV: nome -> valores, de todas as linhas
std::map<std::string, std::vector<std::string>> LeCsv(const std::string &arq)
{
    std::map<std::string, std::vector<std::string>> col;
    std::ifstream in(arq);
    std::string l;
    std::vector<std::string> nomes;
    auto separa = [](const std::string &s) {
        std::vector<std::string> v;
        std::stringstream ss(s);
        std::string t;
        while (std::getline(ss, t, ',')) v.push_back(t);
        return v;
    };
    if (std::getline(in, l)) nomes = separa(l);
    while (std::getline(in, l)) {
        if (l.empty()) continue;
        const auto v = separa(l);
        for (size_t i = 0; i < nomes.size(); i++) col[nomes[i]].push_back(i < v.size() ? v[i] : "");
    }
    return col;
}

// ---- resumo: estatistica de cada caso (como o dados_sat.py: ajuste lognormal por maxima verossimilhanca) e o artigo
void Resumo(const std::string &raiz)
{
    Junta(raiz);
    const std::string linha(132, '=');
    std::cout << linha << "\n MONTE CARLO DO CERON | " << raiz << " | " << Prefixo() << "_fase{1,2,3} | malha " << kMalha
              << ", " << kNunif << " unif. + " << kNref << " adapt., u ate P" << kPmax << "\n" << linha << "\n";
    Print(" %-12s %6s %5s %8s %8s %7s %7s %6s %8s | %7s %7s %6s %8s %8s", "caso", "S", "n_f", "P_f", "CoV(P_f)", "media",
          "desvio", "CoV", "det", "media", "desvio", "CoV", "P_f", "det");
    Print(" %-12s %54s %8s | %47s", "", "este trabalho", "", "artigo (Tabelas 3 a 6)");
    REAL detRef = 0.;
    for (const auto &c : Casos()) {
        const std::string arq = Pasta(raiz, c) + "/mc_resumo.csv";
        REAL det = NAN;
        const Caso &cd = c.detProprio ? c : Casos().front();
        const auto dcol = LeCsv(PastaDet(raiz, cd) + "/resumo.csv");
        if (dcol.count("FS_min") && !dcol.at("FS_min").empty()) det = std::atof(dcol.at("FS_min")[0].c_str());
        if (c.id == "ref") detRef = det;
        if (!fs::exists(arq)) continue;
        const auto col = LeCsv(arq);
        std::vector<REAL> g;
        int64_t fora = 0;
        for (size_t i = 0; i < col.at("FS_min").size(); i++) {
            const REAL v = std::atof(col.at("FS_min")[i].c_str());
            if (col.at("situacao")[i] == "ok" && v > 0.) g.push_back(v);
            else fora++;
        }
        if (g.empty()) continue;
        const REAL S = g.size();
        const int64_t nf = std::count_if(g.begin(), g.end(), [](REAL v) { return v < 1.; });
        const REAL pf = nf / S, covpf = nf > 0 ? std::sqrt((1. - pf) / (S * pf)) : NAN;
        REAL ml = 0., sl = 0.;
        for (REAL v : g) ml += std::log(v) / S;
        for (REAL v : g) sl += (std::log(v) - ml) * (std::log(v) - ml) / S;
        sl = std::sqrt(sl);
        const REAL media = std::exp(ml + 0.5 * sl * sl), cov = std::sqrt(std::exp(sl * sl) - 1.);
        const REAL detArt = c.det > 0. ? c.det : 1.336;
        Print(" %-12s %6.0f %5lld %7.2f%% %7.2f%% %7.3f %7.3f %5.1f%% %8.4f | %7.3f %7.3f %5.1f%% %7.2f%% %8.3f%s",
              c.id.c_str(), S, (long long)nf, 100. * pf, 100. * covpf, media, media * cov, 100. * cov, det, c.mu, c.sigma,
              c.cov, c.pf, detArt, fora ? ("  (" + std::to_string(fora) + " fora)").c_str() : "");
    }
    std::cout << linha << "\n media e desvio: os da lognormal ajustada (como no artigo); det: o Gamma com as medias (o da"
                          " referencia nos casos de CoV e s, que nao o mudam)\n";
    if (std::isfinite(detRef))
        Print(" deterministico da referencia: %.4f (artigo 1,336, %+.1f%%): o vies da malha, que desloca a media de cada caso",
              detRef, 100. * (detRef / 1.336 - 1.));
    std::cout << linha << std::endl;
}

int main(int argc, char *argv[])
{
    const std::string modo = argc > 1 ? argv[1] : "";
    try {
        if (modo == "casos") { // os casos (todos ou os da lista): tambem confere a lista antes de submeter
            for (const auto &c : Seleciona(argc > 2 ? argv[2] : "todos"))
                Print("%-12s fase %d  %s%s%s", c.id.c_str(), c.fase, kBase.c_str(), c.sufixo.c_str(),
                      c.detProprio ? "  (com deterministico proprio)" : "");
            return 0;
        }
        if (modo == "rodar" && argc >= 7) {
            const int k = std::atoi(argv[5]), K = std::atoi(argv[6]);
            const int64_t a = std::atoll(argv[3]), b = std::atoll(argv[4]);
            if (K < 1 || k < 0 || k >= K || a < 0 || b < a) throw std::runtime_error("argumentos invalidos");
            return Rodar(argv[2], a, b, k, K, Seleciona(argc > 7 ? argv[7] : "todos"));
        }
        if (modo == "junta" && argc >= 3) {
            Junta(argv[2]);
            return 0;
        }
        if (modo == "resumo" && argc >= 3) {
            Resumo(argv[2]);
            return 0;
        }
    } catch (const std::exception &e) {
        std::cerr << "erro: " << e.what() << "\n";
        return 2;
    }
    std::cerr << "Uso: SlopeMonteCarloCluster casos [casos = todos]\n"
                 "     SlopeMonteCarloCluster rodar <raiz> <primeira> <ultima> <k> <K> [casos = todos]\n"
                 "     SlopeMonteCarloCluster junta <raiz>\n"
                 "     SlopeMonteCarloCluster resumo <raiz>\n";
    return 1;
}
