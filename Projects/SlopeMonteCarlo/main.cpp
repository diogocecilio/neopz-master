// SlopeMonteCarlo: Monte Carlo do talude de Vargas Ceron et al. (2025), Secao 6, na configuracao de referencia (Tabela 2,
// Fig. 25), com a biblioteca Common:
//   - talude de H = 5 m e beta = 45 graus, rebaixamento total (h_w = H), malha inicial ceron_h5_a45_g15 (126 triangulos);
//   - campos aleatorios lognormais, independentes entre si: c' (10 kPa, CoV 30 %), phi' (30 graus, CoV 10 %) e k_v
//     (1e-6 m4/kN/s, CoV 60 %), com Lx = 20 m e Ly = 2 m, por Karhunen-Loeve com M = 2000 termos (RandomFields);
//   - em cada realizacao (SlopeAnalysis): percolacao permanente com o k_v da realizacao, Gamma pelo GIM com arc-length e
//     as pressoes congeladas, e refinamento adaptativo da zona plastica. O material e o nao saturado (saturado em qualquer
//     pressao), o unico com condutividade por ponto. O Gamma da realizacao e o do ultimo nivel.
// Ruptura: Gamma < 1 (Eq. 59). P_f = n_f / S (Eq. 67) e CoV(P_f) = sqrt((1 - P_f) / (S P_f)) (Eq. 68).
// Artigo: S = 10000, P_f = 11,50 % (CoV 2,77 %); Gamma lognormal com media 1,353, mediana 1,313 e CoV 23,54 %.
//
// Uso: SlopeMonteCarlo [nreal = 100] [primeira = 0] [nunif = 0] [nref = 4] [pmax = 2]
//      SlopeMonteCarlo resumo [nunif = 0] [nref = 4] [pmax = 2]
// Roda as realizacoes primeira .. primeira + nreal - 1 e grava uma linha por realizacao em
// mc_u<nunif>h<nref>p<pmax>_<primeira>_<ultima>.csv, na pasta corrente. A realizacao -1 e o deterministico (as medias),
// rodado pelo bloco que comeca em 0. Retomavel: as realizacoes que ja estao em algum mc_u..h..p.._*.csv da pasta sao
// puladas, e a estatistica do fim (e a do modo resumo) junta todos eles. Assim varios processos dividem as realizacoes
// (rodar_mc.sh) e o Monte Carlo pode ser estendido depois com mais blocos.
#include "SlopeAnalysis.h"
#include "SlopeFields.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>

using namespace slope;

const uint64_t kSemente = 20250416; // a do SlopeRandomField: com o mesmo retangulo e M, a realizacao n e a mesma de la

// ---- leitura dos resultados: realizacao -> Gamma, de todos os CSV desta configuracao na pasta corrente
std::map<int64_t, REAL> ReadResults(const std::string &config)
{
    std::map<int64_t, REAL> res;
    const std::string prefixo = "mc_" + config + "_";
    for (const auto &e : std::filesystem::directory_iterator(".")) {
        const std::string nome = e.path().filename().string();
        if (nome.rfind(prefixo, 0) != 0 || e.path().extension() != ".csv") continue;
        std::ifstream f(e.path());
        std::string linha;
        std::getline(f, linha); // cabecalho
        while (std::getline(f, linha)) {
            std::stringstream ss(linha);
            std::string real, gamma;
            if (std::getline(ss, real, ',') && std::getline(ss, gamma, ','))
                res[std::stoll(real)] = std::stod(gamma);
        }
    }
    return res;
}

// ---- estatistica do Monte Carlo (Eqs. 67 e 68) e ajuste lognormal (por maxima verossimilhanca, como no artigo)
void PrintStatistics(const std::map<int64_t, REAL> &res, const std::string &config)
{
    std::vector<REAL> g, lg;
    for (const auto &[r, v] : res)
        if (r >= 0) g.push_back(v);
    const std::string linha(96, '=');
    std::cout << "\n" << linha << "\n MONTE CARLO | " << config << " | " << g.size() << " realizacoes nos mc_" << config
              << "_*.csv da pasta\n" << linha << "\n";
    if (res.count(-1)) Print(" deterministico (medias): Gamma = %.5f", res.at(-1));
    if (g.empty()) return;
    const REAL S = g.size();
    const int64_t nf = std::count_if(g.begin(), g.end(), [](REAL v) { return v < 1.; });
    const REAL Pf = nf / S, covPf = nf > 0 ? std::sqrt((1. - Pf) / (S * Pf)) : 0.;
    REAL media = 0., dp = 0.;
    for (REAL v : g) media += v / S;
    for (REAL v : g) dp += (v - media) * (v - media) / std::max(S - 1., 1.);
    dp = std::sqrt(dp);
    std::vector<REAL> o(g);
    std::sort(o.begin(), o.end());
    const REAL mediana = o.size() % 2 ? o[o.size() / 2] : 0.5 * (o[o.size() / 2 - 1] + o[o.size() / 2]);
    for (REAL v : g)
        if (v > 0.) lg.push_back(std::log(v)); // Gamma = 0: sem equilibrio inicial (conta em P_f, fica fora do ajuste)
    REAL mu = 0., s = 0.;
    for (REAL v : lg) mu += v / lg.size();
    for (REAL v : lg) s += (v - mu) * (v - mu) / lg.size();
    s = std::sqrt(s);
    const REAL PfLn = 0.5 * std::erfc(mu / (s * std::sqrt(2.))); // P(Gamma < 1) da lognormal
    Print(" %-26s %7s %6s %8s %9s %8s %8s %8s", "", "S", "n_f", "P_f", "CoV(P_f)", "media", "mediana", "CoV");
    Print(" %-26s %7.0f %6lld %7.2f%% %8.2f%% %8.4f %8.4f %7.2f%%", "este (amostra)", S, (long long)nf, 100. * Pf,
          100. * covPf, media, mediana, 100. * dp / media);
    Print(" %-26s %7s %6s %7.2f%% %9s %8.4f %8.4f %7.2f%%", "este (ajuste lognormal)", "", "", 100. * PfLn, "",
          std::exp(mu + 0.5 * s * s), std::exp(mu), 100. * std::sqrt(std::exp(s * s) - 1.));
    Print(" %-26s %7d %6d %7.2f%% %8.2f%% %8.4f %8.4f %7.2f%%", "artigo (Fig. 25, lognormal)", 10000, 1150, 11.50, 2.77,
          1.353, 1.313, 23.54);
    Print(" ajuste lognormal: mu_ln = %.4f, sigma_ln = %.4f; menor Gamma %.4f, maior %.4f", mu, s, o.front(), o.back());
    std::cout << linha << std::endl;
}

int main(int argc, char *argv[])
{
    const bool resumo = argc > 1 && std::string(argv[1]) == "resumo";
    const int a = resumo ? 2 : 3; // posicao de nunif nos argumentos
    const int nreal = !resumo && argc > 1 ? std::atoi(argv[1]) : 100;
    const int64_t primeira = !resumo && argc > 2 ? std::atoll(argv[2]) : 0;
    const int nunif = argc > a ? std::atoi(argv[a]) : 0;    // refinamentos uniformes da malha inicial
    const int nref = argc > a + 1 ? std::atoi(argv[a + 1]) : 4; // niveis de refinamento adaptativo
    const int pmax = argc > a + 2 ? std::atoi(argv[a + 2]) : 2; // ordem maxima de u (2: so h)
    if (nreal < 0 || primeira < 0 || nunif < 0 || nref < 0 || pmax < 2) {
        std::cerr << "Uso: SlopeMonteCarlo [nreal = 100] [primeira = 0] [nunif = 0] [nref = 4] [pmax = 2]\n"
                     "     SlopeMonteCarlo resumo [nunif = 0] [nref = 4] [pmax = 2]\n";
        return 1;
    }
    const std::string config = "u" + std::to_string(nunif) + "h" + std::to_string(nref) + "p" + std::to_string(pmax);
    if (resumo) {
        PrintStatistics(ReadResults(config), config);
        return 0;
    }

    // ---- passo 1: malha inicial e o talude do artigo (o padrao de Data: c' = 10 kPa, phi' = 30 graus, gamma = 20 kN/m3,
    // k_v = 1e-6 m4/kN/s, H = 5 m, beta = 45 graus, crista em (10, 10), h_w = H)
    const std::string arquivo = std::string(MESH_DIR) + "/ceron_h5_a45_g15.msh";
    FlatMesh m0;
    {
        QuietCout q; // o leitor do gmsh imprime cada grupo fisico
        TPZGeoMesh *gmesh = ReadGMesh(arquivo);
        m0 = FromGeoMesh(gmesh);
        delete gmesh;
    }
    const Data d;

    // ---- passo 2: os modos de Karhunen-Loeve, uma vez, no retangulo da malha (0 a 25 m por 0 a 10 m)
    const FieldStats stats; // Tabela 2
    RandomFields campos(stats, m0);
    Print("Monte Carlo do talude do Ceron | malha %s (%zu triangulos) | %d unif. + %d adapt., u ate P%d",
          arquivo.substr(arquivo.rfind('/') + 1).c_str(), m0.tri.size(), nunif, nref, pmax);
    Print("campos: c' %.0f kPa (CoV %.0f%%), phi' %.0f graus (CoV %.0f%%), k_v (CoV %.0f%%), Lx = %.0f m, Ly = %.0f m | "
          "KL: M = %d termos, %.1f%% da variancia", stats.c, 100. * stats.covC, stats.phiDeg, 100. * stats.covPhi,
          100. * stats.covK, stats.Lx, stats.Ly, campos.M(), 100. * campos.VarianceFraction());

    // ---- passo 3: o arquivo do bloco; as realizacoes que ja estao em algum CSV desta configuracao sao puladas
    const int64_t ultima = primeira + nreal - 1;
    const std::string nomeCsv = "mc_" + config + "_" + std::to_string(primeira) + "_" + std::to_string(ultima) + ".csv";
    std::set<int64_t> feitas;
    for (const auto &[r, v] : ReadResults(config)) feitas.insert(r);
    const bool novo = !std::filesystem::exists(nomeCsv);
    std::ofstream csv(nomeCsv, std::ios::app);
    csv << std::setprecision(8);
    if (novo) {
        csv << "real,gamma";
        for (int k = 0; k <= nref; k++) csv << ",gamma_nivel" << k;
        csv << ",triangulos,equacoes,passos,segundos" << std::endl;
    }
    Print("arquivo %s | %zu realizacoes ja feitas na pasta", nomeCsv.c_str(), feitas.size());

    // ---- passo 4: as realizacoes. Cada uma sorteia os campos e roda a SlopeAnalysis com eles; o campo e avaliado na
    // posicao de cada ponto de integracao, entao todos os niveis de refinamento veem a mesma realizacao
    std::vector<int64_t> lista;
    if (primeira == 0) lista.push_back(-1); // o deterministico, com as medias
    for (int64_t r = primeira; r <= ultima; r++) lista.push_back(r);
    int64_t n = 0, nf = 0;
    const auto t0 = Clock::now();
    for (int64_t r : lista) {
        if (feitas.count(r)) continue;
        const auto tr = Clock::now();
        SlopeAnalysis an(d, MaterialKind::Unsat, WaterModel::Drawdown, FactorMethod::GIMArcLength);
        an.SetRefinement({nunif, nref, pmax, 0.1});
        if (r >= 0) {
            campos.Sample(kSemente, r);
            an.SetPointProperties([&campos](const TPZVec<REAL> &x, TPZElastoPlasticMem &mem) {
                const FieldValues v = campos.At(x[0], x[1]);
                mem.m_elastoplastic_state.fmatprop = {v.c, v.phi, v.phi}; // {c, phi, psi}: fluxo associado
            });
            an.SetPointPermeability([&campos](const TPZVec<REAL> &x) { return campos.At(x[0], x[1]).kmult; });
        }
        std::vector<LevelSummary> niveis;
        if (r < 0) {
            std::cout << "\n--- deterministico (medias) ---" << std::endl;
            an.SetTarget(1.336);
            niveis = an.Run(m0);
        } else {
            QuietCout q; // uma linha por realizacao, abaixo
            niveis = an.Run(m0);
        }
        const LevelSummary &fim = niveis.back();
        csv << r << "," << fim.value;
        for (const auto &l : niveis) csv << "," << l.value;
        csv << "," << fim.ntri << "," << fim.neq << "," << fim.steps << "," << Seconds(tr) << std::endl;
        std::string valores;
        for (const auto &l : niveis) {
            char b[16];
            std::snprintf(b, sizeof b, " %.3f", l.value);
            valores += b;
        }
        if (r >= 0) n++, nf += fim.value < 1.;
        Print("real %6lld | Gamma = %.4f%s (niveis:%s) | %lld triang. | %.1f s | bloco: P_f = %lld/%lld | %.0f s",
              (long long)r, fim.value, fim.value < 1. ? " RUPTURA" : "", valores.c_str(), (long long)fim.ntri,
              Seconds(tr), (long long)nf, (long long)n, Seconds(t0));
    }
    PrintStatistics(ReadResults(config), config);
    return 0;
}
