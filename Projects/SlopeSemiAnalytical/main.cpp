// Rebaixamento rapido total (h_w = H) do talude do relatorio (malha ceron_h5_a45_g15, H = 5 m, 45 graus) por GIM com
// arc-length e refinamento hp, com tres modelos da agua:
//   1. rebaix:        percolacao permanente do MEF, pressoes congeladas (x lambda)
//   2. rebaixforca:   a mesma percolacao como forca de corpo -grad u com o peso submerso (verifica a forca de corpo:
//                     tem que dar o mesmo Gamma de 1)
//   3. semianalitico: o campo semi-analitico de Vargas Ceron et al. (2025), Sec. 4.1, como forca de corpo
// Alvo: Gamma = 1,336 (Tabelas 5 e 6 do artigo).
// Uso: SlopeSemiAnalytical [nunif = 1] [nref = 3] [pmax = 3] [materiais: ambos | up | unsat]
#include "SemiAnalyticalSeepage.h"
#include "SlopeAnalysis.h"

#include <cstdlib>

using namespace slope;

int main(int argc, char *argv[])
{
    const int nunif = argc > 1 ? std::atoi(argv[1]) : 1, nref = argc > 2 ? std::atoi(argv[2]) : 3;
    const int pmax = argc > 3 ? std::atoi(argv[3]) : 3;
    const std::string quais = argc > 4 ? argv[4] : "ambos";
    std::vector<MaterialKind> materiais;
    if (quais != "unsat") materiais.push_back(MaterialKind::UP);
    if (quais != "up") materiais.push_back(MaterialKind::Unsat);
    const REAL alvo = 1.336;

    TPZGeoMesh *gmesh = ReadGMesh(std::string(MESH_DIR) + "/ceron_h5_a45_g15.msh");
    const FlatMesh m0 = FromGeoMesh(gmesh);
    delete gmesh;
    const Data d; // o talude do relatorio, rebaixamento total

    const std::string config = std::to_string(nunif) + " unif. + " + std::to_string(nref) + " adapt., u ate P" +
                               std::to_string(pmax);
    const std::string linha(100, '=');
    std::cout << linha << "\n Rebaixamento total, talude do relatorio (H = 5 m) | GIM arc-length | " << config
              << " | alvo Gamma = " << alvo << "\n" << linha << std::endl;
    const SemiAnalyticalSeepage campo(d, 10. * d.H);
    campo.Print();

    const std::vector<std::pair<WaterModel, std::string>> modelos = {
        {WaterModel::Drawdown, "rebaix"}, {WaterModel::DrawdownForce, "rebaixforca"},
        {WaterModel::ForceField, "semianalitico"}};
    std::vector<SummaryRun> runs;
    const auto t0 = Clock::now();
    for (const auto &[agua, nome] : modelos)
        for (MaterialKind kind : materiais) {
            std::cout << "\n--- " << nome << " (" << Name(agua) << "), material " << Name(kind) << " ---" << std::endl;
            SlopeAnalysis a(d, kind, agua, FactorMethod::GIMArcLength);
            a.SetRefinement({nunif, nref, pmax, 0.1});
            a.SetTarget(alvo);
            a.SetOutputPrefix("gim_" + nome + "_ceron");
            if (agua == WaterModel::ForceField)
                a.SetSeepageForce([&campo](const TPZVec<REAL> &x) { return campo.Force(x); });
            runs.push_back({nome + " " + Name(kind), alvo, a.Run(m0)});
        }
    PrintSummary("rebaixamento total, H = 5 m | GIM arc-length | " + config, runs, Seconds(t0));
    return 0;
}
