// Teste dos quatro metodos do fator (GIM e SRM, por arc-length e por bisseccao) no talude seco de H = 10 m e 45 graus
// (Fig. 8.30 de de Souza Neto, 75 m de largura), com c = 10 kPa, phi = 30 graus e gamma = 20 kN/m3: o talude de Cho da
// Sec. 5.3.2 de Vargas Ceron et al. (2025), Gamma = 1,777 e FS = 1,203.
// Uso: SlopeCommonMetodos [nunif = 0] [nref = 2] [pmax = 3] [material: up | unsat] [props]
// props: impoe em cada ponto de integracao as proprias medias {c, phi, psi} (o gancho dos campos aleatorios): os fatores
// tem que ser os mesmos
#include "SlopeAnalysis.h"

#include <cstdlib>

using namespace slope;

int main(int argc, char *argv[])
{
    const int nunif = argc > 1 ? std::atoi(argv[1]) : 0, nref = argc > 2 ? std::atoi(argv[2]) : 2;
    const int pmax = argc > 3 ? std::atoi(argv[3]) : 3;
    const MaterialKind kind = argc > 4 && std::string(argv[4]) == "unsat" ? MaterialKind::Unsat : MaterialKind::UP;
    TPZGeoMesh *gmesh = ReadGMesh(MESH_TALUDE10);
    const FlatMesh m0 = FromGeoMesh(gmesh);
    delete gmesh;
    Data d;
    d.H = 10., d.xCrista = 35., d.yCrista = 40., d.hw = 10.;

    std::vector<SummaryRun> runs;
    const auto t0 = Clock::now();
    for (FactorMethod f : {FactorMethod::GIMArcLength, FactorMethod::GIMBisection, FactorMethod::SRMArcLength,
                           FactorMethod::SRMBisection}) {
        const bool gim = f == FactorMethod::GIMArcLength || f == FactorMethod::GIMBisection;
        const REAL alvo = gim ? 1.777 : 1.203;
        std::cout << "\n=== " << Name(f) << ", material " << Name(kind) << " (alvo " << alvo << ") ===" << std::endl;
        SlopeAnalysis a(d, kind, WaterModel::Dry, f);
        a.SetRefinement({nunif, nref, pmax, 0.1});
        a.SetTarget(alvo);
        if (argc > 5 && std::string(argv[5]) == "props")
            a.SetPointProperties([&d](const TPZVec<REAL> &, TPZElastoPlasticMem &mem) {
                mem.m_elastoplastic_state.fmatprop.Resize(3);
                mem.m_elastoplastic_state.fmatprop[0] = d.c;
                mem.m_elastoplastic_state.fmatprop[1] = d.phi;
                mem.m_elastoplastic_state.fmatprop[2] = d.phi;
            });
        runs.push_back({Name(f), alvo, a.Run(m0)});
    }
    PrintSummary("talude seco H = 10 m, material " + std::string(Name(kind)) + ", " + std::to_string(nunif) +
                     " unif. + " + std::to_string(nref) + " adapt., u ate P" + std::to_string(pmax),
                 runs, Seconds(t0));
    return 0;
}
