// Talude do relatorio (Vargas Ceron et al. 2025) passo a passo, com a biblioteca Common (slope::SlopeAnalysis):
//   passo 1: malha geometrica do gmsh;
//   passo 2: malhas computacionais dos dois materiais (TPZMatPoroElastoPlasticUP e TPZMatUnsaturatedUPwPa);
//   passo 3 (seco) e 4 (rebaixamento total): GIM com arc-length e refinamento adaptativo hp, nos dois materiais;
//   passo 5: chuva de Rahardjo et al. (2007), solo f50,-5 e alfa = 45 graus: transiente nao saturado e FS ao longo do tempo.
// Uso: SlopeRandomFieldRefinement [ceron | talude10 | ceron10H] [seco | rebaix] [nunif = 1] [nref = 3] [pmax = 3]
//      SlopeRandomFieldRefinement rahardjo [I_r = 36 mm/h] [nunif = 0] [nref = 0] [pmax = 3]
#include "SlopeAnalysis.h"

#include "TPZVTKGeoMesh.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>

using namespace slope;

// ---- passo 5: chuva de Rahardjo et al. (2007), Serie A, solo f50,-5 e alfa = 45 graus (Fig. 4c2)
// Talude da Fig. 2 (malha rahardjo_a45 do LearningMohrCoulomb2): H_s = 10 m, pe em (40, 20), crista em (30, 30), 70 m de
// largura, laterais de 30 m e 20 m. Resistencia da Eq. 5: c' = 10 kPa, phi' = phi^b = 26 graus; gamma = 20 kN/m3; E e nu
// nao estao no artigo (2e4 kPa e 0,3; nao mudam o FS). Solo, chuva e nivel d'agua inicial: RainData (o padrao e o f50,-5).
// O FS do artigo e o de Bishop simplificado (SLOPE/W) com as pressoes do SEEP/W; o equivalente aqui e o SRM: em cada
// instante, as pressoes do transiente ficam congeladas e c' e tan(phi') (e com eles a parcela da succao) sao divididos por F.
// Referencias do FS x t para cada I_r (mm/h), em t (h): o artigo, lido da Fig. 4c2 (+-0,03), e o driver antigo do
// LearningMohrCoulomb2 (cluster/fonte/rahardjo.cpp, variante iii = as mesmas hipoteses; SRM por bisseccao, tolerancia 0,01)
const std::map<int, std::map<int, REAL>> kArtigo = {
    {9, {{0, 2.12}, {24, 1.78}}}, {36, {{0, 2.12}, {24, 1.05}, {48, 1.23}}}, {80, {{0, 2.12}, {24, 0.83}}}};
const std::map<int, std::map<int, REAL>> kCluster = {
    {9, {{0, 2.11875}, {6, 2.01875}, {12, 1.94375}, {18, 1.875}, {24, 1.8125}, {36, 1.85}, {48, 1.875}}},
    {36, {{0, 2.11875}, {6, 1.75}, {12, 1.425}, {18, 1.20625}, {24, 1.0625}, {36, 1.1875}, {48, 1.28125}}},
    {80, {{0, 2.11875}, {6, 1.39375}, {12, 1.175}, {18, 1.05}, {24, 0.85}, {36, 1.0625}, {48, 1.175}}}};

// Valor de referencia em t para a chuva I, ou "-"
static std::string Reference(const std::map<int, std::map<int, REAL>> &tabela, int I, REAL t)
{
    char buf[16] = "        -";
    const long h = std::lround(t);
    if (tabela.count(I) && std::fabs(t - h) < 1.e-6 && tabela.at(I).count(h))
        std::snprintf(buf, sizeof buf, "%9.3f", tabela.at(I).at(h));
    return buf;
}

int Passo5Chuva(int argc, char *argv[])
{
    const REAL intensidade = argc > 2 ? std::atof(argv[2]) : 36.; // I_r (mm/h); 36 = 1 k_s do f50,-5
    const int nunif = argc > 3 ? std::atoi(argv[3]) : 0;
    const int nref = argc > 4 ? std::atoi(argv[4]) : 0;
    const int pmax = argc > 5 ? std::atoi(argv[5]) : 3;
    if (!(intensidade >= 0.) || nunif < 0 || nref < 0 || pmax < 2) {
        std::cerr << "Uso: SlopeRandomFieldRefinement rahardjo [I_r = 36 mm/h] [nunif = 0] [nref = 0] [pmax = 3]\n";
        return 1;
    }
    Data d; // resistencia e peso do solo; a agua vem toda da RainData
    d.c = 10.;
    d.phi = 26. * M_PI / 180.;
    d.gamma = 20.;
    d.H = 10., d.xCrista = 30., d.yCrista = 30.;
    RainData chuva; // solo f50,-5, 24 h de chuva, nivel inicial 5 m abaixo do pe com 7 graus, succao ate 75 kPa
    chuva.intensity = intensidade;

    const std::string arquivo = std::string(MESH_DIR) + "/rahardjo_a45.msh";
    std::cout << "Caso rahardjo: lendo " << arquivo << "\n";
    TPZGeoMesh *gmesh = ReadGMesh(arquivo);
    const FlatMesh m0 = FromGeoMesh(gmesh);
    delete gmesh;
    Print("solo f50,-5: k_s = %.0e m/s, Fredlund-Xing a = %.0f kPa (n = m = 1), k_w = k_s Theta^%.0f, phi^b = %.0f graus, "
          "n = %.2f", chuva.ks, chuva.aFX, chuva.power, chuva.phib * 180. / M_PI, chuva.n);
    Print("nivel inicial: y = %.2f m em x = 0 e %.2f m em x = 70; succao ate %.0f kPa", chuva.WaterLevel(0.),
          chuva.WaterLevel(70.), chuva.suctionMax);

    const int I = int(std::lround(intensidade));
    const REAL alvo = kArtigo.count(I) && kArtigo.at(I).count(24) ? kArtigo.at(I).at(24) : 0.; // FS_min do artigo
    char titulo[160];
    std::snprintf(titulo, sizeof titulo, "rahardjo f50,-5, alfa = 45 | I_r = %g mm/h por %g h | %d unif. + %d adapt., u ate P%d",
                  intensidade, chuva.duration, nunif, nref, pmax);
    const std::string prefixo = "chuva_rahardjo_I" + std::to_string(I) + "_u" + std::to_string(nunif) + "h" +
                                std::to_string(nref) + "p" + std::to_string(pmax);
    const std::string linha(100, '=');
    std::cout << "\n" << linha << "\n SRM ao longo da chuva | " << titulo << " | alvo FS_min (artigo, t = 24 h) = " << alvo
              << "\n" << linha << std::endl;
    const auto t0 = Clock::now();
    SlopeAnalysis a(d, MaterialKind::Unsat, WaterModel::Rain, FactorMethod::SRMArcLength);
    a.SetRain(chuva);
    a.SetRefinement({nunif, nref, pmax, 0.1}); // o refinamento segue a zona plastica do instante de menor FS
    a.SetTarget(alvo);
    a.SetOutputPrefix(prefixo); // <prefixo>_unsat_nivel<k>_chuva.csv (FS x t) e _t<h>h.vtk (p_w e S_w)
    const std::vector<SummaryRun> runs = {{"unsat SRM", alvo, a.Run(m0)}};
    PrintSummary(titulo, runs, Seconds(t0));

    // FS x t de cada nivel, com as referencias
    const auto &niveis = runs[0].levels;
    std::cout << "\n FS x t | artigo: lido da Fig. 4c2 (+-0,03) | cluster: driver antigo (LearningMohrCoulomb2)\n";
    std::string cab = "   t (h)";
    for (const auto &l : niveis) cab += "   nivel " + std::to_string(l.level);
    std::cout << cab << "   artigo  cluster\n";
    for (size_t i = 0; i < niveis[0].history.size(); i++) {
        const REAL t = niveis[0].history[i][0];
        char buf[32];
        std::snprintf(buf, sizeof buf, "%8.1f", t);
        std::string lin = buf;
        for (const auto &l : niveis) {
            if (i >= l.history.size()) std::snprintf(buf, sizeof buf, "%10s", "-");
            else if (l.history[i][1] <= 0.) std::snprintf(buf, sizeof buf, "%10s", "< 1"); // SRM sem o peso total
            else std::snprintf(buf, sizeof buf, "%10.4f", l.history[i][1]);
            lin += buf;
        }
        std::cout << lin << Reference(kArtigo, I, t) << Reference(kCluster, I, t) << "\n";
    }
    return 0;
}



// Casos (argumento do programa): malha, geometria e alvos. Seco: o 1,777 da Secao 5.3.2, que e do talude de H = 10 m
// da Fig. 19 (o texto diz 5 m); como Gamma = N c / (gamma H), com H = 5 m o alvo e 2 x 1,777 = 3,554. Rebaixamento
// total (h_w = H): Gamma = 1,336 (Tabelas 5 e 6); no talude10 so 1,336 / 2 (Gamma ~ 1/H), indicativo
struct Case {
    std::string mesh;
    REAL H, xCrista, yCrista, targetDry, targetDrawdown;
};
const std::map<std::string, Case> kCases = {
    {"ceron", {std::string(MESH_DIR) + "/ceron_h5_a45_g15.msh", 5., 10., 10., 3.554, 1.336}},          // relatorio
    {"talude10", {std::string(PROJECT_DIR) + "/malhas/talude_h10_grossa.msh", 10., 35., 40., 1.777, 0.668}}, // Fig. 8.30
    {"ceron10H", {std::string(PROJECT_DIR) + "/malhas/ceron_h5_a45_10H.msh", 5., 50., 55., 3.554, 1.336}}};  // 10H

int main(int argc, char *argv[])
{
    if (argc > 1 && std::string(argv[1]) == "rahardjo") return Passo5Chuva(argc, argv); // passo 5: chuva
    const std::string caso = argc > 1 ? argv[1] : "ceron";
    const std::string modo = argc > 2 ? argv[2] : "rebaix";
    const int nunif = argc > 3 ? std::atoi(argv[3]) : 1; // refinamentos uniformes da malha inicial
    const int nref = argc > 4 ? std::atoi(argv[4]) : 3;  // niveis de refinamento adaptativo
    const int pmax = argc > 5 ? std::atoi(argv[5]) : 3;  // ordem maxima de u no refinamento p
    if (!kCases.count(caso) || (modo != "seco" && modo != "rebaix") || nunif < 0 || nref < 0 || pmax < 2) {
        std::cerr << "Uso: SlopeRandomFieldRefinement [ceron | talude10 | ceron10H] [seco | rebaix] [nunif = 1] "
                     "[nref = 3] [pmax = 3]\n     SlopeRandomFieldRefinement rahardjo [I_r = 36 mm/h] [nunif = 0] "
                     "[nref = 0] [pmax = 3]\n";
        return 1;
    }
    const Case &cs = kCases.at(caso);
    const bool rebaix = modo == "rebaix";
    const REAL alvo = rebaix ? cs.targetDrawdown : cs.targetDry;

    // ---- passo 1: malha geometrica
    const std::string nome = cs.mesh.substr(cs.mesh.rfind('/') + 1);
    std::cout << "Caso " << caso << ": lendo " << cs.mesh << "\n";


    TPZGeoMesh *gmesh = ReadGMesh(cs.mesh);

    const std::string vtk = "malha_" + nome.substr(0, nome.rfind('.')) + ".vtk";
    std::ofstream out(vtk);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); // cada elemento com o seu id de material
    std::cout << "Gravado " << vtk << " (abra no ParaView)\n";
    DuplicateBoundaryBC(gmesh); // gemeos -12 e -14, antes das malhas computacionais

    // ---- passo 2: malhas computacionais com os dois materiais
    Data d;
    d.H = cs.H, d.xCrista = cs.xCrista, d.yCrista = cs.yCrista;
    d.hw = cs.H; // rebaixamento total
    std::cout << "\nMalha u-p (TPZMatPoroElastoPlasticUP):\n";
    TPZMultiphysicsCompMesh *malhaUP = CreateUPCompMesh(gmesh, 2, d);

    std::cout << "\nMalha u-p_w-p_a (TPZMatUnsaturatedUPwPa):\n";
    TPZMultiphysicsCompMesh *malhaUnsat = CreateUnsatCompMesh(gmesh, 2, d, d.gamma);

    EraseMesh(malhaUP);
    EraseMesh(malhaUnsat);

    // ---- passo 3 (seco) e 4 (rebaixamento): GIM com arc-length e refinamento adaptativo hp
    const FlatMesh m0 = FromGeoMesh(gmesh); // sem os gemeos; os refinamentos uniformes ficam com a SlopeAnalysis
    delete gmesh;
    const std::string config = std::to_string(nunif) + " unif. + " + std::to_string(nref) + " adapt., u ate P" +
                               std::to_string(pmax);
    const std::string titulo = "caso " + caso + " | " + (rebaix ? "rebaixamento (h_w = H)" : "seco") + " | " + config;
    const std::string prefixo = "gim_" + modo + "_" + caso + "_u" + std::to_string(nunif) + "h" + std::to_string(nref) +
                                "p" + std::to_string(pmax);
    const std::string linha(100, '=');
    std::cout << "\n" << linha << "\n GIM | " << titulo << " | alvo Gamma = " << alvo << "\n" << linha << std::endl;
    std::vector<SummaryRun> runs;
    const auto t0 = Clock::now();
    for (MaterialKind kind : {MaterialKind::UP, MaterialKind::Unsat}) {
        std::cout << "\n--- material " << Name(kind) << " ---" << std::endl;
        SlopeAnalysis a(d, kind, rebaix ? WaterModel::Drawdown : WaterModel::Dry, FactorMethod::GIMArcLength);
        a.SetRefinement({nunif, nref, pmax, 0.1});
        a.SetTarget(alvo);
        a.SetOutputPrefix(prefixo); // arquivos <prefixo>_<up|unsat>... para o figuras_gim.py
        runs.push_back({Name(kind), alvo, a.Run(m0)});
    }
    PrintSummary(titulo, runs, Seconds(t0));
    return 0;
}
