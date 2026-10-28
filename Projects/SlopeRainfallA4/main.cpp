// Caso A4 (contraexemplo ao "NA inicial é secundário" de Rahardjo et al. 2007) com sucção inicial SEM limite.
// Tudo explícito e fixo no código: malha geométrica, malhas computacionais (u, p_w, p_a e multifísica), material,
// matriz estrutural, solver, análise, condições de contorno, transiente com chuva, fator de segurança (SRM) e saída
// VTK para o ParaView (abra rain_A4.vtk.series: a série de arquivos tem o tempo em horas).
//
// Modelo do artigo: percolação desacoplada (rigid skeleton nas equações de massa) e ar passivo (p_a = 0).
// Solo f50, alpha = 45 graus, Hs = 10 m, NA a Hw = 10 m abaixo do pé, chuva de 36 mm/h por 24 h, 48 h simuladas.
#include "../SlopeRainfall/UnsatSlope.h"

#include <chrono>
#include <iomanip>
#include <limits>

#ifndef SLOPERAINFALL_MESH_DIR
#define SLOPERAINFALL_MESH_DIR "."
#endif

int main() {
    // ------------------------------------------------------------------ dados do caso (fixos)
    const REAL HOUR = 3600.;
    const REAL Hs = 10., alpha = 45., Hw = 10.;          // geometria (Fig. 2 do artigo) e NA abaixo do pé
    const REAL Irain = 36.;                              // mm/h
    const REAL tRain = 24., tTotal = 48.;                // h
    const REAL smax = std::numeric_limits<REAL>::max();  // sucção inicial sem limite (o artigo usa 75 kPa)
    const bool acoplado = false, passivo = true;         // hipóteses do artigo
    const std::vector<REAL> tFS = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 36, 48};   // instantes do FS (h)
    const std::string meshFile = std::string(SLOPERAINFALL_MESH_DIR) + "/talude_H10_a45.msh";
    const std::string vtkName = "rain_A4.vtk";
    const int nthreads = 0;                              // threads da montagem (0: serial)

    const REAL xb = 3. * Hs, xc = xb + Hs / std::tan(alpha * M_PI / 180.), L = xc + 3. * Hs;
    auto waterTable = [&](REAL x) { return (2. * Hs - Hw) - std::tan(7. * M_PI / 180.) * (x - xc); };

    // parâmetros do solo f50 (Tabela 1 do artigo) e hipóteses de talude3f
    UnsatParams par;
    par.hydro.fA = 50.;          // Fredlund-Xing a (kPa), n = m = 1, k_rw = S^4, Brooks-Corey para o ar
    par.ks = 1.e-5;              // m/s
    // E = 20 MPa, nu = 0,3, c' = 10 kPa, phi' = psi = 26 graus, gamma = 20 kN/m3, chi = 1 (padrões de UnsatParams)

    // ------------------------------------------------------------------ 1. malha geométrica (gmsh 2.2, triângulos P1)
    // materiais: 1 domínio; -1 base, -2 lateral esquerda, -3 lateral direita, -4 crista, -5 face, -6 terreno do pé
    TPZGeoMesh *gmesh = ReadMsh(meshFile);
    std::cout << "malha geométrica: " << gmesh->NElements() << " elementos, " << gmesh->NNodes() << " nós\n";

    // ------------------------------------------------------------------ 2. malhas computacionais
    // espaços atômicos H1: deslocamento P2 (2 estados), p_w P1, p_a P1, com material nulo no domínio e contornos
    auto atomic = [&](int nstate, int order) {
        auto *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDimModel(2);
        cmesh->SetDefaultOrder(order);
        cmesh->SetAllCreateFunctionsContinuous();
        auto *nullmat = new TPZNullMaterial<STATE>(1, 2, nstate);
        cmesh->InsertMaterialObject(nullmat);
        TPZFMatrix<STATE> v1(nstate, nstate, 0.);
        TPZManVector<STATE, 2> v2(nstate, 0.);
        for (int id = -6; id <= -1; ++id) cmesh->InsertMaterialObject(nullmat->CreateBC(nullmat, id, 0, v1, v2));
        cmesh->AutoBuild();
        return cmesh;
    };
    TPZManVector<TPZCompMesh *, 3> meshvec = {atomic(2, 2), atomic(1, 1), atomic(1, 1)};

    // material u - p_w - p_a (Mohr-Coulomb na tensão de Bishop)
    auto *mat = new TUnsat(1);
    mat->SetPlasticModel(MohrCoulomb(par));
    mat->SetHydraulics(par.hydro);
    mat->SetPorosity(par.n);
    mat->SetBiot(1., 0.);                                // alpha = 1, grãos incompressíveis
    mat->SetWaterBulkModulus(par.Kw);
    mat->SetConductivity(par.ks, par.gw, par.muRatio);
    mat->SetAir(par.ga, par.patm);
    mat->SetUnitWeight(par.gamma);
    mat->SetBishop(TUnsat::EChiConstant, par.chi);       // chi = tan(phi_b)/tan(phi') = 1
    mat->SetLumpedStorage(true);
    mat->SetIntegrationOrder(2);                         // 3 pontos por triângulo
    mat->SetHydroMechanicalCoupling(acoplado);

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);
    {
        TPZFMatrix<STATE> v1(2, 2, 0.);
        TPZManVector<STATE, 2> v2(2, 0.);
        std::set<int> bcids;
        for (int id = -6; id <= -1; ++id) {   // Dirichlet impostos pela análise (filtro de equações)
            mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, v1, v2));
            bcids.insert(id);
        }
        TPZManVector<int, 3> active(3, 1);
        mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {1}, bcids);
    }
    std::cout << "malha multifísica: " << mphys->NEquations() << " equações\n";

    // ------------------------------------------------------------------ 3. análise, matriz estrutural e solver
    UnsatAnalysis an(mphys, nthreads);                   // Newton com o filtro de equações (Dirichlet eliminados)
#ifdef PZ_USING_MKL
    TPZSpStructMatrix<STATE> smat(mphys);                // CSR + Pardiso
#else
    TPZSkylineNSymStructMatrix<STATE> smat(mphys);
#endif
    smat.SetNumThreads(nthreads);
    an.SetStructuralMatrix(smat);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    an.SetSolver(solver);
    an.SetTolerances(1.e-7, 1.e-6, 40, 25);

    // ------------------------------------------------------------------ 4. condições de contorno
    std::map<int64_t, REAL> fixed = an.DisplacementConstraints();   // base fixa, laterais em rolete
    std::set<int64_t> surf, left, right;
    std::map<int64_t, REAL> supply;                                  // área de chuva por vértice (projeção horizontal)
    for (int64_t i = 0; i < gmesh->NElements(); ++i) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->Dimension() != 1) continue;
        const int id = gel->MaterialId();
        const int64_t a = gel->NodeIndex(0), b = gel->NodeIndex(1);
        if (id == ECrest || id == EFace || id == EToe) {
            surf.insert(a);
            surf.insert(b);
            const REAL dx = std::fabs(gel->NodePtr(1)->Coord(0) - gel->NodePtr(0)->Coord(0));
            supply[a] += dx / 2.;
            supply[b] += dx / 2.;
        } else if (id == ELeft) left.insert(a), left.insert(b);
        else if (id == ERight) right.insert(a), right.insert(b);
    }
    const REAL hl = waterTable(0.), hr = waterTable(L);
    for (auto n : left)
        if (gmesh->NodeVec()[n].Coord(1) <= hl + 1.e-6) fixed[an.WaterEq().at(n)] = par.gw * (hl - gmesh->NodeVec()[n].Coord(1));
    for (auto n : right)
        if (gmesh->NodeVec()[n].Coord(1) <= hr + 1.e-6) fixed[an.WaterEq().at(n)] = par.gw * (hr - gmesh->NodeVec()[n].Coord(1));
    if (passivo)
        for (auto &e : an.AirEq()) fixed[e.second] = 0.;
    else
        for (auto n : surf) fixed[an.AirEq().at(n)] = 0.;

    // ------------------------------------------------------------------ 5. estado inicial
    // pressões hidrostáticas a partir do NA, sem limite de sucção; p_a = 0
    TPZFMatrix<STATE> &X = an.X();
    REAL suctionMax = 0.;
    for (auto &e : an.WaterEq()) {
        const TPZGeoNode &nd = gmesh->NodeVec()[e.first];
        X(e.second, 0) = std::max(par.gw * (waterTable(nd.Coord(0)) - nd.Coord(1)), -smax);
        suctionMax = std::max(suctionMax, -X(e.second, 0));
    }
    for (auto &e : an.AirEq()) X(e.second, 0) = 0.;
    for (auto &f : fixed) X(f.first, 0) = f.second;
    an.Restore();
    std::cout << "sucção inicial máxima: " << suctionMax << " kPa\n";

    auto frozenFixed = [&](const TPZFMatrix<STATE> &x) {   // u restrito e todas as pressões nos valores de x
        std::map<int64_t, REAL> fx = an.DisplacementConstraints();
        for (int64_t i = 0; i < an.NEq(); ++i)
            if (an.IsPressure(i)) fx[i] = x.GetVal(i, 0);
        return fx;
    };
    auto geostatic = [&](REAL F) {   // peso em 4 incrementos com as pressões congeladas, resistência c/F, tan(phi)/F
        mat->PlasticModel().SetStrengthReductionFactor(F);
        const auto fx = frozenFixed(an.X());
        for (REAL lam : {0.25, 0.5, 0.75, 1.0}) {
            mat->SetGravityFactor(lam);
            if (!an.Solve(0., fx, {}, true).conv) return false;
            an.Accept();
        }
        mat->SetGravityFactor(1.);
        return true;
    };
    // desacoplado: a mecânica do transiente não muda as pressões; resistência aumentada (F = 0,25) só para manter o
    // equilíbrio. A estabilidade vem das análises de FS.
    const REAL Ftransient = acoplado ? 1. : 0.25;
    if (!geostatic(Ftransient)) {
        std::cout << "estado geostático não convergiu\n";
        return 1;
    }
    an.ResetDisplacements();
    an.Accept();

    // ------------------------------------------------------------------ 6. fator de segurança (SRM)
    const REAL F0 = 0.6;
    auto factorOfSafety = [&](REAL prev) {
        const TPZAdmChunkVector<TPZUnsaturatedMem> saved(*mat->GetMemory());
        const TPZFMatrix<STATE> xsaved(an.X());
        // estado sem tensões com as pressões do instante, geostático com a resistência aumentada (F0)
        auto mem = mat->GetMemory();
        const TPZUnsaturatedMem def = mat->GetDefaultMemory();
        for (int64_t i = 0; i < mem->NElements(); ++i) (*mem)[i] = def;
        for (int64_t i = 0; i < an.NEq(); ++i)
            if (!an.IsPressure(i)) an.X()(i, 0) = 0.;
        an.Restore();
        REAL fs = -1.;
        if (geostatic(F0)) {
            const auto fx = frozenFixed(an.X());
            std::vector<REAL> xs = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0};
            if (prev > 0.) {
                xs.clear();
                for (REAL d : {-0.25, -0.08, -0.02, 0.03, 0.1, 0.4})
                    if (prev + d > F0) xs.push_back(prev + d);
            }
            // busca do colapso: segue xs e, na primeira falha, bissecção até 0,01
            REAL xok = F0, xfail = 1.e30;
            while (true) {
                REAL F;
                if (xfail < 1.e29 && xfail - xok > 0.01) F = 0.5 * (xok + xfail);
                else if (xfail < 1.e29) F = xfail;
                else {
                    F = xok + 0.5;
                    for (REAL v : xs)
                        if (v > xok + 1.e-12) { F = v; break; }
                }
                if (F > 20.) break;
                mat->PlasticModel().SetStrengthReductionFactor(F);
                const bool ok = an.Solve(0., fx, {}, true).conv;
                if (ok) {
                    an.Accept();
                    if (F == xfail) xfail = 1.e30;
                    xok = F;
                } else {
                    an.Restore();
                    if (xfail < 1.e29 && F == xfail) break;
                    xfail = F;
                }
            }
            fs = xok;
        }
        *mat->GetMemory() = saved;   // de volta ao transiente
        an.X() = xsaved;
        an.Restore();
        mat->SetGravityFactor(1.);
        mat->PlasticModel().SetStrengthReductionFactor(Ftransient);
        return fs;
    };

    // ------------------------------------------------------------------ 7. pós-processamento (ParaView)
    TPZStack<std::string> scal, vec;
    for (auto s : {"WaterPressure", "AirPressure", "Suction", "Saturation", "BishopPressure", "RelativePermeabilityWater"})
        scal.Push(s);
    vec.Push("Displacement");
    an.DefineGraphMesh(2, scal, vec, vtkName);
    std::vector<std::pair<int, REAL>> frames;   // (passo do VTK, tempo em h)
    auto writeVTK = [&](REAL th) {
        const int k = int(frames.size());
        an.Restore();
        an.SetStep(k);
        an.PostProcess(0, 2);
        frames.push_back({k, th});
        std::ofstream ser("rain_A4.vtk.series");   // série com o tempo físico para o ParaView
        ser << "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n";
        for (size_t i = 0; i < frames.size(); ++i)
            ser << "    { \"name\": \"rain_A4.scal_vec." << frames[i].first << ".vtk\", \"time\": " << frames[i].second
                << " }" << (i + 1 < frames.size() ? "," : "") << "\n";
        ser << "  ]\n}\n";
    };
    std::ofstream fsout("fs_A4.csv");
    fsout << "t_h,fs_srm,ponded\n";
    REAL fsPrev = -1.;
    std::set<int64_t> ponded;
    auto report = [&](REAL th) {
        const auto t0 = std::chrono::steady_clock::now();
        fsPrev = factorOfSafety(fsPrev);
        fsout << th << "," << fsPrev << "," << ponded.size() << std::endl;
        std::cout << "t = " << std::setw(5) << th << " h   FS(SRM) = " << fsPrev << "   vértices alagados = "
                  << ponded.size() << "   ("
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)" << std::endl;
    };
    writeVTK(0.);
    report(0.);

    // ------------------------------------------------------------------ 8. transiente: chuva unilateral
    // p_w <= 0 na superfície, infiltração <= chuva, complementaridade (conjunto ativo de vértices alagados)
    std::vector<int64_t> atm(surf.begin(), surf.end());
    auto unilateralStep = [&](REAL dt, REAL rate) {
        std::set<int64_t> pond = ponded;
        bool conv = false;
        for (int pass = 0; pass < 12; ++pass) {
            std::map<int64_t, REAL> fx(fixed), fq;
            for (auto n : pond) fx[an.WaterEq().at(n)] = 0.;
            for (auto n : atm)
                if (!pond.count(n) && rate > 0.) fq[an.WaterEq().at(n)] = rate * supply[n];
            conv = an.Solve(dt, fx, fq, false).conv;
            if (!conv && rate > 0. && pond.size() < atm.size()) {   // chuva intensa em solo seco: começa alagado
                an.Restore();
                pond.insert(atm.begin(), atm.end());
                continue;
            }
            if (!conv) break;
            const TPZFMatrix<STATE> &x = mphys->Solution();
            std::set<int64_t> add, drop;
            for (auto n : atm)
                if (!pond.count(n) && x.GetVal(an.WaterEq().at(n), 0) > 1.e-6) add.insert(n);
            TPZFMatrix<STATE> Rw = an.FullResidual();   // Rw/dt = vazão necessária nos vértices alagados
            for (auto n : pond)
                if (Rw.GetVal(an.WaterEq().at(n), 0) / dt > rate * supply[n] * (1. + 1.e-9) + 1.e-14) drop.insert(n);
            if (add.empty() && drop.empty()) break;
            for (auto n : add) pond.insert(n);
            for (auto n : drop) pond.erase(n);
        }
        if (!conv) {
            an.Restore();
            return false;
        }
        an.Accept();
        ponded = pond;
        return true;
    };
    std::function<bool(REAL, REAL, int)> advance = [&](REAL dt, REAL rate, int depth) {   // bissecção do passo
        if (unilateralStep(dt, rate)) return true;
        if (depth >= 5) return false;
        return advance(dt / 2., rate, depth + 1) && advance(dt / 2., rate, depth + 1);
    };

    REAL t = 0.;   // h
    while (t < tTotal - 1.e-9) {
        const bool raining = t < tRain - 1.e-9;
        REAL dt = raining ? 0.25 : (t < tRain + 24. - 1.e-9 ? 1. : 4.);   // 15 min na chuva, 1 h, depois 4 h
        dt = std::min(dt, tTotal - t);
        if (raining) dt = std::min(dt, tRain - t);
        const REAL rate = raining ? Irain / 1000. / HOUR : 0.;          // m/s
        if (!advance(dt * HOUR, rate, 0)) {
            std::cout << "passo sem convergência em t = " << t << " h\n";
            break;
        }
        t += dt;
        if (std::fabs(t - std::round(t)) < 1.e-6) writeVTK(t);          // VTK a cada hora
        for (REAL tf : tFS)
            if (tf > 0. && std::fabs(t - tf) < 1.e-6) report(t);
    }
    std::cout << "fim: VTK em rain_A4.scal_vec.*.vtk (abra rain_A4.vtk.series no ParaView), FS em fs_A4.csv\n";

    // ------------------------------------------------------------------ limpeza
    return 0;   // as malhas são liberadas no fim do processo (a análise referencia a malha multifísica)
}
