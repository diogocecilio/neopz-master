// Rapid drawdown of the reservoir in front of the slope of SlopeMohrCoulomb (saturated soil):
//  1. coupled Biot u-p poro-elastoplastic analysis with the Mohr-Coulomb model of the paper
//     (TPZMatPoroElastoPlasticUP<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>>): steady state with the reservoir
//     at the crest, instantaneous (undrained) drawdown to the toe, consolidation (transient seepage);
//  2. at selected instants the factor of safety is computed with the procedure of SlopeMohrCoulomb
//     (SlopeAnalysis.h: gravity increase and strength reduction, refinement of the plastic zone), drained,
//     with the pore pressure frozen and applied as the seepage body force b = gamma_sat g - grad p+
//     (p+ = max(p, 0): suction neglected).
//
// Usage: SlopeDrawdown [nref=<n>] [ref=<uniform refinements of the u-p and FS meshes, default 1>]
//                      [L=<reservoir level after the drawdown in [30, 40], default 30 (toe)>]
#include "../SlopeMohrCoulomb/SlopeAnalysis.h"
#include "../SlopeMohrCoulomb/SlopeModel.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZPoroElastoPlasticUPAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzgeoelbc.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <map>
#include <string>
#include <vector>

using TPoroMC = TPZMatPoroElastoPlasticUP<TMCVoigt, TPZElastoPlasticMem>;
using TUPAnalysis = TPZPoroElastoPlasticUPAnalysis;

/// Water: unit weight, hydraulic conductivity (time unit: day) and reservoir levels (crest y = 40, toe y = 30)
struct Water {
    REAL gammaw = 10., kh = 1.e-2, H0 = 40., H1 = 30.;
};

/// Toe ground (-3) and face (-6) need a pore pressure and a water load: twin boundary elements -13 and -16
// void AddWaterLoadElements2(TPZGeoMesh *gmesh) {
//     const int64_t nel = gmesh->NElements();
//     for (int64_t i = 0; i < nel; i++) {
//         TPZGeoEl *gel = gmesh->Element(i);
//         if (gel && !gel->HasSubElement() && (gel->MaterialId() == -3 || gel->MaterialId() == -6))
//             TPZGeoElBC(gel, gel->NSides() - 1, gel->MaterialId() - 10);
//     }
// }
void AddWaterLoadElements(TPZGeoMesh *gmesh) {
    const int64_t nel = gmesh->NElements();
    for (int64_t i = 0; i < nel; ++i) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == -3 || id == -4 || id == -6)
            TPZGeoElBC(gel, gel->NSides() - 1, id - 10);
        // -3 -> -13; -4 -> -14; -6 -> -16.
    }
}

/// Atomic H1 space of the u-p mesh with a null material: nstate = 2 displacement, 1 pore pressure
TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids) {
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(1, 2, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

/// Biot u-p mesh of the slope (plane strain, Taylor-Hood: u P2 and p P1, porder = 2).
/// Base fixed, sides on rollers, base and sides impermeable. Reservoir level H = H0 - lambda (H0 - H1), lambda =
/// load factor of the step (0: crest, 1: drawn down; bisection interpolates it). Ground surface (toe -3, crest -4,
/// face -6): p = gamma_w (H - y)+, i.e. drained above the water (the phreatic surface stays at the crest, as in
/// the rapid drawdown of Ceron et al. 2025); water load t = -gamma_w (H - y)+ n on the twins -13 and -16.
template <class TPlastic>
TPZMultiphysicsCompMesh *CreateCompMesh2(TPZGeoMesh *gmesh, int porder, const TPlastic &model, const Soil &s,
                                        const Water &w) {
    using B = TPZMatPoroElastoPlasticUPBase;
    if (porder != 2) DebugStop(); // PoreField and the nodal Dirichlet values assume a P1 pressure
    const std::set<int> bcids = {-1, -2, -3, -4, -5, -6, -13, -16};
    TPZManVector<TPZCompMesh *, 2> meshvec = {CreateAtomicMesh(gmesh, 2, porder, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids)};
    auto *mat = new TPZMatPoroElastoPlasticUP<TPlastic, TPZElastoPlasticMem>(1, B::EPlaneStrain);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., 0.);                  // incompressible grains and water
    mat->SetPermeability(w.kh / w.gammaw); // mobility
    TPZManVector<REAL, 3> b = {0., -s.gamma, 0.}, gw = {0., -w.gammaw, 0.};
    mat->SetBodyForce(b); // saturated unit weight
    mat->SetFluidWeight(gw);
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);
    auto pw = [mat, w](const TPZVec<REAL> &x) {
        const REAL H = w.H0 - mat->LoadFactor() * (w.H0 - w.H1);
        return w.gammaw * std::max<REAL>(H - x[1], 0.);
    };
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    mphys->InsertMaterialObject(mat->CreateBC(mat, -1, B::EDirichletU, val1, val2));
    for (int id : {-3, -4, -6}) {
        auto *bc = mat->CreateBC(mat, id, B::EDirichletP, val1, val2);
        bc->SetForcingFunctionBC([pw](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) { v[0] = pw(x); });
        mphys->InsertMaterialObject(bc);
    }
    const REAL normal[2][2] = {{0., 1.}, {M_SQRT1_2, M_SQRT1_2}}; // outward normals of -13 (toe) and -16 (face)
    for (int i = 0; i < 2; i++) {
        auto *bc = mat->CreateBC(mat, i ? -16 : -13, B::ENeumannUFixed, val1, val2);
        const REAL nx = normal[i][0], ny = normal[i][1];
        bc->SetForcingFunctionBC([pw, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) {
            v[0] = -pw(x) * nx;
            v[1] = -pw(x) * ny;
        });
        mphys->InsertMaterialObject(bc);
    }
    val1(0, 0) = 1.; // u_x = 0
    for (int id : {-2, -5}) mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    TPZManVector<int, 2> active(2, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {1}, bcids);
    return mphys;
}

/// Pore pressure of the u-p analysis (linear on each triangle of the u-p mesh) at any point of the slope
class PoreField { // the FS mesh must be the u-p mesh or a refinement of it (grad p constant per element)
    struct Tri {
        REAL x0[2], inv[4], p0, dp[2]; // x = x0 + J xi, inv = J^-1, p = p0 + dp . xi
    };
    std::vector<Tri> fTri;

public:
    PoreField() = default; // dry slope

    /// linear interpolation of the vertex values p(node) on the volume triangles of gmesh
    PoreField(TPZGeoMesh *gmesh, const std::function<REAL(int64_t)> &pnode) {
        for (int64_t i = 0; i < gmesh->NElements(); i++) {
            TPZGeoEl *gel = gmesh->Element(i);
            if (!gel || gel->HasSubElement() || gel->MaterialId() != 1) continue;
            REAL x[3][2], p[3];
            for (int k = 0; k < 3; k++) {
                x[k][0] = gel->NodePtr(k)->Coord(0);
                x[k][1] = gel->NodePtr(k)->Coord(1);
                p[k] = pnode(gel->NodeIndex(k));
            }
            const REAL a = x[1][0] - x[0][0], b = x[2][0] - x[0][0], c = x[1][1] - x[0][1], d = x[2][1] - x[0][1];
            const REAL det = a * d - b * c;
            fTri.push_back({{x[0][0], x[0][1]}, {d / det, -b / det, -c / det, a / det}, p[0], {p[1] - p[0], p[2] - p[0]}});
        }
    }

    /// p and grad p at x (false outside the mesh)
    bool Evaluate(const TPZVec<REAL> &x, REAL &p, REAL grad[2]) const {
        thread_local size_t last = 0; // consecutive integration points lie in the same triangle
        const size_t n = fTri.size();
        for (size_t k = 0; k < n; k++) {
            const size_t i = (last + k) % n;
            const Tri &t = fTri[i];
            const REAL dx = x[0] - t.x0[0], dy = x[1] - t.x0[1];
            const REAL xi = t.inv[0] * dx + t.inv[1] * dy, eta = t.inv[2] * dx + t.inv[3] * dy;
            if (xi < -1.e-10 || eta < -1.e-10 || xi + eta > 1. + 1.e-10) continue;
            last = i;
            p = t.p0 + t.dp[0] * xi + t.dp[1] * eta;
            grad[0] = t.inv[0] * t.dp[0] + t.inv[2] * t.dp[1]; // J^-T grad_xi p
            grad[1] = t.inv[1] * t.dp[0] + t.inv[3] * t.dp[1];
            return true;
        }
        return false;
    }

    /// grad p+ = grad max(p, 0) (suction neglected); zero outside the u-p mesh and for the dry slope
    void GradPositive(const TPZVec<REAL> &x, REAL grad[2]) const {
        REAL p;
        if (!Evaluate(x, p, grad) || p <= 0.) grad[0] = grad[1] = 0.;
    }
};

/// Seepage body force of the stability analysis, b = lambda (gamma_sat g - grad p+): lambda is the gravity factor
/// that SlopeAnalysis sets through the body force of the material (0, -lambda gamma_sat, 0)
template <class TPlastic>
void SetSeepageForce(TPZCompMesh *cmesh, const PoreField &pf, REAL gamma) {
    auto *mat = dynamic_cast<TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem> *>(cmesh->FindMaterial(1));
    if (!mat) DebugStop();
    mat->SetForcingFunction(
        [mat, &pf, gamma](const TPZVec<REAL> &x, TPZVec<STATE> &f) {
            const REAL lambda = -mat->GetBodyForce()[1] / gamma;
            REAL g[2];
            pf.GradPositive(x, g);
            f[0] = -lambda * g[0];
            f[1] = -lambda * (gamma + g[1]);
            f[2] = 0.;
        },
        0);
}

struct State {
    std::string name;
    REAL t; ///< time after the drawdown (day); < 0 before
    PoreField pf;
};

/// FS by gravity increase and by strength reduction (as in SlopeMohrCoulomb) for a frozen pore pressure
template <class TPlastic>
std::pair<REAL, REAL> FactorOfSafety(const TPlastic &model, const Soil &s, const State &st, int ref, int nref) {
    TPZGeoMesh *gmesh = TriGMesh(ref); // same level as the u-p mesh: grad p+ is constant in each element
    TPZCompMesh *cmesh = CreateCMesh(gmesh, 2, model, s);
    SetSeepageForce<TPlastic>(cmesh, st.pf, s.gamma);
    SlopeAnalysis<TPlastic> slope(cmesh);
    REAL fsGI = 0., fsSRM = 0.;
    for (int k = 0;; k++) {
        fsGI = slope.GravityIncrease();
        slope.MarkPlasticZone(0.1);
        fsSRM = slope.StrengthReduction();
        slope.MarkPlasticZone(0.1);
        std::cout << "[" << st.name << "] refinement " << k << ": " << cmesh->NEquations() << " equations, FS GI "
                  << fsGI << ", FS SRM " << fsSRM << "\n";
        if (k >= nref) break;
        slope.Refine();
    }
    slope.PostPlasticity("drawdown_" + st.name + ".vtk");
    delete cmesh;
    delete gmesh;
    return {fsGI, fsSRM};
}

/// u-p analysis: pore pressure before the drawdown, right after it (undrained) and during the consolidation
std::vector<State> PorePressureStates(const TMCVoigt &model, const Soil &s, const Water &w, REAL Tc, int ref) {
    TPZGeoMesh *gmesh = TriGMesh(ref);
    AddWaterLoadElements(gmesh);
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh2(gmesh, 2, model, s, w); // no renumbering: p after u (ELU)
    std::vector<State> states;
    {
        TUPAnalysis an(mphys, dynamic_cast<TPoroMC *>(mphys->FindMaterial(1)));
        TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
        an.SetStructuralMatrix(skyl);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        an.SetSolver(step);
        an.SetPredictor(false);
        an.SetVerbose(1);
        TPZStack<std::string> scal, vec;
        scal.Push("PorePressure");
        vec.Push("Displacement");
        an.DefineGraphMesh(2, scal, vec, "drawdown_up.vtk");
        auto record = [&](const std::string &name, REAL t) {
            states.push_back({name, t, PoreField(gmesh, [&an](int64_t n) { return an.NodalValue(n, 1, 0); })});
            an.SetStep(int(states.size()) - 1);
            an.PostProcess(1, 2);
        };
        // reservoir at the crest: gravity in one drained step (steady state); rapid drawdown at t = 0: undrained
        // step (Dt = 0, lambda 0 -> 1); consolidation: T = cv t / H^2 = 1e-3 ... 10 (4 steps per decade) and steady
        // seepage with the phreatic surface kept at the crest (the rapid-drawdown field of Ceron et al.)
        using TS = TUPAnalysis::TLoadState;
        bool ok = an.Run({TS(0., 0., 0.)}, nullptr, TS(-1.e4 * Tc, 0., 0.));
        if (ok) record("crest", -1.);
        ok = ok && an.Run({TS(0., 1., 0.)}, nullptr, TS(0., 0., 0.));
        if (ok) record("T0", 0.);
        const std::map<int, std::string> out = {{8, "T0.1"}, {12, "T1"}, {16, "T10"}, {17, "steady"}};
        std::vector<TS> steps;
        for (int j = 0; j <= 16; j++) steps.emplace_back(std::pow(10., -3. + 0.25 * j) * Tc, 1., 0.);
        steps.emplace_back(1.e4 * Tc, 1., 0.);
        ok = ok && an.Run(steps, [&](int k, const TS &st) { if (out.count(k - 1)) record(out.at(k - 1), st.fTime); },
                          TS(0., 1., 0.));
        if (!ok) std::cout << "u-p analysis stopped\n";
    }
    // reference: drained after the drawdown without recharge (water table at the reservoir level, hydrostatic)
    states.push_back({"drained", -1., PoreField(gmesh, [&](int64_t n) {
                          return w.gammaw * std::max<REAL>(w.H1 - gmesh->NodeVec()[n].Coord(1), 0.); })});
    TPZManVector<TPZCompMesh *, 2> atomic(mphys->MeshVector());
    delete mphys;
    for (auto *m : atomic) delete m;
    delete gmesh;
    return states;
}
/*
int main(int argc, char *argv[]) {
    int nref = 2, ref = 1;
    Soil s;
    s.nu = 0.3; // drained Poisson ratio of the skeleton (the undrained response comes from the coupling)
    Water w;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "nref=", 5)) nref = std:: max(0, atoi(argv[i] + 5));
        else if (!strncmp(argv[i], "ref=", 4)) ref = std::max(0, atoi(argv[i] + 4));
        else if (!strncmp(argv[i], "L=", 2)) w.H1 = atof(argv[i] + 2);
    }
    if (w.H1 < 30. || w.H1 > w.H0) { // no water load on the crest, no reservoir below the toe
        std::cout << "L must be in [30, 40]\n";
        return 1;
    }
    const TMCVoigt model = ModelVoigt(s);
    const REAL M = s.E * (1. - s.nu) / ((1. + s.nu) * (1. - 2. * s.nu)); // oedometric modulus
    const REAL Tc = 100. / (w.kh / w.gammaw * M);                     // t = T H^2 / cv, H = 10 m
    std::vector<State> states = {{"dry", -1., PoreField()}};
    for (auto &st : PorePressureStates(model, s, w, Tc, ref)) states.push_back(st);
    std::vector<std::pair<REAL, REAL>> fs;
    for (auto &st : states) fs.push_back(FactorOfSafety(model, s, st, ref, nref));
    std::cout << "\nstate  t (day)  T = cv t / H^2  FS(gravity increase)  FS(strength reduction)\n";
    for (size_t i = 0; i < states.size(); i++) {
        std::cout << "  " << states[i].name << "  ";
        if (states[i].t < 0.) std::cout << "-  -";
        else std::cout << states[i].t << "  " << states[i].t / Tc;
        std::cout << "  " << fs[i].first << "  " << fs[i].second << "\n";
    }
    return 0;
}*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZVTKGeoMesh.h"
#include <algorithm>
#include <fstream>
#include <string>



template <class TPlastic>
TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int porder,
                                      const TPlastic &model, const Soil &s,
                                      const Water &w) {
    using B = TPZMatPoroElastoPlasticUPBase;
    if (porder != 2) DebugStop();

    const std::set<int> bcids = {-1, -2, -3, -4, -5, -6, -13, -14, -16};
    TPZManVector<TPZCompMesh *, 2> meshvec = {
        CreateAtomicMesh(gmesh, 2, porder, bcids),
        CreateAtomicMesh(gmesh, 1, porder - 1, bcids)
    };

    auto *mat = new TPZMatPoroElastoPlasticUP<TPlastic, TPZElastoPlasticMem>(
        1, B::EPlaneStrain);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., 0.);
    mat->SetPermeability(w.kh / w.gammaw);

    TPZManVector<REAL, 3> body = {0., -s.gamma, 0.};
    TPZManVector<REAL, 3> fluid = {0., -w.gammaw, 0.};
    mat->SetBodyForce(body);
    mat->SetFluidWeight(fluid);

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);

    auto pw = [mat, w](const TPZVec<REAL> &x) {
        const REAL H = w.H0 - mat->LoadFactor() * (w.H0 - w.H1);
        return w.gammaw * std::max<REAL>(H - x[1], 0.);
    };

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    mphys->InsertMaterialObject(
        mat->CreateBC(mat, -1, B::EDirichletU, val1, val2));

    // Condições hidráulicas, incluindo o topo -4.
    for (int id : {-3, -4, -6}) {
        auto *bc = mat->CreateBC(mat, id, B::EDirichletP, val1, val2);
        bc->SetForcingFunctionBC(
            [pw](const TPZVec<REAL> &x, TPZVec<STATE> &v,
                 TPZFMatrix<STATE> &) { v[0] = pw(x); });
        mphys->InsertMaterialObject(bc);
    }

    // Trações exercidas pela água.
    const REAL normal[2][2] = {{0., 1.}, {M_SQRT1_2, M_SQRT1_2}};
    for (int i = 0; i < 2; ++i) {
        auto *bc = mat->CreateBC(
            mat, i ? -16 : -13, B::ENeumannUFixed, val1, val2);
        const REAL nx = normal[i][0], ny = normal[i][1];
        bc->SetForcingFunctionBC(
            [pw, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v,
                         TPZFMatrix<STATE> &) {
                v[0] = -pw(x) * nx;
                v[1] = -pw(x) * ny;
            });
        mphys->InsertMaterialObject(bc);
    }

    // Carga vertical no topo: t_y = -20 * lambda.
    TPZManVector<STATE, 2> topLoad = {0., -20.};
    mphys->InsertMaterialObject(mat->CreateBC(mat, -14, B::ENeumannU, val1, topLoad));

    // Apoios laterais.
    val1(0, 0) = 1.;
    for (int id : {-2, -5})
        mphys->InsertMaterialObject(
            mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));

    TPZManVector<int, 2> active(2, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {1}, bcids);
    return mphys;
}
// 1. Colocar antes de Consolidation10Years.
void PostPlasticityUP(TPZMultiphysicsCompMesh *mphys,
                      const std::string &filename) {
    {
        TPZPostProcAnalysis post;
        post.SetCompMesh(mphys);

        TPZStack<std::string> scal, vec, variables;
        scal.Push("PorePressure");
        scal.Push("MeanEffectiveStress");
        scal.Push("DeviatoricStress");
        scal.Push("EffectiveStressXX");
        scal.Push("EffectiveStressYY");
        scal.Push("EffectiveStressZZ");
        scal.Push("EffectiveStressXY");
        scal.Push("VolumetricStrain");
        scal.Push("PlasticType");
        vec.Push("Displacement");

        for (const auto &name : scal) variables.Push(name);
        for (const auto &name : vec) variables.Push(name);

        TPZManVector<int, 1> materialIds(1, 1);
        TPZFStructMatrix<STATE> matrix(post.Mesh());
        post.SetStructuralMatrix(matrix);
        post.SetPostProcessVariables(materialIds, variables);

        // Projeta os campos da solução e da memória dos pontos de integração.
        post.TransferSolution();
        post.DefineGraphMesh(2, scal, vec, filename);
        post.PostProcess(0);
    }

    // Restaura as referências da malha usada na consolidação.
    mphys->LoadReferences();

    // // Exporta também a malha geométrica, como em slope.PostPlasticity.
    // std::ofstream meshFile(filename);
    // TPZVTKGeoMesh::PrintGMeshVTK(mphys->Reference(), meshFile, true);
}


#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Manter PostPlasticityUP e as demais funções do arquivo.
bool Consolidation10Years(const TMCVoigt &model, const Soil &s,
                          const Water &w, REAL Tc, int ref) {
    TPZGeoMesh *gmesh = TriGMesh(ref);
    AddWaterLoadElements(gmesh);
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, 2, model, s, w);
    bool ok = false;

    try {
        TUPAnalysis an(mphys, dynamic_cast<TPoroMC *>(mphys->FindMaterial(1)));
        TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
        an.SetStructuralMatrix(skyl);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        an.SetSolver(solver);
        an.SetPredictor(false);
        an.SetVerbose(1);

        using TS = TUPAnalysis::TLoadState;
        int outputStep = 0;
        std::vector<std::pair<std::string, REAL>> vtkFrames;

        auto writeSeries = [&]() {
            std::ofstream file("drawdown.vtk.series");
            if (!file)
                throw std::runtime_error("Erro ao criar drawdown.vtk.series");

            file.imbue(std::locale::classic());
            file << std::setprecision(17);
            file << "{\n"
                 << "  \"file-series-version\": \"1.0\",\n"
                 << "  \"files\": [\n";

            for (size_t i = 0; i < vtkFrames.size(); ++i) {
                file << "    { \"name\": \"" << vtkFrames[i].first
                     << "\", \"time\": " << vtkFrames[i].second << " }";
                if (i + 1 < vtkFrames.size()) file << ",";
                file << "\n";
            }
            file << "  ]\n}\n";
            file.close();
            if (!file)
                throw std::runtime_error("Erro ao gravar drawdown.vtk.series");
        };

        // Reinicia a lista de resultados desta execução.
        writeSeries();

        auto record = [&](const std::string &name, REAL t) {
            std::string basename;
            if (t < 0.) {
                // Estado anterior ao evento em t = 0.
                // -1 é um marcador, não um tempo físico.
                basename = "drawdown_initial";
            } else {
                std::ostringstream filename;
                filename << "drawdown_"
                         << std::setw(6) << std::setfill('0') << outputStep;
                basename = filename.str();
            }

            PostPlasticityUP(mphys, basename + ".vtk");

            // Nome efetivamente gerado pelo escritor do NeoPZ.
            const std::string vtkFile = basename + ".scal_vec.0.vtk";
            std::ifstream check(vtkFile);
            if (!check)
                throw std::runtime_error("VTK nao encontrado: " + vtkFile);

            if (t >= 0.) {
                vtkFrames.emplace_back(vtkFile, t);
                writeSeries();
                ++outputStep;
                std::cout << name << ": t = " << t << " dias"
                          << " (" << t / 365.25 << " anos)\n";
            } else {
                std::cout << name << ": " << vtkFile << "\n";
            }
        };

        // 1. Equilíbrio inicial: reservatório na crista, lambda = 0.
        ok = an.Run({TS(0., 0., 0.)}, nullptr,
                    TS(-1.e4 * Tc, 0., 0.));
        if (ok) record("Inicial", -1.);

        // 2. Evento instantâneo em t = 0: Dt = 0, lambda = 1.
        // Com a carga -14 configurada anteriormente, aplica também
        // a carga vertical de 20 kN/m antes da consolidação.
        if (ok)
            ok = an.Run({TS(0., 1., 0.)}, nullptr, TS(0., 0., 0.));
        if (ok) record("Apos o rebaixamento e a carga", 0.);

        // 3. Consolidação durante 10 anos.
        if (ok) {
            const REAL tFinal = 20. * 365.25; // 3652,5 dias
            const int nSteps = 120;
            const REAL tFirst = std::min<REAL>(1.e-3 * Tc, tFinal / nSteps);

            std::vector<TS> steps;
            steps.reserve(nSteps);
            for (int j = 1; j <= nSteps; ++j) {
                const REAL r = REAL(j - 1) / REAL(nSteps - 1);
                const REAL t = (j == nSteps)
                    ? tFinal
                    : tFirst * std::pow(tFinal / tFirst, r);
                steps.emplace_back(t, 1., 0.);
            }

            ok = an.Run(
                steps,
                [&](int k, const TS &st) {
                    // k = 0 corresponde ao estado já registrado em t = 0.
                    if (k > 0) record("Consolidacao", st.fTime);
                },
                TS(0., 1., 0.)
            );
        }
    } catch (const std::exception &e) {
        std::cerr << "Erro na analise ou na exportacao: " << e.what() << "\n";
        ok = false;
    }

    TPZManVector<TPZCompMesh *, 2> atomic(mphys->MeshVector());
    delete mphys;
    for (auto *m : atomic) delete m;
    delete gmesh;

    if (!ok) std::cout << "Analise interrompida. Consulte as mensagens anteriores.\n";
    return ok;
}

int main(int argc, char *argv[]) {
    const int ref = 2;
    Soil s;
    s.nu = 0.3;
    Water w;

    if (w.H1 < 30. || w.H1 > w.H0) {
        std::cout << "H1 deve estar entre 30 m e H0.\n";
        return 1;
    }

    const TMCVoigt model = ModelVoigt(s);
    const REAL M = s.E * (1. - s.nu) /
                   ((1. + s.nu) * (1. - 2. * s.nu));
    const REAL Tc = 100. / ((w.kh / w.gammaw) * M);

    return Consolidation10Years(model, s, w, Tc, ref) ? 0 : 1;
}
