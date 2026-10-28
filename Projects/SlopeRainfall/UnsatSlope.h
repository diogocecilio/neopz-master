// Homogeneous slope of Rahardjo, Ong, Rezaur & Leong (2007), J. Geotech. Geoenviron. Eng. 133(12), under rainfall,
// with the three-phase material TPZMatUnsaturatedUPwPa (skeleton u, water p_w, air p_a; Mohr-Coulomb in the Bishop
// stress): geometric mesh, multiphysics mesh and the incremental analysis (Newton with the Dirichlet equations
// eliminated, unilateral rainfall condition, factor of safety with the pressures of an instant frozen).
// Port of nsat3.py / caso.py (talude3f, Python reference code).
#ifndef UNSATSLOPE_H
#define UNSATSLOPE_H

#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "TPZGeoLinear.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZStructMatrix.h"
#include "TPZEquationFilter.h"
#include "pzstepsolver.h"
#include "pzgeotriangle.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"
#include "pzmultiphysicselement.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef SLOPERAINFALL_PV
// legacy iterative Mohr-Coulomb (TPZPlasticStepPV): accepts a dilatancy angle psi != phi (executable SlopeRainfallPV)
using TMC = TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;
#else
// closed-form RHW Mohr-Coulomb with consistent tangent: associative (psi = phi)
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
#endif
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;

/// Boundary ids of the mesh (gmsh physical groups 1..6, scripts/malha_talude.py)
enum EBoundary { EBase = -1, ELeft = -2, ERight = -3, ECrest = -4, EFace = -5, EToe = -6 };

/// Reads a gmsh 2.2 ASCII mesh of linear triangles (type 2, material = physical tag > 0: soil layer) and lines
/// (type 1, material -physical tag)
inline TPZGeoMesh *ReadMsh(const std::string &file) {
    std::ifstream in(file);
    if (!in) {
        std::cout << "mesh file " << file << " not found\n";
        DebugStop();
    }
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    std::string line;
    std::map<int64_t, int64_t> idx;
    while (std::getline(in, line)) {
        if (line.rfind("$Nodes", 0) == 0) {
            int64_t n;
            in >> n;
            gmesh->NodeVec().Resize(n);
            for (int64_t i = 0; i < n; ++i) {
                int64_t id;
                REAL x, y, z;
                in >> id >> x >> y >> z;
                idx[id] = i;
                TPZManVector<REAL, 3> co = {x, y, 0.};
                gmesh->NodeVec()[i] = TPZGeoNode(i, co, *gmesh);
            }
        } else if (line.rfind("$Elements", 0) == 0) {
            int64_t m;
            in >> m;
            std::getline(in, line);
            for (int64_t k = 0; k < m; ++k) {
                std::getline(in, line);
                std::istringstream ss(line);
                int64_t id, type, ntag;
                ss >> id >> type >> ntag;
                std::vector<int64_t> tags(ntag);
                for (auto &t : tags) ss >> t;
                std::vector<int64_t> nodes;
                int64_t v;
                while (ss >> v) nodes.push_back(idx.at(v));
                if (type == 2) {
                    // counterclockwise
                    const auto &a = gmesh->NodeVec()[nodes[0]], &b = gmesh->NodeVec()[nodes[1]], &c = gmesh->NodeVec()[nodes[2]];
                    const REAL det = (b.Coord(0) - a.Coord(0)) * (c.Coord(1) - a.Coord(1)) -
                                     (c.Coord(0) - a.Coord(0)) * (b.Coord(1) - a.Coord(1));
                    if (det < 0.) std::swap(nodes[1], nodes[2]);
                    TPZManVector<int64_t, 3> nn = {nodes[0], nodes[1], nodes[2]};
                    new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nn, ntag > 0 && tags[0] > 0 ? int(tags[0]) : 1,
                                                                    *gmesh);
                } else if (type == 1) {
                    TPZManVector<int64_t, 2> nn = {nodes[0], nodes[1]};
                    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nn, -int(tags[0]), *gmesh);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

/// Soil and fluid parameters of the analyses (Table 1 of Rahardjo et al. 2007 and the hypotheses of talude3f)
struct UnsatParams {
    REAL E = 2.e4, nu = 0.3, phi = 26. * M_PI / 180., psi = 26. * M_PI / 180., c = 10., gamma = 20.;
    REAL chi = 1.; ///< phi^b = phi'
    TUnsat::EBishop bishop = TUnsat::EChiConstant; ///< chi constant (on the suction), chi = S_w or chi = S_e
    REAL n = 0.45, Kw = 2.2e6, gw = 9.81, ga = 1.2 * 9.81e-3, muRatio = 1.e-3 / 1.8e-5, patm = 101.325;
    REAL ks = 1.e-5; ///< saturated conductivity (m/s)
    /// anisotropic conductivity (m/s): principal values k_x (direction at theta from the x axis) and k_y; <= 0: k_s
    REAL kx = 0., ky = 0., theta = 0.; ///< theta in degrees, counterclockwise
    TPZUnsaturatedHydraulics hydro;
    UnsatParams() {
        hydro.fRetention = TPZUnsaturatedHydraulics::EFredlundXing;
        hydro.fA = 50.;
        hydro.fN = hydro.fM = 1.;
        hydro.fKrw = TPZUnsaturatedHydraulics::EPowerLaw;
        hydro.fPower = 4.;
        hydro.fKrwMin = 1.e-12;
        hydro.fKra = TPZUnsaturatedHydraulics::EBrooksCorey;
        hydro.fLambdaBC = 3.;
        hydro.fKraMin = 1.e-4;
    }
};

inline TMC MohrCoulomb(const UnsatParams &p) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(p.E, p.nu);
#ifdef SLOPERAINFALL_PV
    TMC pv;
    pv.fYC.SetUp(p.phi, p.psi, p.c, ER);
    pv.fER = ER;
    return pv;
#else
    if (std::fabs(p.psi - p.phi) > 1.e-12) {
        std::cout << "the closed-form Mohr-Coulomb is associative: use SlopeRainfallPV for psi != phi\n";
        DebugStop();
    }
    return TMC(TPZYCMohrCoulombPV2(p.phi, p.psi, p.c, ER), ER);
#endif
}

/// Material ids (> 0) of the triangles of the mesh: the soil layers
inline std::set<int> DomainIds(TPZGeoMesh *gmesh) {
    std::set<int> ids;
    for (int64_t i = 0; i < gmesh->NElements(); ++i) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && gel->Dimension() == 2) ids.insert(gel->MaterialId());
    }
    return ids;
}

/// H1 atomic space with a null material on the domain (each layer) and on all the boundaries
inline TPZCompMesh *AtomicMesh(TPZGeoMesh *gmesh, int nstate, int order) {
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    TPZNullMaterial<STATE> *mat = nullptr;
    for (int id : DomainIds(gmesh)) {
        mat = new TPZNullMaterial<STATE>(id, 2, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id = -6; id <= -1; ++id) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

/// Unsaturated material of a soil layer
inline TUnsat *UnsatMaterial(int id, const UnsatParams &p) {
    auto *mat = new TUnsat(id);
    mat->SetPlasticModel(MohrCoulomb(p));
    mat->SetHydraulics(p.hydro);
    mat->SetPorosity(p.n);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(p.Kw);
    mat->SetConductivity(p.ks, p.gw, p.muRatio);
    mat->SetAnisotropicConductivity(p.kx > 0. ? p.kx : p.ks, p.ky > 0. ? p.ky : p.ks, p.theta * M_PI / 180.);
    mat->SetAir(p.ga, p.patm);
    mat->SetUnitWeight(p.gamma);
    mat->SetBishop(p.bishop, p.chi);
    mat->SetLumpedStorage(true);
    mat->SetIntegrationOrder(2);
    return mat;
}

/// Multiphysics mesh u (P2) - p_w (P1) - p_a (P1) with an unsaturated material per layer (parameters of the layer
/// id in layers, otherwise p; 3-point rule, 3 memory items)
inline TPZMultiphysicsCompMesh *UnsatMesh(TPZGeoMesh *gmesh, const UnsatParams &p,
                                          const std::map<int, UnsatParams> &layers = {}) {
    TPZManVector<TPZCompMesh *, 3> meshvec = {AtomicMesh(gmesh, 2, 2), AtomicMesh(gmesh, 1, 1), AtomicMesh(gmesh, 1, 1)};
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    const std::set<int> ids = DomainIds(gmesh);
    TUnsat *mat = nullptr;
    for (int id : ids) {
        auto it = layers.find(id);
        mat = UnsatMaterial(id, it == layers.end() ? p : it->second);
        mphys->InsertMaterialObject(mat);
    }
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    std::set<int> bcids;
    for (int id = -6; id <= -1; ++id) {
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, val1, val2)); // no contribution
        bcids.insert(id);
    }
    TPZManVector<int, 3> active(3, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, ids, bcids);
    return mphys;
}

/**
 * Incremental analysis of the unsaturated slope. The Dirichlet equations (fixed displacements, prescribed p_w and p_a,
 * ponded vertices) are eliminated with the equation filter; the rainfall enters as nodal inflows of the water
 * equations of the surface vertices (unilateral condition: p_w <= 0, inflow <= rain, complementarity), as in
 * unilateral_step of nsat3.py.
 */
class UnsatAnalysis : public TPZLinearAnalysis {
public:
    struct StepResult {
        bool conv = false;
        int iter = 0;
    };

    UnsatAnalysis(TPZMultiphysicsCompMesh *mesh, int nthreads = 0)
        : TPZLinearAnalysis(mesh, true), fMPhys(mesh) {
        for (auto &m : mesh->MaterialVec())
            if (m.first > 0)
                if (auto *u = dynamic_cast<TUnsat *>(m.second)) fMats.push_back(u);
        if (fMats.empty()) DebugStop();
        fMat = fMats[0];
#ifdef PZ_USING_MKL
        TPZSpStructMatrix<STATE> smat(mesh);
#else
        TPZSkylineNSymStructMatrix<STATE> smat(mesh);
#endif
        smat.SetNumThreads(nthreads);
        SetStructuralMatrix(smat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        fNEq = mesh->NEquations();
        IdentifyEquations();
        fX = mesh->Solution();
    }

    /// material of the first layer, and all the layers
    TUnsat *Material() { return fMat; }
    const std::vector<TUnsat *> &Materials() const { return fMats; }
    TPZMultiphysicsCompMesh *Mesh() { return fMPhys; }

    /// settings applied to the materials of all the layers
    void SetStrengthReductionFactor(REAL F) {
        for (auto *m : fMats) m->PlasticModel().SetStrengthReductionFactor(F);
    }
    void SetGravityFactor(REAL lam) {
        for (auto *m : fMats) m->SetGravityFactor(lam);
    }
    void SetGravityScalesSuction(bool b) {
        for (auto *m : fMats) m->SetGravityScalesSuction(b);
    }
    void SetHydroMechanicalCoupling(bool b) {
        for (auto *m : fMats) m->SetHydroMechanicalCoupling(b);
    }
    /// copies of the memory of the integration points of all the layers
    std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> SaveMemory() {
        std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> s;
        for (auto *m : fMats) s.push_back(*m->GetMemory());
        return s;
    }
    void RestoreMemory(const std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> &s) {
        for (size_t i = 0; i < fMats.size(); ++i) *fMats[i]->GetMemory() = s[i];
    }
    /// memory back to the default (stress-free) state
    void ResetMemory() {
        for (auto *m : fMats) {
            auto mem = m->GetMemory();
            const TPZUnsaturatedMem def = m->GetDefaultMemory();
            for (int64_t i = 0; i < mem->NElements(); ++i) (*mem)[i] = def;
        }
    }

    /// node -> equation of p_w / p_a (vertices) and of the displacement components
    const std::map<int64_t, int64_t> &WaterEq() const { return fEqW; }
    const std::map<int64_t, int64_t> &AirEq() const { return fEqA; }
    int64_t NEq() const { return fNEq; }
    bool IsPressure(int64_t eq) const { return fIsP[eq]; }

    /// converged solution (total u, p_w, p_a)
    TPZFMatrix<STATE> &X() { return fX; }

    /// Constrained displacement equations: base fixed, lateral sides on rollers
    std::map<int64_t, REAL> DisplacementConstraints() const { return fFixU; }

    /// Newton tolerance on the momentum residual (relative to the weight) and on the pressure correction (kPa) and
    /// maximum iterations of the coupled steps and of the frozen-pressure (factor of safety) steps
    void SetTolerances(REAL tolu, REAL tolp, int maxit, int maxitfrozen) {
        fTolU = tolu;
        fTolP = tolp;
        fMaxIt = maxit;
        fMaxItFrozen = maxitfrozen;
    }
    void SetVerbose(int v) { fVerbose = v; }
    int64_t NIterations() const { return fNIter; }
    /// time (s) spent in the assembly of the tangent systems and in the linear solver
    double TimeAssemble() const { return fTAssemble; }
    double TimeSolve() const { return fTSolve; }

    /**
     * One backward-Euler step from the converged state (memory of the material and fX): Newton with the equations of
     * fixed eliminated, water inflow dt * fq added to the water equations. On convergence the solution is loaded in the
     * mesh but neither the memory nor fX are updated (see Accept). frozen: only the momentum balance (all the
     * pressure equations must be in fixed).
     */
    StepResult Solve(REAL dt, const std::map<int64_t, REAL> &fixed, const std::map<int64_t, REAL> &fq, bool frozen,
                     int linesearch = 4) {
        StepResult res;
        for (auto *m : fMats) {
            m->SetTimeStep(dt);
            m->SetFrozenPressures(frozen);
        }
        SetFilter(fixed);
        const int64_t neq = NEq();
        TPZFMatrix<STATE> x(fX);
        for (auto &f : fixed) x(f.first, 0) = f.second;
        // reference of the momentum residual: the weight of the soil
        LoadSolution(x);
        for (auto *m : fMats) m->SetAssembleGravityOnly(true);
        AssembleFullResidual();
        for (auto *m : fMats) m->SetAssembleGravityOnly(false);
        REAL fgrav = 0.;
        {
            TPZFMatrix<STATE> &g = Rhs();
            for (int64_t i = 0; i < neq; ++i)
                if (!fIsP[i]) fgrav += g.GetVal(i, 0) * g.GetVal(i, 0);
        }
        fgrav = std::max(std::sqrt(fgrav), REAL(1.e-12));

        REAL eu, normp;
        auto evaluate = [&](const TPZFMatrix<STATE> &xt) {
            LoadSolution(xt);
            auto ta = std::chrono::steady_clock::now();
            AssembleSystem();
            fTAssemble += std::chrono::duration<double>(std::chrono::steady_clock::now() - ta).count();
            fNIter++;
            TPZFMatrix<STATE> &rhs = Rhs();
            for (auto &q : fq) rhs(q.first, 0) += dt * q.second;
            REAL nu2 = 0., np2 = 0.;
            for (int64_t i = 0; i < neq; ++i) {
                if (fixed.count(i)) continue;
                const REAL r = rhs.GetVal(i, 0);
                (fIsP[i] ? np2 : nu2) += r * r;
            }
            eu = std::sqrt(nu2) / fgrav;
            normp = std::sqrt(np2);
        };
        evaluate(x);
        REAL lastdp = frozen ? 0. : 1.e300;
        const int maxit = frozen ? fMaxItFrozen : fMaxIt;
        for (int it = 1; it <= maxit; ++it) {
            res.iter = it;
            if (fVerbose > 1)
                std::cout << "      it " << it << " eu = " << eu << " |Rp| = " << normp << " |dp| = " << lastdp << "\n";
            if (!std::isfinite(eu) || !std::isfinite(normp) || eu > 1.e6) return res;
            if (eu < fTolU && lastdp < fTolP && it > 1) {
                res.conv = true;
                return res;
            }
            auto ts = std::chrono::steady_clock::now();
            TPZLinearAnalysis::Solve();
            fTSolve += std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
            TPZFMatrix<STATE> dx = Solution();
            lastdp = 0.;
            for (int64_t i = 0; i < neq; ++i)
                if (fIsP[i]) lastdp = std::max(lastdp, std::fabs(dx.GetVal(i, 0)));
            const REAL mref = normp;
            REAL lam = 1.;
            // no line search once the pressure correction is below the tolerance (at a steady state the residual
            // is at the round-off level and cannot decrease)
            const bool ls = linesearch > 0 && dt > 0. && !frozen && lastdp >= fTolP;
            for (int k = 0;; ++k) {
                TPZFMatrix<STATE> xt(x);
                for (int64_t i = 0; i < neq; ++i) xt(i, 0) += lam * dx.GetVal(i, 0);
                evaluate(xt);
                if (!ls || k == linesearch || normp < mref) {
                    x = xt;
                    break;
                }
                lam *= 0.5;
            }
            if (lam < 1.) lastdp = 1.e300;
        }
        return res;
    }

    /// Updates the memory with the solution loaded in the mesh and stores it as the converged state
    void Accept() {
        for (auto *m : fMats) m->SetUpdateMem(true);
        AssembleFullResidual();
        for (auto *m : fMats) m->SetUpdateMem(false);
        fX = fMPhys->Solution();
    }

    /// Restores the converged state in the mesh
    void Restore() { LoadSolution(fX); }

    /// Water residual (all equations, nodal inflows not included) of the solution loaded in the mesh: R_w / dt is the
    /// inflow needed at the vertices with prescribed p_w
    TPZFMatrix<STATE> FullResidual() {
        AssembleFullResidual();
        TPZFMatrix<STATE> r = Rhs();
        r *= -1.;
        return r;
    }

    void LoadSolution(const TPZFMatrix<STATE> &sol) override {
        TPZLinearAnalysis::LoadSolution(sol);
        fMPhys->LoadSolutionFromMultiPhysics();
    }
    using TPZLinearAnalysis::LoadSolution;

    /// u = 0 keeping the effective stress (eps_p <- eps_p - eps_t)
    void ResetDisplacements() {
        for (auto *m : fMats) m->ResetDisplacementReference();
        for (int64_t i = 0; i < NEq(); ++i)
            if (!fIsP[i]) fX(i, 0) = 0.;
        LoadSolution(fX);
    }

    /// Value of the p_w (space 1) or p_a (space 2) of a node in a solution vector
    REAL Nodal(const TPZFMatrix<STATE> &x, int64_t node, int space) const {
        return x.GetVal((space == 1 ? fEqW : fEqA).at(node), 0);
    }

private:
    /// Assembly of the tangent and the residual; while the equation filter is unchanged the sparse structure is
    /// copied from a template (a clone of the first matrix, before its decomposition) instead of rebuilt
    void AssembleSystem() {
        auto *solver = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
        if (!fTemplate) {
            Assemble();
            fTemplate = dynamic_cast<TPZMatrix<STATE> *>(solver->Matrix()->Clone());
            return;
        }
        TPZMatrix<STATE> *m = dynamic_cast<TPZMatrix<STATE> *>(fTemplate->Clone());
        m->Zero();
        fRhs.Redim(NEq(), 1);
        fStructMatrix->Assemble(*m, fRhs, fGuiInterface);
        solver->SetMatrix(m);
    }

    void AssembleFullResidual() {
        const TPZEquationFilter saved(fStructMatrix->EquationFilter());
        fStructMatrix->EquationFilter().Reset();
        AssembleResidual();
        fStructMatrix->EquationFilter() = saved;
    }

    void SetFilter(const std::map<int64_t, REAL> &fixed) {
        std::set<int64_t> key;
        for (auto &f : fixed) key.insert(f.first);
        if (key == fFilterKey && fStructMatrix->EquationFilter().IsActive()) return;
        fFilterKey = key;
        TPZStack<int64_t> act;
        for (int64_t i = 0; i < NEq(); ++i)
            if (!key.count(i)) act.Push(i);
        TPZEquationFilter &filter = fStructMatrix->EquationFilter();
        filter.Reset();
        filter.SetActiveEquations(act);
        // the matrix of the solver has the size of the previous filter: discard it
        auto *solver = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
        if (solver) solver->ResetMatrix();
        fTemplate = nullptr;
    }

    int64_t Position(int64_t connect) const {
        const int64_t seq = fMPhys->ConnectVec()[connect].SequenceNumber();
        return fMPhys->Block().Position(seq);
    }

    void IdentifyEquations() {
        TPZVec<TPZCompMesh *> &mv = fMPhys->MeshVector();
        const int64_t off[3] = {0, mv[0]->NConnects(), mv[0]->NConnects() + mv[1]->NConnects()};
        const int64_t neq = NEq();
        fIsP.assign(neq, false);
        for (int s = 1; s < 3; ++s)
            for (int64_t ic = 0; ic < mv[s]->NConnects(); ++ic) {
                const TPZConnect &c = fMPhys->ConnectVec()[off[s] + ic];
                if (c.SequenceNumber() < 0) continue;
                const int64_t pos = Position(off[s] + ic);
                for (int k = 0; k < fMPhys->Block().Size(c.SequenceNumber()); ++k) fIsP[pos + k] = true;
            }
        for (int s = 0; s < 3; ++s) {
            TPZCompMesh *cm = mv[s];
            for (int64_t iel = 0; iel < cm->NElements(); ++iel) {
                TPZCompEl *cel = cm->Element(iel);
                if (!cel || !cel->Reference()) continue;
                TPZGeoEl *gel = cel->Reference();
                const int id = gel->MaterialId();
                if (s > 0 && id > 0)
                    for (int i = 0; i < 3; ++i)
                        (s == 1 ? fEqW : fEqA)[gel->NodeIndex(i)] = Position(off[s] + cel->ConnectIndex(i));
                if (s == 0 && (id == EBase || id == ELeft || id == ERight)) {
                    for (int i = 0; i < cel->NConnects(); ++i) {
                        const int64_t ci = off[0] + cel->ConnectIndex(i);
                        const int64_t pos = Position(ci);
                        const int size = fMPhys->Block().Size(fMPhys->ConnectVec()[ci].SequenceNumber());
                        for (int k = 0; k < size; ++k)
                            if (id == EBase || k % 2 == 0) fFixU[pos + k] = 0.;
                    }
                }
            }
        }
    }

    TPZMultiphysicsCompMesh *fMPhys;
    TUnsat *fMat;
    std::vector<TUnsat *> fMats;
    std::map<int64_t, int64_t> fEqW, fEqA;
    std::map<int64_t, REAL> fFixU;
    std::vector<bool> fIsP;
    std::set<int64_t> fFilterKey;
    TPZAutoPointer<TPZMatrix<STATE>> fTemplate;
    TPZFMatrix<STATE> fX;
    REAL fTolU = 1.e-7, fTolP = 1.e-6;
    int fMaxIt = 40, fMaxItFrozen = 25, fVerbose = 0;
    int64_t fNIter = 0, fNEq = 0;
    double fTAssemble = 0., fTSolve = 0.;
};

#endif
