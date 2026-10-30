// Incremental analysis of the u - p_w - p_a model: backward Euler step solved by Newton with the full Jacobian,
// Dirichlet equations eliminated by an equation filter, nodal water inflows, memory of the converged state.
// Port of Projects/SlopeRainfall/UnsatSlope.h (branch claude/relaxed-volta-1r4gog), with the constraints and the
// boundary nodes taken from the boundary ids of the mesh.
#ifndef UNSATANALYSIS_H
#define UNSATANALYSIS_H

#include "UnsatModel.h"
#include "TPZLinearAnalysis.h"
#include "TPZSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZEquationFilter.h"
#include "pzstepsolver.h"

#include <chrono>
#include <iostream>

class UnsatAnalysis : public TPZLinearAnalysis {
public:
    struct StepResult {
        bool conv = false;
        int iter = 0;
        REAL eu = 0., normp = 0.; ///< last relative momentum residual and norm of the pressure residual
    };

    /// fixU: boundary id -> (fix u_x, fix u_y) for all the displacement connects of the boundary
    UnsatAnalysis(TPZMultiphysicsCompMesh *mesh, const std::map<int, std::array<bool, 2>> &fixU)
        : TPZLinearAnalysis(mesh, true), fMPhys(mesh) {
        for (auto &m : mesh->MaterialVec())
            if (m.first > 0)
                if (auto *u = dynamic_cast<TUnsat *>(m.second)) fMats.push_back(u);
        if (fMats.empty()) DebugStop();
#ifdef PZ_USING_MKL
        TPZSpStructMatrix<STATE> smat(mesh);
#else
        TPZSkylineNSymStructMatrix<STATE> smat(mesh);
#endif
        SetStructuralMatrix(smat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        fNEq = mesh->NEquations();
        IdentifyEquations(fixU);
        fX = mesh->Solution();
    }

    const std::vector<TUnsat *> &Materials() const { return fMats; }
    TPZMultiphysicsCompMesh *Mesh() { return fMPhys; }
    TPZGeoMesh *GMesh() { return fMPhys->Reference(); }

    void SetStrengthReductionFactor(REAL F) {
        for (auto *m : fMats) m->PlasticModel().SetStrengthReductionFactor(F);
    }
    void SetGravityFactor(REAL lam) {
        for (auto *m : fMats) m->SetGravityFactor(lam);
    }
    void SetGravityScalesSuction(bool b) {
        for (auto *m : fMats) m->SetGravityScalesSuction(b);
    }

    std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> SaveMemory() {
        std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> s;
        for (auto *m : fMats) s.push_back(*m->GetMemory());
        return s;
    }
    void RestoreMemory(const std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> &s) {
        for (size_t i = 0; i < fMats.size(); ++i) {
            auto &mem = *fMats[i]->GetMemory();
            for (int64_t k = 0; k < mem.NElements(); ++k) mem[k] = s[i][k];
        }
    }
    void ResetMemory() {
        for (auto *m : fMats) {
            auto mem = m->GetMemory();
            const TPZUnsaturatedMem def = m->GetDefaultMemory();
            for (int64_t i = 0; i < mem->NElements(); ++i) (*mem)[i] = def;
        }
    }

    /// node -> equation of p_w / p_a, node -> first displacement equation (vertices)
    const std::map<int64_t, int64_t> &WaterEq() const { return fEqW; }
    const std::map<int64_t, int64_t> &AirEq() const { return fEqA; }
    const std::map<int64_t, int64_t> &DispEq() const { return fEqU; }
    int64_t NEq() const { return fNEq; }
    bool IsPressure(int64_t eq) const { return fIsP[eq]; }
    /// displacement constraints of the boundaries (equation -> 0)
    const std::map<int64_t, REAL> &DisplacementConstraints() const { return fFixU; }
    /// constrained displacement equations of a boundary id, component comp (0: x, 1: y)
    std::set<int64_t> ConstrainedEqs(int id, int comp) const {
        auto it = fFixById.find(id);
        return it == fFixById.end() ? std::set<int64_t>() : it->second[comp];
    }
    /// all the displacement equations at their values in x (rigid skeleton)
    std::map<int64_t, REAL> AllDisplacements(const TPZFMatrix<STATE> &x) const {
        std::map<int64_t, REAL> f;
        for (int64_t i = 0; i < fNEq; ++i)
            if (!fIsP[i]) f[i] = x.GetVal(i, 0);
        return f;
    }

    /// vertices of the boundary elements with the id, with their tributary length (half of each edge; horizontal
    /// projection if horizontal = true)
    std::map<int64_t, REAL> BoundaryNodes(int id, bool horizontal = false) const {
        std::map<int64_t, REAL> nodes;
        TPZGeoMesh *g = fMPhys->Reference();
        for (int64_t i = 0; i < g->NElements(); ++i) {
            TPZGeoEl *gel = g->Element(i);
            if (!gel || gel->Dimension() != 1 || gel->MaterialId() != id) continue;
            const TPZGeoNode &a = *gel->NodePtr(0), &b = *gel->NodePtr(1);
            const REAL L = horizontal ? std::fabs(b.Coord(0) - a.Coord(0))
                                      : std::hypot(b.Coord(0) - a.Coord(0), b.Coord(1) - a.Coord(1));
            nodes[gel->NodeIndex(0)] += L / 2.;
            nodes[gel->NodeIndex(1)] += L / 2.;
        }
        return nodes;
    }
    REAL Coord(int64_t node, int d) const { return fMPhys->Reference()->NodeVec()[node].Coord(d); }

    void SetTolerances(REAL tolu, REAL tolp, int maxit, int maxitfrozen) {
        fTolU = tolu;
        fTolP = tolp;
        fMaxIt = maxit;
        fMaxItFrozen = maxitfrozen;
    }
    void SetVerbose(int v) { fVerbose = v; }

    /// converged solution (total u, p_w, p_a)
    TPZFMatrix<STATE> &X() { return fX; }

    /**
     * One backward-Euler step from the converged state: Newton with the equations of fixed eliminated and the
     * water inflow dt * fq (m3/s per unit thickness) added to the water equations. On convergence the solution is
     * loaded in the mesh; the memory and fX are updated only by Accept.
     */
    StepResult Solve(REAL dt, const std::map<int64_t, REAL> &fixed, const std::map<int64_t, REAL> &fq, bool frozen,
                     int linesearch = 4) {
        StepResult res;
        for (auto *m : fMats) {
            m->SetTimeStep(dt);
            m->SetFrozenPressures(frozen);
        }
        SetFilter(fixed);
        TPZFMatrix<STATE> x(fX);
        for (auto &f : fixed) x(f.first, 0) = f.second;
        LoadSolution(x);
        for (auto *m : fMats) m->SetAssembleGravityOnly(true);
        AssembleFullResidual();
        for (auto *m : fMats) m->SetAssembleGravityOnly(false);
        REAL fref = 0.;
        {
            TPZFMatrix<STATE> &g = Rhs();
            for (int64_t i = 0; i < fNEq; ++i)
                if (!fIsP[i]) fref += std::pow(g.GetVal(i, 0), 2);
        }
        fref = std::max(std::sqrt(fref), fForceRef);
        REAL eu = 0., normp = 0.;
        auto evaluate = [&](const TPZFMatrix<STATE> &xt) {
            LoadSolution(xt);
            AssembleSystem();
            TPZFMatrix<STATE> &rhs = Rhs();
            for (auto &q : fq) rhs(q.first, 0) += dt * q.second;
            REAL nu2 = 0., np2 = 0.;
            for (int64_t i = 0; i < fNEq; ++i) {
                if (fixed.count(i)) continue;
                const REAL r = rhs.GetVal(i, 0);
                (fIsP[i] ? np2 : nu2) += r * r;
            }
            eu = std::sqrt(nu2) / fref;
            normp = std::sqrt(np2);
        };
        evaluate(x);
        REAL lastdp = frozen ? 0. : 1.e300;
        const int maxit = frozen ? fMaxItFrozen : fMaxIt;
        for (int it = 1; it <= maxit; ++it) {
            res.iter = it;
            res.eu = eu;
            res.normp = normp;
            if (fVerbose) std::cout << "      it " << it << " eu = " << eu << " |Rp| = " << normp << " |dp| = " << lastdp << "\n";
            if (!std::isfinite(eu) || !std::isfinite(normp) || eu > 1.e6) return res;
            if (eu < fTolU && lastdp < fTolP && it > 1) {
                res.conv = true;
                return res;
            }
            TPZLinearAnalysis::Solve();
            TPZFMatrix<STATE> dx = Solution();
            lastdp = 0.;
            for (int64_t i = 0; i < fNEq; ++i)
                if (fIsP[i]) lastdp = std::max(lastdp, std::fabs(dx.GetVal(i, 0)));
            const REAL mref = normp;
            REAL lam = 1.;
            const bool ls = linesearch > 0 && dt > 0. && !frozen && lastdp >= fTolP;
            for (int k = 0;; ++k) {
                TPZFMatrix<STATE> xt(x);
                for (int64_t i = 0; i < fNEq; ++i) xt(i, 0) += lam * dx.GetVal(i, 0);
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
    void Restore() { LoadSolution(fX); }

    /// Residual R (all equations, no nodal inflow) of x with the step dt (the vector of NeoPZ is -R)
    TPZFMatrix<STATE> Residual(const TPZFMatrix<STATE> &x, REAL dt, bool frozen) {
        for (auto *m : fMats) {
            m->SetTimeStep(dt);
            m->SetFrozenPressures(frozen);
        }
        LoadSolution(x);
        AssembleFullResidual();
        TPZFMatrix<STATE> r = Rhs();
        r *= -1.;
        return r;
    }
    /// Product of the tangent dR/dX of x with d (all equations; no filter)
    TPZFMatrix<STATE> TangentTimes(const TPZFMatrix<STATE> &x, REAL dt, bool frozen, const TPZFMatrix<STATE> &d) {
        for (auto *m : fMats) {
            m->SetTimeStep(dt);
            m->SetFrozenPressures(frozen);
        }
        SetFilter({});
        LoadSolution(x);
        AssembleSystem();
        auto *solver = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
        TPZFMatrix<STATE> jd;
        solver->Matrix()->Multiply(d, jd);
        return jd;
    }

    /// Water (or air) residual of the solution loaded in the mesh: R_w / dt at a vertex with p_w prescribed is the
    /// inflow there
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

    /// u = 0 keeping the effective stress
    void ResetDisplacements() {
        for (auto *m : fMats) m->ResetDisplacementReference();
        for (int64_t i = 0; i < fNEq; ++i)
            if (!fIsP[i]) fX(i, 0) = 0.;
        LoadSolution(fX);
    }

    REAL Nodal(const TPZFMatrix<STATE> &x, int64_t node, int space) const {
        return x.GetVal((space == 1 ? fEqW : fEqA).at(node), 0);
    }
    REAL NodalDisplacement(const TPZFMatrix<STATE> &x, int64_t node, int comp) const {
        return x.GetVal(fEqU.at(node) + comp, 0);
    }
    /// reference of the momentum residual when the weight is null (kN)
    void SetForceReference(REAL f) { fForceRef = f; }

private:
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
        auto *solver = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
        if (solver) solver->ResetMatrix();
        fTemplate = nullptr;
    }
    int64_t Position(int64_t connect) const {
        const int64_t seq = fMPhys->ConnectVec()[connect].SequenceNumber();
        return fMPhys->Block().Position(seq);
    }
    void IdentifyEquations(const std::map<int, std::array<bool, 2>> &fixU) {
        TPZVec<TPZCompMesh *> &mv = fMPhys->MeshVector();
        const int64_t off[3] = {0, mv[0]->NConnects(), mv[0]->NConnects() + mv[1]->NConnects()};
        fIsP.assign(fNEq, false);
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
                if (s == 0 && id > 0)
                    for (int i = 0; i < 3; ++i) fEqU[gel->NodeIndex(i)] = Position(off[0] + cel->ConnectIndex(i));
                auto it = fixU.find(id);
                if (s == 0 && id < 0 && it != fixU.end())
                    for (int i = 0; i < cel->NConnects(); ++i) {
                        const int64_t ci = off[0] + cel->ConnectIndex(i);
                        const int64_t pos = Position(ci);
                        const int size = fMPhys->Block().Size(fMPhys->ConnectVec()[ci].SequenceNumber());
                        for (int k = 0; k < size; ++k)
                            if (it->second[k % 2]) fFixU[pos + k] = 0., fFixById[id][k % 2].insert(pos + k);
                    }
            }
        }
    }

    TPZMultiphysicsCompMesh *fMPhys;
    std::vector<TUnsat *> fMats;
    std::map<int64_t, int64_t> fEqW, fEqA, fEqU;
    std::map<int64_t, REAL> fFixU;
    std::map<int, std::array<std::set<int64_t>, 2>> fFixById;
    std::vector<bool> fIsP;
    std::set<int64_t> fFilterKey;
    TPZAutoPointer<TPZMatrix<STATE>> fTemplate;
    TPZFMatrix<STATE> fX;
    REAL fTolU = 1.e-7, fTolP = 1.e-6, fForceRef = 1.e-12;
    int fMaxIt = 40, fMaxItFrozen = 25, fVerbose = 0;
    int64_t fNEq = 0;
};

#endif
