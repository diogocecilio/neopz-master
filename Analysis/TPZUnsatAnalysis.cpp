/**
 * @file TPZUnsatAnalysis.cpp
 * @brief Implementation of TPZUnsatAnalysis.
 */
#include "TPZUnsatAnalysis.h"

#include "TPZEquationFilter.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <cmath>
#include <iostream>

TPZUnsatAnalysis::TPZUnsatAnalysis(TPZMultiphysicsCompMesh *mesh, const std::map<int, std::array<bool, 2>> &fixU)
    : TPZLinearAnalysis(mesh, true), fMPhys(mesh) {
    for (auto &m : mesh->MaterialVec())
        if (m.first > 0)
            if (auto *u = dynamic_cast<TMat *>(m.second)) fMats.push_back(u);
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

void TPZUnsatAnalysis::SetStrengthReductionFactor(REAL F) {
    for (auto *m : fMats) m->PlasticModel().SetStrengthReductionFactor(F);
}

void TPZUnsatAnalysis::SetGravityFactor(REAL lam) {
    for (auto *m : fMats) m->SetGravityFactor(lam);
}

void TPZUnsatAnalysis::SetGravityScalesSuction(bool b) {
    for (auto *m : fMats) m->SetGravityScalesSuction(b);
}

std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> TPZUnsatAnalysis::SaveMemory() {
    std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> s;
    for (auto *m : fMats) s.push_back(*m->GetMemory());
    return s;
}

void TPZUnsatAnalysis::RestoreMemory(const std::vector<TPZAdmChunkVector<TPZUnsaturatedMem>> &s) {
    for (size_t i = 0; i < fMats.size(); ++i) {
        auto &mem = *fMats[i]->GetMemory();
        for (int64_t k = 0; k < mem.NElements(); ++k) mem[k] = s[i][k];
    }
}

void TPZUnsatAnalysis::ResetMemory() {
    for (auto *m : fMats) {
        auto mem = m->GetMemory();
        const TPZUnsaturatedMem def = m->GetDefaultMemory();
        for (int64_t i = 0; i < mem->NElements(); ++i) (*mem)[i] = def;
    }
}

std::set<int64_t> TPZUnsatAnalysis::ConstrainedEqs(int id, int comp) const {
    auto it = fFixById.find(id);
    return it == fFixById.end() ? std::set<int64_t>() : it->second[comp];
}

std::map<int64_t, REAL> TPZUnsatAnalysis::AllDisplacements(const TPZFMatrix<STATE> &x) const {
    std::map<int64_t, REAL> f;
    for (int64_t i = 0; i < fNEq; ++i)
        if (!fIsP[i]) f[i] = x.GetVal(i, 0);
    return f;
}

std::map<int64_t, REAL> TPZUnsatAnalysis::FrozenFixed(const TPZFMatrix<STATE> &x) const {
    std::map<int64_t, REAL> f = fFixU;
    for (int64_t i = 0; i < fNEq; ++i)
        if (fIsP[i]) f[i] = x.GetVal(i, 0);
    return f;
}

std::map<int64_t, REAL> TPZUnsatAnalysis::BoundaryNodes(int id, bool horizontal) const {
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

TPZUnsatAnalysis::StepResult TPZUnsatAnalysis::Solve(REAL dt, const std::map<int64_t, REAL> &fixed,
                                                     const std::map<int64_t, REAL> &fq, bool frozen, int linesearch) {
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
        // a singular tangent (e.g. at the collapse of a factor of safety search) makes the LU stop with an
        // exception (DebugStop): the step is then reported as not converged
        try {
            TPZLinearAnalysis::Solve();
        } catch (...) {
            return res;
        }
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

void TPZUnsatAnalysis::Accept() {
    for (auto *m : fMats) m->SetUpdateMem(true);
    AssembleFullResidual();
    for (auto *m : fMats) m->SetUpdateMem(false);
    fX = fMPhys->Solution();
}

TPZFMatrix<STATE> TPZUnsatAnalysis::Residual(const TPZFMatrix<STATE> &x, REAL dt, bool frozen) {
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

TPZFMatrix<STATE> TPZUnsatAnalysis::TangentTimes(const TPZFMatrix<STATE> &x, REAL dt, bool frozen,
                                                 const TPZFMatrix<STATE> &d) {
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

TPZFMatrix<STATE> TPZUnsatAnalysis::FullResidual() {
    AssembleFullResidual();
    TPZFMatrix<STATE> r = Rhs();
    r *= -1.;
    return r;
}

void TPZUnsatAnalysis::LoadSolution(const TPZFMatrix<STATE> &sol) {
    TPZLinearAnalysis::LoadSolution(sol);
    fMPhys->LoadSolutionFromMultiPhysics();
}

void TPZUnsatAnalysis::ResetDisplacements() {
    for (auto *m : fMats) m->ResetDisplacementReference();
    for (int64_t i = 0; i < fNEq; ++i)
        if (!fIsP[i]) fX(i, 0) = 0.;
    LoadSolution(fX);
}

void TPZUnsatAnalysis::AssembleSystem() {
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

void TPZUnsatAnalysis::AssembleFullResidual() {
    const TPZEquationFilter saved(fStructMatrix->EquationFilter());
    fStructMatrix->EquationFilter().Reset();
    AssembleResidual();
    fStructMatrix->EquationFilter() = saved;
}

void TPZUnsatAnalysis::SetFilter(const std::map<int64_t, REAL> &fixed) {
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

int64_t TPZUnsatAnalysis::Position(int64_t connect) const {
    const int64_t seq = fMPhys->ConnectVec()[connect].SequenceNumber();
    return fMPhys->Block().Position(seq);
}

void TPZUnsatAnalysis::IdentifyEquations(const std::map<int, std::array<bool, 2>> &fixU) {
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
