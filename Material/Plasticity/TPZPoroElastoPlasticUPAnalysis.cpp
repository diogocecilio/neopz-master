/**
 * @file TPZPoroElastoPlasticUPAnalysis.cpp
 * @brief Implementation of the incremental u-p Newton driver (see TPZPoroElastoPlasticUPAnalysis.h).
 */

#include "TPZPoroElastoPlasticUPAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicselement.h"
#include "TPZBndCondT.h"
#include "TPZStructMatrix.h"
#include "TPZEquationFilter.h"
#include "pzstack.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "pzerror.h"
#include <cmath>
#include <iomanip>

TPZPoroElastoPlasticUPAnalysis::TPZPoroElastoPlasticUPAnalysis(TPZMultiphysicsCompMesh *mesh,
                                                               TPZMatPoroElastoPlasticUPBase *mat,
                                                               std::ostream &out)
    : TPZLinearAnalysis(mesh, false, out), fMPhys(mesh), fMat(mat) {
    if (!mesh || !mat) DebugStop();
    fConverged = mesh->Solution();
}

void TPZPoroElastoPlasticUPAnalysis::SetControlledDisplacement(int bcid, int component) {
    auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(fMPhys->FindMaterial(bcid));
    if (!bc || !TPZMatPoroElastoPlasticUPBase::IsDirichletU(bc->Type())) {
        std::cout << __PRETTY_FUNCTION__ << " boundary condition " << bcid
                  << " does not exist or is not a Dirichlet condition of the displacement" << std::endl;
        DebugStop();
    }
    fControlled.push_back({bcid, component});
}

void TPZPoroElastoPlasticUPAnalysis::LoadSolution(const TPZFMatrix<STATE> &sol) {
    TPZLinearAnalysis::LoadSolution(sol);
    fMPhys->LoadSolutionFromMultiPhysics();
}

void TPZPoroElastoPlasticUPAnalysis::IdentifyEquations() {
    const int64_t neq = fMPhys->NEquations();
    fConstrained.assign(neq, false);
    fIsPressure.assign(neq, false);
    fNodeConnects.clear();
    TPZVec<TPZCompMesh *> &meshvec = fMPhys->MeshVector();
    const int64_t nconU = meshvec[0]->NConnects();
    const int dim = fMat->DimensionU();
    TPZBlock &block = fMPhys->Block();
    // pressure equations: connects of the second space
    for (int64_t ic = nconU; ic < fMPhys->NConnects(); ++ic) {
        TPZConnect &c = fMPhys->ConnectVec()[ic];
        const int64_t seq = c.SequenceNumber();
        if (seq < 0) continue;
        const int64_t pos = block.Position(seq);
        for (int k = 0; k < block.Size(seq); ++k) fIsPressure[pos + k] = true;
    }
    // node -> vertex connects of each space (same connect indices in the multiphysics mesh)
    for (int space = 0; space < 2; ++space) {
        TPZCompMesh *cmesh = meshvec[space];
        const int64_t offset = (space == 0) ? 0 : nconU;
        for (int64_t iel = 0; iel < cmesh->NElements(); ++iel) {
            TPZCompEl *cel = cmesh->Element(iel);
            if (!cel || !cel->Reference()) continue;
            TPZGeoEl *gel = cel->Reference();
            const int ncorner = gel->NCornerNodes();
            if (cel->NConnects() < ncorner) continue;
            for (int i = 0; i < ncorner; ++i) {
                const int64_t node = gel->NodeIndex(i);
                auto it = fNodeConnects.find(node);
                if (it == fNodeConnects.end()) it = fNodeConnects.insert({node, {-1, -1}}).first;
                if (space == 0) it->second.first = offset + cel->ConnectIndex(i);
                else it->second.second = offset + cel->ConnectIndex(i);
            }
        }
    }
    // edges of the displacement mesh (hierarchical edge function and its two vertex functions)
    fEdgeEquations.clear();
    {
        TPZCompMesh *umesh = meshvec[0];
        std::set<int64_t> visited;
        for (int64_t iel = 0; iel < umesh->NElements(); ++iel) {
            auto *intel = dynamic_cast<TPZInterpolatedElement *>(umesh->Element(iel));
            if (!intel || !intel->Reference()) continue;
            TPZGeoEl *gel = intel->Reference();
            for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side) {
                if (gel->SideDimension(side) != 1) continue;
                const int64_t ic = intel->ConnectIndex(intel->MidSideConnectLocId(side));
                if (visited.count(ic)) continue;
                visited.insert(ic);
                const int64_t seq = fMPhys->ConnectVec()[ic].SequenceNumber();
                if (seq < 0 || block.Size(seq) != dim) continue; // only edges with one (quadratic) function
                const int64_t a = gel->SideNodeIndex(side, 0), b = gel->SideNodeIndex(side, 1);
                const int64_t ia = intel->ConnectIndex(gel->SideNodeLocIndex(side, 0));
                const int64_t ib = intel->ConnectIndex(gel->SideNodeLocIndex(side, 1));
                (void)a;
                (void)b;
                const int64_t sa = fMPhys->ConnectVec()[ia].SequenceNumber(), sb = fMPhys->ConnectVec()[ib].SequenceNumber();
                fEdgeEquations.push_back({block.Position(seq), block.Position(sa), block.Position(sb)});
            }
        }
    }
    // constrained equations: connects of the boundary elements with Dirichlet conditions
    for (int64_t iel = 0; iel < fMPhys->NElements(); ++iel) {
        TPZCompEl *cel = fMPhys->Element(iel);
        auto *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel || !cel->Material()) continue;
        auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(cel->Material());
        if (!bc) continue;
        const int type = bc->Type();
        const bool dirU = TPZMatPoroElastoPlasticUPBase::IsDirichletU(type);
        const bool dirP = TPZMatPoroElastoPlasticUPBase::IsDirichletP(type);
        if (!dirU && !dirP) continue;
        const int ncu = mfel->Element(0) ? mfel->Element(0)->NConnects() : 0;
        for (int i = 0; i < cel->NConnects(); ++i) {
            const bool isU = (i < ncu);
            if ((isU && !dirU) || (!isU && !dirP)) continue;
            TPZConnect &c = cel->Connect(i);
            const int64_t seq = c.SequenceNumber();
            const int64_t pos = block.Position(seq);
            const int size = block.Size(seq);
            if (!isU) {
                for (int k = 0; k < size; ++k) fConstrained[pos + k] = true;
                continue;
            }
            for (int d = 0; d < dim; ++d) {
                if (type == TPZMatPoroElastoPlasticUPBase::EDirichletUDirectional && bc->Val1().GetVal(d, d) == 0.)
                    continue;
                for (int k = d; k < size; k += dim) fConstrained[pos + k] = true;
            }
        }
    }
    fEquationsIdentified = true;
}

void TPZPoroElastoPlasticUPAnalysis::ApplyControlledDisplacement(REAL uc) {
    for (auto &ctrl : fControlled) {
        auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(fMPhys->FindMaterial(ctrl.first));
        TPZManVector<STATE, 3> v2(bc->Val2());
        if (v2.size() < 3) {
            const int n0 = v2.size();
            v2.Resize(3);
            for (int i = n0; i < 3; ++i) v2[i] = 0.;
        }
        v2[ctrl.second] = uc;
        bc->SetVal2(v2);
    }
}

void TPZPoroElastoPlasticUPAnalysis::ImposeDirichletValues(TPZFMatrix<STATE> &sol) {
    const int dim = fMat->DimensionU();
    TPZBlock &block = fMPhys->Block();
    // nodal values of the mid-edge nodes (kept for the free components, as in the nodal basis)
    std::vector<std::array<STATE, 3>> mid(fEdgeEquations.size());
    if (fNodalNorm) {
        for (size_t ie = 0; ie < fEdgeEquations.size(); ++ie) {
            const auto &e = fEdgeEquations[ie];
            for (int d = 0; d < dim; ++d)
                mid[ie][d] = sol.GetVal(e[0] + d, 0) + 0.5 * (sol.GetVal(e[1] + d, 0) + sol.GetVal(e[2] + d, 0));
        }
    }
    for (int64_t iel = 0; iel < fMPhys->NElements(); ++iel) {
        TPZCompEl *cel = fMPhys->Element(iel);
        auto *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel || !cel->Material()) continue;
        auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(cel->Material());
        if (!bc || bc->HasForcingFunctionBC()) continue;
        const int type = bc->Type();
        const bool dirU = TPZMatPoroElastoPlasticUPBase::IsDirichletU(type);
        const bool dirP = TPZMatPoroElastoPlasticUPBase::IsDirichletP(type);
        if (!dirU && !dirP) continue;
        const TPZVec<STATE> &v2 = bc->Val2();
        const int ncorner = cel->Reference()->NCornerNodes();
        const int ncu = mfel->Element(0) ? mfel->Element(0)->NConnects() : 0;
        for (int i = 0; i < cel->NConnects(); ++i) {
            const bool isU = (i < ncu);
            if ((isU && !dirU) || (!isU && !dirP)) continue;
            const int local = isU ? i : i - ncu;
            const bool vertex = local < ncorner;
            TPZConnect &c = cel->Connect(i);
            const int64_t seq = c.SequenceNumber();
            const int64_t pos = block.Position(seq);
            const int size = block.Size(seq);
            if (!isU) {
                for (int k = 0; k < size; ++k) sol(pos + k, 0) = vertex ? v2[0] : 0.;
                continue;
            }
            for (int d = 0; d < dim; ++d) {
                if (type == TPZMatPoroElastoPlasticUPBase::EDirichletUDirectional && bc->Val1().GetVal(d, d) == 0.)
                    continue;
                const STATE val = (d < v2.size()) ? v2[d] : 0.;
                for (int k = d; k < size; k += dim) sol(pos + k, 0) = vertex ? val : 0.;
            }
        }
    }
    if (fNodalNorm) {
        for (size_t ie = 0; ie < fEdgeEquations.size(); ++ie) {
            const auto &e = fEdgeEquations[ie];
            for (int d = 0; d < dim; ++d) {
                if (fConstrained[e[0] + d]) continue;
                sol(e[0] + d, 0) = mid[ie][d] - 0.5 * (sol.GetVal(e[1] + d, 0) + sol.GetVal(e[2] + d, 0));
            }
        }
    }
}

void TPZPoroElastoPlasticUPAnalysis::ToNodalBasis(TPZFMatrix<STATE> &v) const {
    const int dim = fMat->DimensionU();
    for (auto &e : fEdgeEquations) {
        for (int d = 0; d < dim; ++d) {
            const STATE re = v.GetVal(e[0] + d, 0);
            v(e[1] + d, 0) -= 0.5 * re;
            v(e[2] + d, 0) -= 0.5 * re;
        }
    }
}

void TPZPoroElastoPlasticUPAnalysis::ApplyEquationFilter() {
    if (!fEliminateDirichlet || !fStructMatrix) return;
    TPZEquationFilter &filter = fStructMatrix->EquationFilter();
    if (filter.IsActive()) return;
    TPZStack<int64_t> active;
    for (int64_t eq = 0; eq < (int64_t)fConstrained.size(); ++eq)
        if (!fConstrained[eq]) active.Push(eq);
    filter.Reset();
    filter.SetActiveEquations(active);
}

void TPZPoroElastoPlasticUPAnalysis::AssembleFullResidual() {
    if (!fStructMatrix) DebugStop();
    const TPZEquationFilter saved(fStructMatrix->EquationFilter());
    fStructMatrix->EquationFilter().Reset();
    AssembleResidual();
    fStructMatrix->EquationFilter() = saved;
}

void TPZPoroElastoPlasticUPAnalysis::ExternalForces(TPZFMatrix<STATE> &fext) {
    fMat->SetAssembleExternalForcesOnly(true);
    AssembleFullResidual();
    fMat->SetAssembleExternalForcesOnly(false);
    fext = Rhs();
}

bool TPZPoroElastoPlasticUPAnalysis::SolveStep(const TLoadState &s0, const TLoadState &s1,
                                               const TPZFMatrix<STATE> &guess, std::vector<REAL> &residuals) {
    if (!fEquationsIdentified) IdentifyEquations();
    ApplyEquationFilter();
    residuals.clear();
    const int64_t neq = fMPhys->NEquations();
    fMat->SetTimeStep(s1.fTime - s0.fTime);
    fMat->SetLoadFactor(s1.fLambda);
    ApplyControlledDisplacement(s1.fUc);

    TPZFMatrix<STATE> fext;
    ExternalForces(fext);
    if (fNodalNorm) ToNodalBasis(fext);
    const REAL ref = std::max(Norm(fext), REAL(1.));

    // initial guess: displacement from guess, pressure of the converged step, prescribed values
    TPZFMatrix<STATE> sol(fConverged);
    for (int64_t eq = 0; eq < neq; ++eq)
        if (!fIsPressure[eq]) sol(eq, 0) = guess.GetVal(eq, 0);
    ImposeDirichletValues(sol);

    for (int it = 1; it <= fMaxIt; ++it) {
        LoadSolution(sol);
        fMat->ResetFailedProjections();
        Assemble();
        if (fMat->NFailedProjections() > 0) {
            if (fVerbose > 1) std::cout << "      local projection failed at " << fMat->NFailedProjections()
                                        << " integration points" << std::endl;
            return false;
        }
        TPZFMatrix<STATE> rhs(Rhs());
        if (fNodalNorm) ToNodalBasis(rhs);
        REAL norm2 = 0.;
        for (int64_t eq = 0; eq < neq; ++eq)
            if (!fConstrained[eq]) norm2 += rhs.GetVal(eq, 0) * rhs.GetVal(eq, 0);
        const REAL nr = std::sqrt(norm2) / ref;
        residuals.push_back(nr);
        if (fVerbose > 1) std::cout << "      iter " << it << "  |R|/|F| = " << std::scientific << std::setprecision(3)
                                    << nr << std::defaultfloat << std::endl;
        if (!std::isfinite(nr)) return false;
        if (nr < fTol && it > 1) return true;
        Solve();
        const TPZFMatrix<STATE> &du = Solution();
        for (int64_t eq = 0; eq < neq; ++eq) sol(eq, 0) += du.GetVal(eq, 0);
    }
    return false;
}

void TPZPoroElastoPlasticUPAnalysis::AcceptSolution() {
    // assembling the residual with fUpdateMem set stores the state of the converged solution at each point
    auto *withmem = dynamic_cast<TPZMatWithMemBase *>(fMPhys->FindMaterial(fMat->MaterialId()));
    if (!withmem) DebugStop();
    withmem->SetUpdateMem(true);
    AssembleResidual();
    withmem->SetUpdateMem(false);
}

bool TPZPoroElastoPlasticUPAnalysis::AdvanceStep(const TLoadState &s0, const TLoadState &s1,
                                                 const TPZFMatrix<STATE> &dupred, int level, int &nbisect) {
    if (!fEquationsIdentified) IdentifyEquations();
    const int64_t neq = fMPhys->NEquations();
    // initial guess of the displacement: converged solution plus the predicted increment
    TPZFMatrix<STATE> guess(fConverged);
    for (int64_t eq = 0; eq < neq; ++eq)
        if (!fIsPressure[eq]) guess(eq, 0) += dupred.GetVal(eq, 0);
    std::vector<REAL> res;
    if (SolveStep(s0, s1, guess, res)) {
        AcceptSolution();
        fConverged = fMPhys->Solution();
        fCurrent = s1;
        TStepLog log;
        log.fState = s1;
        log.fResiduals = res;
        log.fLevel = level;
        fStepLog.push_back(log);
        return true;
    }
    LoadSolution(fConverged);
    if (level >= fMaxBisections) return false;
    nbisect++;
    TLoadState mid((s0.fTime + s1.fTime) / 2., (s0.fLambda + s1.fLambda) / 2., (s0.fUc + s1.fUc) / 2.);
    if (fVerbose > 0) std::cout << "    bisection level " << level + 1 << std::endl;
    TPZFMatrix<STATE> half(dupred);
    half *= 0.5;
    const TPZFMatrix<STATE> un(fConverged);
    if (!AdvanceStep(s0, mid, half, level + 1, nbisect)) return false;
    TPZFMatrix<STATE> dmid(fConverged);
    dmid -= un;
    return AdvanceStep(mid, s1, dmid, level + 1, nbisect);
}

bool TPZPoroElastoPlasticUPAnalysis::Run(const std::vector<TLoadState> &steps,
                                         const std::function<void(int, const TLoadState &)> &monitor,
                                         const TLoadState &start) {
    if (!fEquationsIdentified) IdentifyEquations();
    fConverged = fMPhys->Solution();
    fCurrent = start;
    LoadSolution(fConverged);
    if (monitor) monitor(0, start);
    const int64_t neq = fMPhys->NEquations();
    TPZFMatrix<STATE> du(neq, 1, 0.);
    TLoadState prev = start;
    for (size_t k = 0; k < steps.size(); ++k) {
        TPZFMatrix<STATE> dupred(neq, 1, 0.);
        if (fPredictor && k > 0) dupred = du;
        const TPZFMatrix<STATE> un(fConverged);
        const size_t nlog0 = fStepLog.size();
        int nbisect = 0;
        if (!AdvanceStep(prev, steps[k], dupred, 0, nbisect)) {
            std::cout << "TPZPoroElastoPlasticUPAnalysis: increment " << k + 1 << " failed" << std::endl;
            return false;
        }
        du = fConverged;
        du -= un;
        prev = steps[k];
        if (fVerbose > 0) {
            int its = 0;
            for (size_t l = nlog0; l < fStepLog.size(); ++l) its += fStepLog[l].fResiduals.size();
            std::cout << "  step " << k + 1 << "/" << steps.size() << ": t = " << steps[k].fTime
                      << "  lambda = " << steps[k].fLambda << "  uc = " << steps[k].fUc << "  evaluations = " << its
                      << "  bisections = " << nbisect << std::endl;
        }
        if (monitor) monitor(int(k + 1), steps[k]);
    }
    return true;
}

void TPZPoroElastoPlasticUPAnalysis::ComputeReactionVector(TPZFMatrix<STATE> &Ru) {
    // assemble all the materials except the Dirichlet conditions
    std::set<int> ids;
    for (auto &it : fMPhys->MaterialVec()) {
        auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(it.second);
        if (bc && (TPZMatPoroElastoPlasticUPBase::IsDirichletU(bc->Type()) ||
                   TPZMatPoroElastoPlasticUPBase::IsDirichletP(bc->Type())))
            continue;
        ids.insert(it.first);
    }
    const std::set<int> previous = fStructMatrix->MaterialIds();
    fStructMatrix->SetMaterialIds(ids);
    AssembleFullResidual();
    fStructMatrix->SetMaterialIds(previous);
    Ru = Rhs();
    Ru *= -1.;
}

std::set<int64_t> TPZPoroElastoPlasticUPAnalysis::NodesOfMaterials(const std::set<int> &matids) const {
    std::set<int64_t> nodes;
    TPZGeoMesh *gmesh = fMPhys->Reference();
    for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel || gel->HasSubElement() || !matids.count(gel->MaterialId())) continue;
        for (int i = 0; i < gel->NCornerNodes(); ++i) nodes.insert(gel->NodeIndex(i));
    }
    return nodes;
}

int64_t TPZPoroElastoPlasticUPAnalysis::NodeEquation(int64_t node, int space, int component) const {
    auto it = fNodeConnects.find(node);
    if (it == fNodeConnects.end()) return -1;
    const int64_t ic = (space == 0) ? it->second.first : it->second.second;
    if (ic < 0) return -1;
    const TPZConnect &c = fMPhys->ConnectVec()[ic];
    const int64_t seq = c.SequenceNumber();
    if (seq < 0) return -1;
    return fMPhys->Block().Position(seq) + component;
}

REAL TPZPoroElastoPlasticUPAnalysis::NodalValue(int64_t node, int space, int component) const {
    const int64_t eq = NodeEquation(node, space, component);
    if (eq < 0) DebugStop();
    TPZFMatrix<STATE> &sol = fMPhys->Solution();
    return sol.GetVal(eq, 0);
}

REAL TPZPoroElastoPlasticUPAnalysis::Reaction(const std::set<int> &bcids, int component,
                                              const std::set<int> &exclude) {
    if (!fEquationsIdentified) IdentifyEquations();
    // Nodes of the (serendipity) Lagrange basis: vertices and mid-edge nodes of the boundary elements.
    // The reaction of a set S of Lagrange nodes is r(v), v = sum_{n in S} N_n e_d. With the hierarchical
    // basis of NeoPZ (vertex functions phi_c, edge functions phi_e equal to 1 at the mid-edge node),
    // N_e = phi_e and N_c = phi_c - 1/2 sum_{e containing c} phi_e, so that v has coefficient [c in S] on
    // the vertex functions and [e in S] - ([c1 in S] + [c2 in S])/2 on the edge functions.
    TPZGeoMesh *gmesh = fMPhys->Reference();
    auto edgesOf = [&](const std::set<int> &ids, std::set<std::pair<int64_t, int64_t>> &edges,
                       std::set<int64_t> &nodes) {
        for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel || gel->HasSubElement() || !ids.count(gel->MaterialId())) continue;
            for (int i = 0; i < gel->NCornerNodes(); ++i) nodes.insert(gel->NodeIndex(i));
            for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side) {
                if (gel->SideDimension(side) != 1) continue;
                int64_t a = gel->SideNodeIndex(side, 0), b = gel->SideNodeIndex(side, 1);
                edges.insert({std::min(a, b), std::max(a, b)});
            }
        }
    };
    std::set<std::pair<int64_t, int64_t>> edgesS, edgesEx;
    std::set<int64_t> nodesS, nodesEx;
    edgesOf(bcids, edgesS, nodesS);
    edgesOf(exclude, edgesEx, nodesEx);
    for (auto n : nodesEx) nodesS.erase(n);
    for (auto &e : edgesEx) edgesS.erase(e);

    TPZFMatrix<STATE> Ru;
    ComputeReactionVector(Ru);
    REAL sum = 0.;
    for (auto n : nodesS) {
        const int64_t eq = NodeEquation(n, 0, component);
        if (eq >= 0) sum += Ru.GetVal(eq, 0);
    }
    // edge functions of the displacement mesh
    TPZCompMesh *umesh = fMPhys->MeshVector()[0];
    const int dim = fMat->DimensionU();
    std::set<std::pair<int64_t, int64_t>> visited;
    for (int64_t iel = 0; iel < umesh->NElements(); ++iel) {
        auto *intel = dynamic_cast<TPZInterpolatedElement *>(umesh->Element(iel));
        if (!intel || !intel->Reference()) continue;
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side) {
            if (gel->SideDimension(side) != 1) continue;
            int64_t a = gel->SideNodeIndex(side, 0), b = gel->SideNodeIndex(side, 1);
            std::pair<int64_t, int64_t> e(std::min(a, b), std::max(a, b));
            if (visited.count(e)) continue;
            visited.insert(e);
            const REAL coef = (edgesS.count(e) ? 1. : 0.) - 0.5 * ((nodesS.count(a) ? 1. : 0.) + (nodesS.count(b) ? 1. : 0.));
            if (coef == 0.) continue;
            const int64_t ic = intel->ConnectIndex(intel->MidSideConnectLocId(side));
            TPZConnect &c = fMPhys->ConnectVec()[ic];
            const int64_t seq = c.SequenceNumber();
            if (seq < 0) continue;
            const int64_t pos = fMPhys->Block().Position(seq);
            const int size = fMPhys->Block().Size(seq);
            for (int k = component; k < size; k += dim) sum += coef * Ru.GetVal(pos + k, 0);
        }
    }
    return sum;
}
