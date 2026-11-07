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

namespace {
bool gUseSkyline = false; // a single copy, inside libpz
}

TPZUnsatAnalysis::TPZUnsatAnalysis(TPZMultiphysicsCompMesh *mesh, const std::map<int, std::array<bool, 2>> &fixU)
    : TPZLinearAnalysis(mesh, true), fMPhys(mesh) {
    for (auto &m : mesh->MaterialVec())
        if (m.first > 0)
            if (auto *u = dynamic_cast<TMat *>(m.second)) fMats.push_back(u);
    if (fMats.empty()) DebugStop();
    // libpz itself is compiled with USING_MKL; PZ_USING_MKL is only given to the targets that link to it, so testing
    // PZ_USING_MKL alone here always selected the skyline matrix
#if defined(PZ_USING_MKL) || defined(USING_MKL)
    if (!gUseSkyline) {
        TPZSpStructMatrix<STATE> smat(mesh);
        SetStructuralMatrix(smat);
    } else
#endif
    {
        TPZSkylineNSymStructMatrix<STATE> smat(mesh);
        SetStructuralMatrix(smat);
    }
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    SetSolver(step);
    fNEq = mesh->NEquations();
    IdentifyEquations(fixU);
    fX = mesh->Solution();
}

void TPZUnsatAnalysis::SetUseSkyline(bool skyline) { gUseSkyline = skyline; }

bool TPZUnsatAnalysis::UsesPardiso() {
#if defined(PZ_USING_MKL) || defined(USING_MKL)
    return !gUseSkyline;
#else
    return false;
#endif
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
        for (int64_t i = 0; i < mem->NElements(); ++i) {
            // the properties of the point (random fields: strength {c, phi, psi} and the conductivity multiplier)
            // survive; the state (stresses, strains, pressures) goes back to the default
            TPZUnsaturatedMem &it = (*mem)[i];
            const TPZVec<REAL> prop = it.m_elastoplastic_state.fmatprop;
            const REAL kmult = it.fKMult;
            it = def;
            it.m_elastoplastic_state.fmatprop = prop;
            it.fKMult = kmult;
        }
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

TPZFMatrix<STATE> TPZUnsatAnalysis::LoadVector(const TPZFMatrix<STATE> &x) {
    SetGravityFactor(0.);
    TPZFMatrix<STATE> fbar = Residual(x, 0., true);
    SetGravityFactor(1.);
    fbar -= Residual(x, 0., true);
    return fbar;
}

TPZUnsatAnalysis::ArcLengthResult TPZUnsatAnalysis::ArcLengthStep(const std::map<int64_t, REAL> &fixed,
                                                                  const TPZFMatrix<STATE> &fbar, REAL lambdan, REAL l,
                                                                  const TPZFMatrix<STATE> &dun, REAL dlambdan,
                                                                  const ArcLengthOptions &options) {
    REAL fnorm = 0.;
    for (int64_t i = 0; i < fNEq; ++i)
        if (!fixed.count(i)) fnorm += fbar.GetVal(i, 0) * fbar.GetVal(i, 0);
    // the factor does not move the yield surface; for the predictor the surface itself is shrunk a little
    return ArcLengthIterate(fixed, lambdan, l, dun, dlambdan, options, std::sqrt(fnorm), 0.,
                            [this](REAL lam, bool predictor) {
                                SetGravityFactor(lam);
                                SetStrengthReductionFactor(predictor ? 1. + 1.e-6 : 1.);
                            },
                            [&fbar](REAL, const TPZFMatrix<STATE> &, TPZFMatrix<STATE> &f) { f = fbar; });
}

TPZUnsatAnalysis::ArcLengthResult TPZUnsatAnalysis::ArcLengthStepSRM(const std::map<int64_t, REAL> &fixed, REAL Fn, REAL l,
                                                                     const TPZFMatrix<STATE> &dun, REAL dFn,
                                                                     const ArcLengthOptions &options) {
    // reference of the residual: the weight, as in Solve
    LoadSolution(fX);
    for (auto *m : fMats) m->SetAssembleGravityOnly(true);
    AssembleFullResidual();
    for (auto *m : fMats) m->SetAssembleGravityOnly(false);
    REAL fref = 0.;
    {
        TPZFMatrix<STATE> &g = Rhs();
        for (int64_t i = 0; i < fNEq; ++i)
            if (!fIsP[i]) fref += g.GetVal(i, 0) * g.GetVal(i, 0);
    }
    fref = std::max(std::sqrt(fref), fForceRef);
    return ArcLengthIterate(fixed, Fn, l, dun, dFn, options, fref, 1.e-6, [this](REAL F, bool) { SetStrengthReductionFactor(F); },
                            [this](REAL F, const TPZFMatrix<STATE> &rhsF, TPZFMatrix<STATE> &f) {
                                // rhs = -R, so -dR/dF = (rhs(F + h) - rhs(F)) / h
                                const REAL h = 1.e-6 * F;
                                SetStrengthReductionFactor(F + h);
                                AssembleFullResidual();
                                SetStrengthReductionFactor(F);
                                f = Rhs();
                                f -= rhsF;
                                f *= 1. / h;
                            });
}

TPZUnsatAnalysis::ArcLengthResult TPZUnsatAnalysis::ArcLengthIterate(
    const std::map<int64_t, REAL> &fixed, REAL lambdan, REAL l, const TPZFMatrix<STATE> &dun, REAL dlambdan,
    const ArcLengthOptions &options, REAL fref, REAL predictorShift, const std::function<void(REAL, bool)> &setFactor,
    const std::function<void(REAL, const TPZFMatrix<STATE> &, TPZFMatrix<STATE> &)> &tangentLoad) {
    ArcLengthResult res;
    res.lambda = lambdan;
    res.du.Redim(fNEq, 1);
    for (auto *m : fMats) {
        m->SetTimeStep(0.);
        m->SetFrozenPressures(true);
    }
    SetFilter(fixed);
    std::vector<int64_t> free;
    for (int64_t i = 0; i < fNEq; ++i)
        if (!fixed.count(i)) free.push_back(i);
    auto dot = [&free](const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) {
        REAL s = 0.;
        for (const int64_t i : free) s += a.GetVal(i, 0) * b.GetVal(i, 0);
        return s;
    };
    if (!(fref > 0.)) return res;

    TPZFMatrix<STATE> x(fX);
    for (auto &f : fixed) x(f.first, 0) = f.second;
    TPZFMatrix<STATE> &Du = res.du; // (i) Delta u^(0) = 0, lambda^(0) = lambda_n
    REAL lambda = lambdan;
    const REAL alpha2 = options.alpha * options.alpha;
    const int maxit = options.maxit > 0 ? options.maxit : fMaxItFrozen;
    const bool hasPrevious = dun.Rows() == fNEq;
    int kfirst = 1;
    if (options.secantPredictor && hasPrevious) {
        // the previous increment scaled to the length l; not extrapolated by more than 100 times
        const REAL previous2 = dot(dun, dun) + alpha2 * dlambdan * dlambdan;
        if (previous2 > 0. && l * l < 1.e4 * previous2) {
            const REAL scale = l / std::sqrt(previous2);
            for (const int64_t i : free) {
                x(i, 0) += scale * dun.GetVal(i, 0);
                Du(i, 0) = scale * dun.GetVal(i, 0);
            }
            lambda += scale * dlambdan;
            res.lambda = lambda;
            kfirst = 2;
        }
    }
    for (int k = kfirst; k <= maxit; ++k) {
        res.iter = k - kfirst + 1;
        // a factor that is not positive has no meaning, and F <= 0 stops the program in ChangeLocalMatParameters
        if (!(lambda > 0.)) return res;
        // (ii) K_T(u^(k-1)); the assembled right-hand side is -r^(k-1)
        const bool predictor = k == 1 && options.plasticPredictor;
        const REAL lambdak = predictor ? lambda * (1. + predictorShift) : lambda;
        setFactor(lambdak, predictor);
        LoadSolution(x);
        AssembleSystem();
        const TPZFMatrix<STATE> rhsk = Rhs();
        // (ix) of the previous iteration: at least one correction after the predictor
        res.eu = std::sqrt(dot(rhsk, rhsk)) / fref;
        if (fVerbose)
            std::cout << "      arc it " << k << " lambda = " << lambdak << " |r|/ref = " << res.eu
                      << " |Du| = " << std::sqrt(dot(Du, Du)) << " l = " << l << "\n";
        if (!std::isfinite(res.eu) || res.eu > 1.e6) return res;
        if (k > 1 && res.eu < fTolU) {
            res.conv = true;
            return res;
        }
        // the vector of the tangential solution, while the mesh holds the iterate: a linear solve loads its result in the
        // mesh, and a residual assembled after it would be evaluated at the correction
        TPZFMatrix<STATE> fbar;
        tangentLoad(lambdak, rhsk, fbar);
        // (iii) K_T dustar = -r and K_T dubar = fbar, with the same factorization. A singular tangent makes the LU stop
        // with an exception (DebugStop): the increment is then reported as not converged
        TPZFMatrix<STATE> dustar, dubar;
        try {
            TPZFMatrix<STATE> &rhs = Rhs();
            rhs = rhsk;
            TPZLinearAnalysis::Solve();
            dustar = Solution();
            rhs = fbar;
            TPZLinearAnalysis::Solve();
            dubar = Solution();
        } catch (...) {
            return res;
        }
        // (iv) iterative load factor
        const REAL a = dot(dubar, dubar) + alpha2;
        if (a == 0.) res.vertical = true;
        if (!(a > 0.) || !std::isfinite(a)) return res;
        const REAL Dlambda = lambda - lambdan;
        REAL dlambda;
        if (k == 1) {
            // predictor, Eq. (4.120), with the sign of the secant criterion (4.123)
            const REAL sign = hasPrevious && dot(dun, dubar) + alpha2 * dlambdan < 0. ? -1. : 1.;
            dlambda = sign * l / std::sqrt(a);
        } else {
            // constraint (4.116)-(4.117), with w = Delta u^(k-1) + dustar and the term of alpha
            TPZFMatrix<STATE> w(Du);
            w += dustar;
            const REAL b = 2. * (dot(w, dubar) + alpha2 * Dlambda), c = dot(w, w) + alpha2 * Dlambda * Dlambda - l * l;
            const REAL disc = b * b - 4. * a * c;
            if (disc < 0.) {
                res.noRoot = true;
                return res;
            }
            // the two roots without cancellation
            const REAL q = -0.5 * (b + (b >= 0. ? 1. : -1.) * std::sqrt(disc));
            const REAL r1 = q / a, r2 = q != 0. ? c / q : r1;
            // Eq. (4.118): the root with the smallest angle between Delta u^(k-1) and Delta u^(k)
            const REAL direction = dot(dubar, Du) + alpha2 * Dlambda;
            dlambda = r1 * direction >= r2 * direction ? r1 : r2;
        }
        // (v)-(vii) corrections of the load factor and of the total and incremental displacements
        lambda += dlambda;
        res.lambda = lambda;
        for (const int64_t i : free) {
            const REAL d = dustar.GetVal(i, 0) + dlambda * dubar.GetVal(i, 0);
            x(i, 0) += d;
            Du(i, 0) += d;
        }
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
