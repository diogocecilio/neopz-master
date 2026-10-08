// Factor of safety of a plane-strain slope with an elastoplastic Mohr-Coulomb material:
//  - GravityIncrease():   FS = largest multiplier of the self weight with equilibrium;
//  - StrengthReduction(): FS = largest F with equilibrium for c/F and tan(phi)/F (SRM).
// Both use the same path-following driver: Newton with the consistent tangent, accept the
// converged state (plastic memory) and halve the parameter step when Newton fails.
#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H

#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/pzelastoplasticanalysis.h"
#include "TPZLinearAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzcmesh.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzgeoelside.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <iostream>
#include <memory>
#include <set>
#include <thread>
#include <vector>

template <class TPlastic>
class SlopeAnalysis {
public:
    using TMat = TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem>;

    explicit SlopeAnalysis(TPZCompMesh *cmesh) : fCMesh(cmesh) {
        fMat = dynamic_cast<TMat *>(cmesh->FindMaterial(1));
        if (!fMat) DebugStop();
        fGravity = fMat->GetBodyForce();
        fNThreads = std::max(1u, std::thread::hardware_concurrency());
        TPZLinearAnalysis renumber(fCMesh, true); // skyline bandwidth (TPZElastoPlasticAnalysis does not renumber)
    }

    /// FS by gravity increase (load factor on the self weight, strength unchanged)
    REAL GravityIncrease() {
        Reset();
        return Continuation([this](REAL l) { SetGravity(l); }, 0., 0.5, 100., "GI");
    }

    /// FS by shear strength reduction (self weight unchanged)
    REAL StrengthReduction() {
        REAL F0 = 1.;
        for (;;) { // equilibrium under full weight is required before reducing the strength
            Reset();
            SetReduction(F0);
            if (Continuation([this](REAL l) { SetGravity(l); }, 0., 0.5, 1., "SRM-load") >= 1.) break;
            F0 *= 0.5;
            if (F0 < 1.e-3) { std::cout << "[SRM] no equilibrium under self weight\n"; return 0.; }
        }
        return Continuation([this](REAL F) { SetReduction(F); }, F0, 0.25 * F0, 100., "SRM");
    }

    /// Marks the elements with sqrt(J2(eps_p)) >= frac * max in the last accepted state
    void MarkPlasticZone(REAL frac) {
        const TPZVec<REAL> ind = PlasticIndicator();
        const REAL vmax = *std::max_element(ind.begin(), ind.end());
        for (int64_t el = 0; el < ind.size(); el++) if (vmax > 0. && ind[el] >= frac * vmax) fMarked.insert(el);
    }

    /// h-refinement of the marked elements (2:1 balanced, p kept). The state must be Reset() afterwards.
    int64_t Refine() {
        fCMesh->LoadReferences(); // the post-processing mesh takes over the geometric references
        std::vector<TPZCompEl *> todivide;
        for (int64_t el : fMarked) todivide.push_back(fCMesh->ElementVec()[el]);
        fMarked.clear();
        while (!todivide.empty()) { // pointers: freed indices are reused by the sub-elements
            for (TPZCompEl *cel : todivide) {
                auto *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                if (!intel) continue;
                const int p = intel->GetPreferredOrder();
                TPZStack<int64_t> sub;
                intel->Divide(cel->Index(), sub, 0);
                for (int64_t s : sub) dynamic_cast<TPZInterpolationSpace *>(fCMesh->ElementVec()[s])->SetPreferredOrder(p);
            }
            todivide = UnbalancedElements();
        }
        fCMesh->AdjustBoundaryElements();
        fCMesh->CleanUpUnconnectedNodes();
        fCMesh->InitializeBlock();
        fAn.reset();
        TPZLinearAnalysis renumber(fCMesh, true); // new connects are appended: restore a small skyline profile
        return fCMesh->NEquations();
    }

    /// VTK with total displacement and plastic/failure indicators of the last accepted state
    void PostProcess(const std::string &file) {
        TPZPostProcAnalysis pp;
        pp.SetCompMesh(fCMesh);
        TPZManVector<int, 1> matids(1, 1);
        TPZManVector<std::string, 3> scal = {"StrainPlasticJ2", "FailureType"}, vec = {"Displacement"};
        TPZManVector<std::string, 3> vars = {"StrainPlasticJ2", "FailureType", "Displacement"};
        TPZFStructMatrix<STATE> str(pp.Mesh());
        pp.SetStructuralMatrix(str);
        pp.SetPostProcessVariables(matids, vars);
        fCMesh->LoadSolution(fAn->CumulativeSolution()); // "Displacement" reads the mesh solution
        pp.TransferSolution();
        fAn->LoadSolution();                             // back to the zero increment
        pp.DefineGraphMesh(2, scal, vec, file);
        pp.PostProcess(0);
    }

    REAL fTolNewton = 1.e-8; ///< ||R|| <= tol * ||F_ext||
    int fMaxNewton = 30;
    REAL fTolFS = 2.e-3;     ///< stop when the parameter step is below tol * parameter

private:
    TPZCompMesh *fCMesh;
    TMat *fMat;
    std::unique_ptr<TPZElastoPlasticAnalysis> fAn;
    TPZManVector<REAL, 3> fGravity;
    REAL fLambda = 1., fNormF = 1.;
    unsigned fNThreads;
    std::set<int64_t> fMarked; ///< elements to refine

    void SetGravity(REAL lambda) {
        fLambda = lambda;
        TPZManVector<REAL, 3> f(fGravity);
        for (auto &v : f) v *= lambda;
        fMat->SetBodyForce(f);
    }

    void SetReduction(REAL F) { fMat->GetPlasticModel().SetStrengthReductionFactor(F); }

    /// Virgin stress-free state (memory of all materials) and a fresh analysis
    void Reset() {
        for (auto &m : fCMesh->MaterialVec())
            if (auto *mem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(m.second)) mem->ResetMemory();
        fAn = std::make_unique<TPZElastoPlasticAnalysis>(fCMesh, std::cout);
        TPZSkylineStructMatrix<STATE> skl(fCMesh);
        skl.SetNumThreads(fNThreads);
        fAn->SetStructuralMatrix(skl);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        fAn->SetSolver(step);
        SetGravity(1.);
        SetReduction(1.);
        fAn->AssembleResidual(); // stress-free state: residual = external load
        fNormF = Norm(fAn->Rhs());
    }

    /// Newton-Raphson for the increment since the last accepted state, with backtracking on ||R||
    bool Newton(int &it) {
        const int64_t neq = fCMesh->NEquations();
        TPZFMatrix<STATE> u(neq, 1, 0.), du, best;
        fAn->LoadSolution(u);
        const REAL tol = fTolNewton * fNormF * std::max<REAL>(fLambda, 1.e-2);
        fAn->Assemble();
        REAL r = Norm(fAn->Rhs());
        const REAL r0 = r;
        for (it = 1; it <= fMaxNewton; it++) {
            fAn->Solve();
            du = fAn->Solution();
            REAL alpha = 1., rbest = std::numeric_limits<REAL>::max();
            for (int ls = 0; ls < 6; ls++, alpha *= 0.5) { // keep the best trial if none decreases ||R||
                TPZFMatrix<STATE> trial(du);
                trial *= alpha;
                trial += u;
                fAn->LoadSolution(trial);
                fAn->AssembleResidual();
                const REAL rt = Norm(fAn->Rhs());
                if (rt < rbest) { rbest = rt; best = trial; }
                if (rt < r) break;
            }
            u = best;
            r = rbest;
            if (getenv("SLOPE_VERBOSE")) std::cout << "    it " << it << " |R|/|F| " << r / fNormF << "\n";
            if (!std::isfinite(r) || r > 1.e3 * r0) return false; // diverging
            fAn->LoadSolution(u);
            if (r <= tol) return true;
            fAn->Assemble();
        }
        return false;
    }

    /// Increase a parameter p (load factor or strength reduction) from an equilibrium state at p
    REAL Continuation(const std::function<void(REAL)> &set, REAL p, REAL dp, REAL pmax, const char *name) {
        while (p < pmax && dp > fTolFS * std::max<REAL>(p, 1.e-2)) {
            const REAL ptrial = std::min(p + dp, pmax);
            set(ptrial);
            int it = 0;
            const bool ok = Newton(it);
            std::cout << "[" << name << "] p = " << ptrial << (ok ? " converged" : " failed") << " (" << it << " it.)\n";
            if (ok) {
                fAn->AcceptSolution(); // commits the plastic memory
                p = ptrial;
                if (it <= 8) dp *= 1.5;
            } else {
                dp *= 0.5;
            }
        }
        set(p);
        TPZFMatrix<STATE> zero(fCMesh->NEquations(), 1, 0.);
        fAn->LoadSolution(zero); // drop the last (failed) trial: the mesh holds the accepted state
        std::cout << "[" << name << "] last converged p = " << p << (p >= pmax ? " (no collapse up to pmax)" : "") << "\n";
        return p;
    }

    /// Element indicator: max over integration points of sqrt(J2(eps_p)) (engineering shear -> tensor)
    TPZVec<REAL> PlasticIndicator() {
        TPZVec<REAL> ind(fCMesh->NElements(), 0.);
        for (int64_t el = 0; el < ind.size(); el++) {
            TPZCompEl *cel = fCMesh->ElementVec()[el];
            if (!cel || cel->Material() != fMat) continue;
            TPZManVector<int64_t> mem;
            cel->GetMemoryIndices(mem);
            for (int64_t m : mem) {
                if (m < 0) continue;
                TPZTensor<REAL> ep = fMat->MemItem(m).m_elastoplastic_state.m_eps_p;
                ep.XY() *= 0.5; ep.XZ() *= 0.5; ep.YZ() *= 0.5;
                ind[el] = std::max(ind[el], std::sqrt(std::max<REAL>(ep.J2(), 0.)));
            }
        }
        return ind;
    }

    /// Volume elements neighbouring an element more than one level finer
    std::vector<TPZCompEl *> UnbalancedElements() {
        std::set<int64_t> need;
        for (int64_t el = 0; el < fCMesh->NElements(); el++) {
            TPZCompEl *cel = fCMesh->ElementVec()[el];
            if (!cel || !cel->Reference() || cel->Material() != fMat) continue;
            TPZGeoEl *gel = cel->Reference();
            for (int s = 0; s < gel->NSides(); s++) {
                TPZGeoElSide side(gel, s);
                if (side.Dimension() != 1) continue;
                TPZCompElSide big = side.LowerLevelCompElementList2(1);
                if (!big || big.Element()->Material() != fMat) continue;
                if (gel->Level() - big.Element()->Reference()->Level() > 1) need.insert(big.Element()->Index());
            }
        }
        std::vector<TPZCompEl *> cels;
        for (int64_t el : need) cels.push_back(fCMesh->ElementVec()[el]);
        return cels;
    }
};

#endif
