// SlopeStability.cpp — ver SlopeStability.h

#include "SlopeStability.h"

#include <cmath>
#include <iostream>
#include <set>
#include <sstream>

#include "TPZBndCondT.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzcompelwithmem.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzpostprocanalysis.h"
#include "pzstepsolver.h"

// ============================================================== material
template <class T>
int TPZMatElastoPlastic2DSeepage<T>::VariableIndex(const std::string &name) const {
    if (name == "Cohesion") return ECohesion;
    if (name == "FrictionAngle") return EFriction;
    if (name == "ExcessPorePressure") return EExcessPorePressure;
    if (name == "SeepageForce") return ESeepageForce;
    if (name == "PlasticStrainNorm") return EPlasticStrainNorm;
    return TBase::VariableIndex(name);
}

template <class T>
int TPZMatElastoPlastic2DSeepage<T>::NSolutionVariables(int var) const {
    switch (var) {
        case ECohesion:
        case EFriction:
        case EExcessPorePressure:
        case EPlasticStrainNorm:
            return 1;
        case ESeepageForce:
            return 3;
        default:
            return TBase::NSolutionVariables(var);
    }
}

template <class T>
void TPZMatElastoPlastic2DSeepage<T>::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) {
    if (var < ECohesion || var > EPlasticStrainNorm) {
        TBase::Solution(data, var, Solout);
        return;
    }
    const TPZElastoPlasticMem &mem = this->MemItem(data.intGlobPtIndex);
    const TPZVec<REAL> &mp = mem.m_elastoplastic_state.fmatprop;
    switch (var) {
        case ECohesion:
            Solout[0] = mp.size() > 0 ? mp[0] : 0.;
            break;
        case EFriction:
            Solout[0] = mp.size() > 1 ? mp[1] * 180. / M_PI : 0.;
            break;
        case EExcessPorePressure:
            Solout[0] = mem.fPorePressure;
            break;
        case ESeepageForce:
            Solout[0] = mem.fdPorePressure.size() > 0 ? -mem.fdPorePressure[0] : 0.;
            Solout[1] = mem.fdPorePressure.size() > 1 ? -mem.fdPorePressure[1] : 0.;
            Solout[2] = 0.;
            break;
        case EPlasticStrainNorm: {
            const TPZTensor<REAL> &ep = mem.m_elastoplastic_state.m_eps_p;
            REAL s = 0.;
            for (int i = 0; i < 6; i++) s += ep[i] * ep[i];
            Solout[0] = std::sqrt(s);
        } break;
    }
}

// ============================================================== modelos
namespace {

void SetupModel(TMohrCoulomb &mc, const TSoil &soil) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(soil.E, soil.nu);
    const REAL phi = soil.phiDeg * M_PI / 180.;
    TPZYCMohrCoulombPV2 yc(phi, phi, soil.c, ER);
    mc.SetPlasticCriterion(yc);
    mc.SetElasticResponse(ER);
    mc.SetStrengthReductionFactor(1.);
}

void SetupModel(TPZModifiedCamClay &mcc, const TSoil &soil) {
    const REAL phi = soil.phiDeg * M_PI / 180.;
    const REAL M = TPZModifiedCamClay::MFromFriction(phi, soil.mapping);
    const REAL pt = soil.c / std::tan(phi);
    const REAL G = soil.E / (2. * (1. + soil.nu)), K = soil.E / (3. * (1. - 2. * soil.nu));
    mcc.SetUp(M, soil.lambda, soil.kappa, 0., soil.v0, 100., 50., pt, 1., TPZModifiedCamClay::ELinear,
              TPZModifiedCamClay::EConstantG, G, soil.nu);
    mcc.SetLinearBulkModulus(K);
    mcc.SetStrengthMapping(soil.mapping);
    mcc.SetStrengthReductionFactor(1.);
}

template <class T>
constexpr bool IsCamClay() {
    return std::is_same_v<T, TPZModifiedCamClay>;
}

} // namespace

// ============================================================== TSlopeFEM
template <class T>
TSlopeFEM<T>::TSlopeFEM(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TSoil &soil, const TSolverOptions &opt)
    : fGMesh(gmesh), fGeo(geo), fSoil(soil), fOpt(opt) {
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(2);
    fCMesh->SetDefaultOrder(opt.porder);
    fCMesh->SetAllCreateFunctionsContinuousWithMem();
    fMat = new TPZMatElastoPlastic2DSeepage<T>(ESoil);
    T model;
    SetupModel(model, soil);
    fMat->SetPlasticityModel(model);
    TPZManVector<REAL, 3> f = {0., -soil.BodyForce(), 0.};
    fMat->SetBodyForce0(f);
    fMat->SetBodyForce(f);
    fCMesh->InsertMaterialObject(fMat);
    TPZFNMatrix<4, STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> fixed = {1., 1.}, roller = {1., 0.};
    fCMesh->InsertMaterialObject(fMat->CreateBC(fMat, EBottom, 3, val1, fixed));
    fCMesh->InsertMaterialObject(fMat->CreateBC(fMat, ELeft, 3, val1, roller));
    fCMesh->InsertMaterialObject(fMat->CreateBC(fMat, ERight, 3, val1, roller));
    std::set<int> matids = {ESoil, EBottom, ELeft, ERight};
    fCMesh->AutoBuild(matids);

    fAn = std::make_unique<TPZLinearAnalysis>(fCMesh, true);
    TPZStepSolver<STATE> step;
    if constexpr (IsCamClay<T>()) {
        TPZSkylineNSymStructMatrix<STATE> skl(fCMesh);
        skl.SetNumThreads(0);
        fAn->SetStructuralMatrix(skl);
        step.SetDirect(ELU);
    } else {
        TPZSkylineStructMatrix<STATE> skl(fCMesh);
        skl.SetNumThreads(0);
        fAn->SetStructuralMatrix(skl);
        step.SetDirect(ELDLt);
    }
    fAn->SetSolver(step);
    BuildPoints();
    if constexpr (IsCamClay<T>()) {
        // σ0 isotrópica provisória até SetInitialStress
        fSigma0.assign(fPoints.size(), TPZTensor<REAL>());
    }
    SetUniformStrength(soil.c, soil.phiDeg * M_PI / 180.);
}

template <class T>
TSlopeFEM<T>::~TSlopeFEM() {
    delete fPost;
    fAn.reset();
    delete fCMesh;
}

template <class T>
void TSlopeFEM<T>::BuildPoints() {
    const int64_t nmem = fMat->GetMemory()->NElements();
    fPoints.assign(nmem, TPoint{-1, TPZManVector<REAL, 3>(2, 0.), TPZManVector<REAL, 3>(3, 0.)});
    for (TPZCompEl *cel : fCMesh->ElementVec()) {
        auto *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel || !cel->Reference() || cel->Reference()->Dimension() != 2) continue;
        if (cel->Material() != fMat) continue;
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        TPZMaterialDataT<STATE> data;
        intel->InitMaterialData(data);
        TPZManVector<REAL, 3> qsi(2, 0.);
        REAL w;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            rule.Point(ip, qsi, w);
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData(data, qsi);
            const int idx = data.intGlobPtIndex;
            fPoints[idx].gel = cel->Reference()->Index();
            fPoints[idx].qsi = qsi;
            fPoints[idx].x = data.x;
        }
    }
}

template <class T>
void TSlopeFEM<T>::SetStrength(const std::vector<REAL> &c, const std::vector<REAL> &phi) {
    auto &memory = *fMat->GetMemory();
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        if (fPoints[i].gel < 0) continue;
        TPZVec<REAL> &mp = memory[i].m_elastoplastic_state.fmatprop;
        if constexpr (IsCamClay<T>()) {
            if (mp.size() < 9) {
                mp.Resize(9);
                mp.Fill(0.);
            }
            mp[0] = c[i];
            mp[1] = phi[i];
        } else {
            mp.Resize(3);
            mp[0] = c[i];
            mp[1] = phi[i];
            mp[2] = phi[i];
        }
        memory[i].m_elastoplastic_state.fmatpropinit = mp;
    }
    if constexpr (IsCamClay<T>()) {
        // p_c0 depende de M e p_t locais: recalcula se σ0 já foi definida
        bool hasSigma0 = false;
        for (auto &s : fSigma0)
            if (s.I1() != 0.) hasSigma0 = true;
        if (hasSigma0) SetInitialStress(fSigma0);
    }
}

template <class T>
void TSlopeFEM<T>::SetUniformStrength(REAL c, REAL phi) {
    std::vector<REAL> cv(fPoints.size(), c), pv(fPoints.size(), phi);
    SetStrength(cv, pv);
}

template <class T>
void TSlopeFEM<T>::SetSeepage(const std::vector<REAL> &u, const std::vector<TPZManVector<REAL, 2>> &gradu) {
    auto &memory = *fMat->GetMemory();
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        if (fPoints[i].gel < 0) continue;
        if (u.empty()) {
            memory[i].fPorePressure = 0.;
            memory[i].fdPorePressure.Resize(0);
        } else {
            memory[i].fPorePressure = u[i];
            memory[i].fdPorePressure.Resize(2);
            memory[i].fdPorePressure[0] = gradu[i][0];
            memory[i].fdPorePressure[1] = gradu[i][1];
        }
    }
}

template <class T>
void TSlopeFEM<T>::SetInitialStress(const std::vector<TPZTensor<REAL>> &sigma0) {
    if constexpr (!IsCamClay<T>()) {
        (void)sigma0;
        return;
    } else {
        fSigma0 = sigma0;
        auto &memory = *fMat->GetMemory();
        const REAL pmin = 1.;  // kPa: σ0 com p' > -1 kPa (tração) é deslocada para p' = -1 kPa
        for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
            if (fPoints[i].gel < 0) continue;
            TPZVec<REAL> &mp = memory[i].m_elastoplastic_state.fmatprop;
            TPZTensor<REAL> s = sigma0[i];
            REAL p, q;
            TPZModifiedCamClay::Invariants(s, p, q);
            if (p > -pmin) {
                const REAL dp = -pmin - p;
                s[_XX_] += dp;
                s[_YY_] += dp;
                s[_ZZ_] += dp;
                p = -pmin;
            }
            const REAL M = TPZModifiedCamClay::MFromFriction(mp[1], fSoil.mapping);
            const REAL pt = mp[0] / std::tan(mp[1]);
            // superfície (β = 1) que passa por σ0: (p - pt + a)² + (q/M)² = a²
            const REAL d = pt - p;  // > 0
            const REAL a = (d * d + (q / M) * (q / M)) / (2. * d);
            const REAL pc = std::max(2. * a - pt, REAL(1.e-3));
            mp[2] = fSoil.OCR * pc;
            for (int k = 0; k < 6; k++) mp[3 + k] = s[k];
            memory[i].m_elastoplastic_state.fmatpropinit = mp;
            memory[i].m_sigma = s;
            fSigma0[i] = s;
        }
    }
}

template <class T>
void TSlopeFEM<T>::Stresses(std::vector<TPZTensor<REAL>> &sigma) const {
    auto &memory = *fMat->GetMemory();
    sigma.resize(fPoints.size());
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) sigma[i] = memory[i].m_sigma;
}

template <class T>
void TSlopeFEM<T>::ResetState() {
    auto &memory = *fMat->GetMemory();
    for (int64_t i = 0; i < memory.NElements(); i++) {
        TPZElastoPlasticMem &m = memory[i];
        m.m_elastoplastic_state.m_eps_t.Zero();
        m.m_elastoplastic_state.m_eps_p.Zero();
        m.m_elastoplastic_state.m_hardening = 0.;
        m.m_elastoplastic_state.m_m_type = 0;
        if constexpr (IsCamClay<T>()) {
            m.m_sigma = (i < (int64_t)fSigma0.size()) ? fSigma0[i] : TPZTensor<REAL>();
        } else {
            m.m_sigma.Zero();
        }
        m.m_u.Fill(0.);
        m.m_plastic_steps = 0;
    }
    for (auto &kv : fCMesh->MaterialVec()) {
        if (kv.second == fMat) continue;
        if (auto *bcmem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(kv.second)) {
            auto &bm = *bcmem->GetMemory();
            for (int64_t i = 0; i < bm.NElements(); i++) bm[i].m_u.Fill(0.);
        }
    }
    fAn->Solution().Zero();
    fAn->LoadSolution();
    fLambdaG = fLambdaS = 0.;
    fEpsPPrev.clear();
}

template <class T>
void TSlopeFEM<T>::SetLoads(REAL lambdaGravity, REAL lambdaSeepage, REAL F) {
    TPZManVector<REAL, 3> f = {0., -lambdaGravity * fSoil.BodyForce(), 0.};
    fMat->SetBodyForce(f);
    fMat->SetSeepageFactor(lambdaSeepage);
    fMat->GetPlasticModel().SetStrengthReductionFactor(F);
}

template <class T>
void TSlopeFEM<T>::ComputeReference() {
    // ||F_γ'|| nas equações livres: diferença entre os resíduos com e sem a força de corpo (estado atual)
    const REAL lg = fLambdaG, ls = fLambdaS, F = fF;
    SetLoads(1., 0., F);
    fAn->Assemble();
    const int64_t neq = fCMesh->NEquations();
    auto K = fAn->MatrixSolver<STATE>().Matrix();
    fFree.resize(neq);
    for (int64_t i = 0; i < neq; i++) fFree[i] = std::fabs(K->GetVal(i, i)) < 1.e-6 * fMat->BigNumber();
    TPZFMatrix<STATE> r1 = fAn->Rhs();
    SetLoads(0., 0., F);
    fAn->AssembleResidual();
    const TPZFMatrix<STATE> &r0 = fAn->Rhs();
    REAL s = 0.;
    for (int64_t i = 0; i < neq; i++)
        if (fFree[i]) s += std::pow(r1.GetVal(i, 0) - r0.GetVal(i, 0), 2);
    fFRef = std::max(std::sqrt(s), REAL(1.e-12));
    SetLoads(lg, ls, F);
}

template <class T>
void TSlopeFEM<T>::Accept() {
    auto &memory = *fMat->GetMemory();
    fEpsPPrev.resize(fPoints.size());
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) fEpsPPrev[i] = memory[i].m_elastoplastic_state.m_eps_p;
    fMat->SetUpdateMem(true);
    fAn->AssembleResidual();
    fMat->SetUpdateMem(false);
    fAn->Solution().Zero();
    fAn->LoadSolution();
}

template <class T>
bool TSlopeFEM<T>::Solve(REAL lambdaGravity, REAL lambdaSeepage, REAL F, int &iterations) {
    if (fFree.empty()) ComputeReference();
    SetLoads(lambdaGravity, lambdaSeepage, F);
    const int64_t neq = fCMesh->NEquations();
    TPZFMatrix<STATE> x(neq, 1, 0.);
    fAn->LoadSolution(x);
    iterations = 0;
    REAL nr = 0., nr0 = -1.;
    std::vector<REAL> hist;
    bool ok = false;
    try {
        for (int it = 1; it <= fOpt.maxIter; it++) {
            fAn->Assemble();
            const TPZFMatrix<STATE> &rhs = fAn->Rhs();
            REAL s = 0.;
            for (int64_t i = 0; i < neq; i++)
                if (fFree[i]) s += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
            nr = std::sqrt(s) / fFRef;
            if (fOpt.verbose > 1) std::cout << "      it " << it << " |R|/|F| = " << nr << "\n";
            if (!std::isfinite(nr)) break;
            if (nr0 < 0.) nr0 = nr;
            if (nr < fOpt.tol) {
                ok = true;
                break;
            }
            if (it > 3 && nr > 1.e3 * std::max(nr0, REAL(1.))) break;  // divergência
            // estagnação: sem redução de pelo menos 1/2 em 3 iterações (o Newton com tangente consistente
            // converge quadraticamente quando o passo é admissível)
            hist.push_back(nr);
            if (fOpt.stagnation && it >= 6 && nr > 0.5 * hist[hist.size() - 4]) break;
            iterations = it;
            fAn->Solve();
            x += fAn->Solution();
            fAn->LoadSolution(x);
        }
    } catch (std::exception &e) {
        if (fOpt.verbose > 1) std::cout << "      exceção: " << e.what() << "\n";
        ok = false;
    }
    if (ok) {
        Accept();
        fLambdaG = lambdaGravity;
        fLambdaS = lambdaSeepage;
        fF = F;
    } else {
        x.Zero();
        fAn->LoadSolution(x);
        SetLoads(fLambdaG, fLambdaS, fF);
    }
    return ok;
}

template <class T>
TFactorResult TSlopeFEM<T>::Follow(const std::function<void(REAL)> &apply, REAL t0, REAL dt0, REAL tmax) {
    TFactorResult r;
    REAL t = t0, dt = dt0;
    r.factor = t0;
    while (true) {
        REAL tn = std::min(t + dt, tmax);
        apply(tn);
        int its = 0;
        const bool ok = Solve(fLambdaGTarget, fLambdaSTarget, fFTarget, its);
        r.iterations += its;
        if (ok) {
            r.steps++;
            t = tn;
            r.factor = t;
            if (fOpt.verbose) std::cout << "    fator " << t << " convergiu (" << its << " it.)\n";
            if (t >= tmax) {
                r.status = "limite_maximo";
                return r;
            }
            if (its <= fOpt.maxIter / 2) dt *= 1.5;
        } else {
            r.cuts++;
            r.upper = tn;
            r.bracketed = true;
            if (fOpt.verbose) std::cout << "    fator " << tn << " sem convergência\n";
            dt *= 0.5;
            if (dt <= fOpt.relTol * std::max(t, REAL(0.1))) {
                r.status = "ok";
                return r;
            }
        }
    }
}

template <class T>
TFactorResult TSlopeFEM<T>::LoadFactor(REAL lambda0) {
    fFTarget = 1.;
    if constexpr (IsCamClay<T>()) {
        // Cam-Clay: σ0 já equilibra γ' (λ = 1); aplica a percolação com λ = 1 e então aumenta as duas cargas.
        // Colapso antes de λ = 1 só indica Γ < 1 (o caminho não permite reduzir a gravidade abaixo de σ0).
        if (lambda0 < 1.) lambda0 = 1.;
        fLambdaGTarget = 1.;
        auto applySeep = [this](REAL s) { fLambdaSTarget = s; };
        TFactorResult r0 = Follow(applySeep, 0., 1., 1.);
        if (r0.factor < 1.) {
            r0.factor = 0.;
            r0.status = "colapso_antes_de_lambda1";
            return r0;
        }
        auto apply = [this](REAL lam) {
            fLambdaGTarget = lam;
            fLambdaSTarget = lam;
        };
        TFactorResult r = Follow(apply, lambda0, fOpt.step0, fOpt.maxFactor);
        r.iterations += r0.iterations;
        r.steps += r0.steps;
        return r;
    }
    auto apply = [this](REAL lam) {
        fLambdaGTarget = lam;
        fLambdaSTarget = lam;
    };
    return Follow(apply, lambda0, fOpt.step0, fOpt.maxFactor);
}

template <class T>
TFactorResult TSlopeFEM<T>::StrengthReduction(REAL F0) {
    // 1) cargas reais com resistência majorada (F0), a partir do estado atual
    int its = 0;
    TFactorResult r;
    fLambdaGTarget = 1.;
    fLambdaSTarget = 1.;
    fFTarget = F0;
    // aplica a força de corpo e a percolação em incrementos com F0 (Cam-Clay: σ0 já equilibra γ')
    auto applyLoad = [this](REAL lam) {
        if constexpr (IsCamClay<T>()) {
            fLambdaGTarget = 1.;
            fLambdaSTarget = lam;
        } else {
            fLambdaGTarget = lam;
            fLambdaSTarget = lam;
        }
    };
    TFactorResult r0 = Follow(applyLoad, 0., 1., 1.);
    if (r0.factor < 1.) {
        r.status = "instavel_em_F0";
        r.factor = F0;
        r.upper = F0;
        return r;
    }
    (void)its;
    fLambdaGTarget = 1.;
    fLambdaSTarget = 1.;
    auto applyF = [this](REAL F) { fFTarget = F; };
    TFactorResult rf = Follow(applyF, F0, fOpt.step0, fOpt.maxFactor);
    rf.iterations += r0.iterations;
    rf.steps += r0.steps;
    return rf;
}

template <class T>
void TSlopeFEM<T>::PlasticIndicator(std::vector<REAL> &byGel, bool increment) const {
    byGel.assign(fGMesh->NElements(), 0.);
    auto &memory = *fMat->GetMemory();
    const bool inc = increment && fEpsPPrev.size() == fPoints.size();
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        if (fPoints[i].gel < 0) continue;
        const TPZTensor<REAL> &ep = memory[i].m_elastoplastic_state.m_eps_p;
        REAL s = 0.;
        for (int k = 0; k < 6; k++) {
            const REAL d = inc ? ep[k] - fEpsPPrev[i][k] : ep[k];
            s += d * d;
        }
        byGel[fPoints[i].gel] = std::max(byGel[fPoints[i].gel], std::sqrt(s));
    }
}

template <class T>
void TSlopeFEM<T>::DefineVTK(const std::string &file) {
    delete fPost;
    fPost = new TPZPostProcAnalysis();
    fPost->SetCompMesh(fCMesh);
    TPZFStructMatrix<STATE> structmatrix(fPost->Mesh());
    structmatrix.SetNumThreads(0);
    fPost->SetStructuralMatrix(structmatrix);
    TPZStack<std::string> scal, vec, tens, vars;
    for (const char *nm : {"Cohesion", "FrictionAngle", "ExcessPorePressure", "PlasticStrainNorm", "StressXX",
                           "StressYY", "SqrtStressJ2"})
        scal.Push(nm);
    vec.Push("DisplacementDoF");
    vec.Push("SeepageForce");
    for (auto &s : scal) vars.Push(s);
    for (auto &s : vec) vars.Push(s);
    TPZVec<int> matids(1, ESoil);
    fPost->SetPostProcessVariables(matids, vars);
    fPost->DefineGraphMesh(2, scal, vec, file);
}

template <class T>
void TSlopeFEM<T>::WriteVTK(int step) {
    if (!fPost) return;
    fPost->TransferSolution();
    fPost->SetStep(step);
    fPost->PostProcess(1);
}

template class TPZMatElastoPlastic2DSeepage<TMohrCoulomb>;
template class TPZMatElastoPlastic2DSeepage<TPZModifiedCamClay>;
template class TSlopeFEM<TMohrCoulomb>;
template class TSlopeFEM<TPZModifiedCamClay>;
