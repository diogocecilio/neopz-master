// CoupledDrawdown.cpp — ver CoupledDrawdown.h

#include "CoupledDrawdown.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

#include "TPZBndCondT.h"
#include "TPZNullMaterial.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzpostprocanalysis.h"
#include "pzstepsolver.h"

namespace {
template <class T>
constexpr bool IsCamClay() {
    return std::is_same_v<T, TPZModifiedCamClay>;
}
} // namespace

template <class T>
TCoupledDrawdown<T>::TCoupledDrawdown(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TSoil &soil,
                                      const TParams &par)
    : fGMesh(gmesh), fGeo(geo), fSoil(soil), fPar(par) {
    fZw = geo.D();
    T model;
    SetupSoilModel(model, soil);
    const int bcids[] = {EBottom, ELeft, ERight, ECrest, EFace, EToe};

    // ---- malha de u (H1 vetorial com memória); o material só serve ao pós-processamento
    fCMeshU = new TPZCompMesh(gmesh);
    fCMeshU->SetDimModel(2);
    fCMeshU->SetDefaultOrder(par.porderU);
    fCMeshU->SetAllCreateFunctionsContinuousWithMem();
    fMatU = new TPZMatElastoPlastic2DSeepage<T>(ESoil);
    fMatU->SetPlasticityModel(model);
    fCMeshU->InsertMaterialObject(fMatU);
    {
        TPZFNMatrix<4, STATE> v1(2, 2, 0.);
        TPZManVector<STATE, 2> v2(2, 0.);
        for (int id : bcids) fCMeshU->InsertMaterialObject(fMatU->CreateBC(fMatU, id, 1, v1, v2));
    }
    fCMeshU->AutoBuild();

    // ---- malha de p (H1 escalar linear)
    fCMeshP = new TPZCompMesh(gmesh);
    fCMeshP->SetDimModel(2);
    fCMeshP->SetDefaultOrder(1);
    fCMeshP->SetAllCreateFunctionsContinuous();
    auto *nullmat = new TPZNullMaterial<STATE>(ESoil, 2, 1);
    fCMeshP->InsertMaterialObject(nullmat);
    {
        TPZFNMatrix<1, STATE> v1(1, 1, 0.);
        TPZManVector<STATE, 1> v2(1, 0.);
        for (int id : bcids) fCMeshP->InsertMaterialObject(nullmat->CreateBC(nullmat, id, 0, v1, v2));
    }
    fCMeshP->AutoBuild();

    // ---- malha multifísica u-p
    fMPhys = new TPZMultiphysicsCompMesh(gmesh);
    fMPhys->SetDimModel(2);
    fMPhys->SetAllCreateFunctionsMultiphysicElem();
    fMat = new TMatUP(ESoil, 2);
    fMat->SetPlasticModel(model);
    fMat->SetAlpha(1.);
    fMat->SetSe(par.Se);
    const REAL gw = soil.gammaW, D = geo.D();
    fMat->SetPermeability(par.alpha * par.kv / gw, par.kv / gw, 0.);  // K/μ = k/γw
    fMat->SetViscosity(1.);
    fMat->SetRhoF(1.);
    fMat->SetHydrostatic([gw, D](const TPZVec<REAL> &x) { return gw * (D - x[1]); });
    fMPhys->InsertMaterialObject(fMat);
    TPZFNMatrix<9, STATE> v1(3, 3, 0.);
    TPZManVector<STATE, 4> zero(4, 0.), rollerX = {1., 0., 0., 0.};
    fMPhys->InsertMaterialObject(fMat->CreateBC(fMat, EBottom, TMatUP::EDirichletU, v1, zero));
    for (int id : {ELeft, ERight}) fMPhys->InsertMaterialObject(fMat->CreateBC(fMat, id, TMatUP::EDirectionalNullU, v1, rollerX));
    // reservatório: p = p_w e t = -p_w n (n exterior) na crista, na face e no pé
    const REAL beta = geo.betaDeg * M_PI / 180.;
    const REAL normals[3][3] = {{REAL(ECrest), 0., 1.}, {REAL(EFace), std::sin(beta), std::cos(beta)}, {REAL(EToe), 0., 1.}};
    for (const auto &nrm : normals) {
        auto *bc = fMat->CreateBC(fMat, int(nrm[0]), TMatUP::ETractionDirichletP, v1, zero);
        const REAL nx = nrm[1], ny = nrm[2];
        bc->SetForcingFunctionBC([this, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) {
            const REAL pw = fGravity * fSoil.gammaW * std::max(fZw - x[1], REAL(0.));
            v[0] = -pw * nx;
            v[1] = -pw * ny;
            v[2] = 0.;
            v[3] = pw;
        });
        fMPhys->InsertMaterialObject(bc);
    }
    TPZManVector<int, 2> active = {1, 1};
    TPZManVector<TPZCompMesh *, 2> meshvec = {fCMeshU, fCMeshP};
    fMPhys->BuildMultiphysicsSpace(active, meshvec);
    // o material u-p usa os índices de memória dos elementos de u (datavec[0].intGlobPtIndex)
    fMat->GetMemory() = fMatU->GetMemory();
    SetGravityFactor(1.);

    // as regras de integração do multifísico e da malha de u devem coincidir (índices de memória)
    for (TPZCompEl *cel : fMPhys->ElementVec()) {
        auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mf || !cel->Reference() || cel->Reference()->Dimension() != 2) continue;
        TPZCompEl *celu = mf->Element(0);
        if (!celu || cel->GetIntegrationRule().NPoints() != celu->GetIntegrationRule().NPoints())
            throw std::runtime_error("TCoupledDrawdown: regras de integração de u e do multifísico diferentes");
    }

    // ---- pontos de integração (índices de memória) e elementos de p por elemento geométrico
    const int64_t nmem = fMatU->GetMemory()->NElements();
    fPoints.assign(nmem, TPoint{-1, TPZManVector<REAL, 3>(2, 0.), TPZManVector<REAL, 3>(3, 0.)});
    for (TPZCompEl *cel : fCMeshU->ElementVec()) {
        auto *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel || !cel->Reference() || cel->Reference()->Dimension() != 2 || cel->Material() != fMatU) continue;
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        TPZMaterialDataT<STATE> data;
        intel->InitMaterialData(data);
        TPZManVector<REAL, 3> qsi(2, 0.);
        REAL w;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            rule.Point(ip, qsi, w);
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData(data, qsi);
            const int64_t idx = data.intGlobPtIndex;
            fPoints[idx].gel = cel->Reference()->Index();
            fPoints[idx].qsi = qsi;
            fPoints[idx].x = data.x;
        }
    }
    fCelP.assign(gmesh->NElements(), nullptr);
    for (TPZCompEl *cel : fCMeshP->ElementVec())
        if (cel && cel->Reference() && cel->Reference()->Dimension() == 2) fCelP[cel->Reference()->Index()] = cel;
    fC.assign(fPoints.size(), soil.c);
    fPhi.assign(fPoints.size(), soil.phiDeg * M_PI / 180.);

    // ---- análise: renumeração de banda, skyline não simétrica, LU
    fAn = std::make_unique<TPZLinearAnalysis>(fMPhys, true);
    TPZSkylineNSymStructMatrix<STATE> str(fMPhys);
    str.SetNumThreads(0);
    fAn->SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    fAn->SetSolver(step);
}

template <class T>
TCoupledDrawdown<T>::~TCoupledDrawdown() {
    delete fPost;
    fAn.reset();
    delete fMPhys;
    delete fCMeshP;
    delete fCMeshU;
}

template <class T>
REAL TCoupledDrawdown<T>::Cv() const {
    const REAL E = fSoil.E, nu = fSoil.nu;
    const REAL M = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu));
    return fPar.kv * M / fSoil.gammaW;
}

template <class T>
REAL TCoupledDrawdown<T>::TimeScale() const {
    return fGeo.H * fGeo.H / Cv();
}

template <class T>
REAL TCoupledDrawdown<T>::WaterLevelAt(REAL Tad) const {
    const REAL s = fPar.Td > 0. ? std::min(Tad / fPar.Td, REAL(1.)) : 1.;
    return fGeo.D() - fPar.hw * std::max(s, REAL(0.));
}

template <class T>
void TCoupledDrawdown<T>::SetStrength(const std::vector<REAL> &c, const std::vector<REAL> &phi) {
    fC = c;
    fPhi = phi;
}

template <class T>
void TCoupledDrawdown<T>::SetInitialStress(const std::vector<TPZTensor<REAL>> &sigma0) {
    fSigma0 = sigma0;
}

template <class T>
void TCoupledDrawdown<T>::SetGravityFactor(REAL f) {
    fGravity = f;
    fMat->SetBodyForce(0., -f * fSoil.gamma, 0.);  // peso específico saturado (tensões totais)
    fMat->SetGravity(0., -f * fSoil.gammaW, 0.);   // ρ_f g com ρ_f = 1: q = -(k/γw)(∇p + γw e_y)
}

template <class T>
void TCoupledDrawdown<T>::InitializeMemory() {
    auto &memory = *fMatU->GetMemory();
    T model;
    SetupSoilModel(model, fSoil);
    const REAL gw = fSoil.gammaW, D = fGeo.D();
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        if (fPoints[i].gel < 0) continue;
        TPZElastoPlasticMem &m = memory[i];
        m.m_elastoplastic_state.CleanUp();
        TPZVec<REAL> &mp = m.m_elastoplastic_state.fmatprop;
        if constexpr (IsCamClay<T>()) {
            if (fSigma0.size() != fPoints.size())
                throw std::runtime_error("TCoupledDrawdown: Cam-Clay exige a tensão inicial (SetInitialStress)");
            mp.Resize(9);
            mp.Fill(0.);
            mp[0] = fC[i];
            mp[1] = fPhi[i];
            TPZTensor<REAL> s = fSigma0[i];
            CamClayInitialState(mp, s, fSoil);
            m.m_sigma = s;
            m.fPorePressure = gw * (D - fPoints[i].x[1]);
        } else {
            mp.Resize(3);
            mp[0] = fC[i];
            mp[1] = fPhi[i];
            mp[2] = fPhi[i];
            m.m_sigma.Zero();
            m.fPorePressure = 0.;
        }
        m.m_elastoplastic_state.fmatpropinit = mp;
        m.fdPorePressure.Resize(0);
        m.m_ER = model.GetElasticResponse();
        m.m_u.Resize(2);
        m.m_u.Fill(0.);
        m.m_plastic_steps = 0;
    }
}

template <class T>
bool TCoupledDrawdown<T>::Step(REAL dt, bool flow, int &iterations) {
    fMat->SetTimeStep(dt);
    fMat->SetFlow(flow);
    fMat->SetAssembleMode(TMatUP::EExternalForces);
    fAn->AssembleResidual();
    fMat->SetAssembleMode(TMatUP::EFull);
    const REAL ref = std::max(REAL(Norm(fAn->Rhs())), REAL(1.));
    const TPZFMatrix<STATE> xn = fAn->Solution();
    TPZFMatrix<STATE> x = xn;
    auto load = [&](const TPZFMatrix<STATE> &v) {
        fAn->LoadSolution(v);
        fMPhys->LoadSolutionFromMultiPhysics();
    };
    // ||R|| nas equações sem penalidade de Dirichlet
    auto freeNorm = [&](const TPZFMatrix<STATE> &rhs) {
        REAL s = 0.;
        for (int64_t i = 0; i < (int64_t)fFree.size(); i++)
            if (fFree[i]) s += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
        return std::sqrt(s) / ref;
    };
    iterations = 0;
    REAL nr = 0., nrMin = 1.e300;
    int sinceMin = 0;
    try {
        for (int it = 1; it <= fPar.maxIter; it++) {
            iterations = it;
            fAn->Assemble();
            const int64_t neq = fMPhys->NEquations();
            if (fFree.empty()) {  // equações sem penalidade de Dirichlet
                auto K = fAn->MatrixSolver<STATE>().Matrix();
                fFree.resize(neq);
                for (int64_t i = 0; i < neq; i++) fFree[i] = std::fabs(K->GetVal(i, i)) < 1.e-6 * TMatUP::BigNumber();
            }
            nr = freeNorm(fAn->Rhs());
            if (fPar.verbose > 1) std::cout << "      it " << it << " |R|/|F| = " << nr << "\n";
            if (!std::isfinite(nr)) break;
            if (it > 1 && nr < fPar.tol) {
                fMat->SetUpdateMem(true);  // aceita o passo: memória dos pontos de integração
                fAn->AssembleResidual();
                fMat->SetUpdateMem(false);
                return true;
            }
            if (it > 1) {  // sem redução do resíduo em 6 iterações: não converge
                if (nr < nrMin) {
                    nrMin = nr;
                    sinceMin = 0;
                } else if (++sinceMin >= 6) {
                    break;
                }
            }
            fAn->Solve();
            const TPZFMatrix<STATE> dx = fAn->Solution();
            if (it == 1) {  // primeira iteração completa (valores de Dirichlet novos)
                x += dx;
                load(x);
                continue;
            }
            // busca linear por bissecção no resíduo das equações livres (pontos que alternam entre carga
            // plástica e descarga elástica fazem o Newton oscilar, sobretudo no Cam-Clay normalmente adensado)
            REAL alpha = 1., best = 1.e300, alphaBest = 1.;
            for (int ls = 0; ls < 4; ls++) {
                TPZFMatrix<STATE> xt(x);
                xt.ZAXPY(alpha, dx);
                load(xt);
                fAn->AssembleResidual();
                const REAL nt = freeNorm(fAn->Rhs());
                if (nt < best) {
                    best = nt;
                    alphaBest = alpha;
                }
                if (nt < nr) break;
                alpha *= 0.5;
            }
            x.ZAXPY(alphaBest, dx);
            load(x);
        }
    } catch (...) {
    }
    fMat->SetUpdateMem(false);
    load(xn);
    if (fPar.verbose) std::cout << "    passo sem convergência (dt = " << dt << " s, |R|/|F| = " << nr << ")\n";
    return false;
}

template <class T>
bool TCoupledDrawdown<T>::Initialize() {
    InitializeMemory();
    fFree.clear();
    fZw = fGeo.D();
    fTime = 0.;
    fAn->Solution().Zero();
    fAn->LoadSolution();
    fMPhys->LoadSolutionFromMultiPhysics();
    const REAL dtDrained = 1.e6 * TimeScale();  // passo drenado: p hidrostática
    int its = 0;
    const int nG = IsCamClay<T>() ? 1 : std::max(fPar.nGravity, 1);
    for (int k = 1; k <= nG; k++) {
        SetGravityFactor(REAL(k) / nG);
        if (!Step(dtDrained, true, its)) return false;
        if (fPar.verbose) std::cout << "    peso próprio " << k << "/" << nG << ": " << its << " it.\n";
    }
    // deslocamentos de referência nos pontos de integração
    const auto &memory = *fMatU->GetMemory();
    fU0.assign(fPoints.size(), {0., 0.});
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        const TPZVec<REAL> &u = memory[i].m_u;
        if (fPoints[i].gel >= 0 && u.size() >= 2) fU0[i] = {u[0], u[1]};
    }
    fDtAfter = fPar.Td * TimeScale() / std::max(fPar.nDrawdown, 1);
    return true;
}

template <class T>
bool TCoupledDrawdown<T>::AdvanceTo(REAL Ttarget, REAL dtMax) {
    const REAL ts = TimeScale(), td = fPar.Td * ts, tEnd = Ttarget * ts;
    REAL cut = 1.;
    int ncuts = 0;
    while (fTime < tEnd * (1. - 1.e-12)) {
        const REAL dtNom = (fTime < td * (1. - 1.e-12)) ? td / std::max(fPar.nDrawdown, 1) : fDtAfter;
        REAL dt = std::min({dtNom * cut, tEnd - fTime, dtMax * ts});
        if (fTime < td && fTime + dt > td) dt = td - fTime;  // não atravessa o fim do rebaixamento
        const REAL zwOld = fZw;
        fZw = WaterLevelAt((fTime + dt) / ts);
        int its = 0;
        if (Step(dt, true, its)) {
            fTime += dt;
            if (fPar.verbose)
                std::cout << "    T = " << fTime / ts << " z_w = " << fZw << " (" << its << " it.)\n";
            if (fTime >= td * (1. - 1.e-12)) fDtAfter = std::max(fDtAfter, dt) * fPar.growth;
            if (cut < 1.) cut = std::min(1., 2. * cut);
            ncuts = 0;
        } else {
            fZw = zwOld;
            cut = 0.5 * dt / dtNom;  // metade do passo tentado
            if (++ncuts > fPar.maxCuts) return false;
        }
    }
    return true;
}

template <class T>
void TCoupledDrawdown<T>::ExcessPorePressure(int64_t gel, const TPZVec<REAL> &qsi, REAL &u,
                                             TPZManVector<REAL, 2> &gradu) const {
    auto *intel = dynamic_cast<TPZInterpolationSpace *>(fCelP[gel]);
    if (!intel) throw std::runtime_error("TCoupledDrawdown: elemento sem pressão");
    TPZMaterialDataT<STATE> data;
    TPZManVector<REAL, 3> q(qsi), x(3);
    intel->ComputeSolution(q, data, false);
    fGMesh->Element(gel)->X(q, x);
    const REAL gw = fSoil.gammaW;
    u = data.sol[0][0] - gw * (fGeo.D() - x[1]);
    gradu.Resize(2);
    for (int d = 0; d < 2; d++) gradu[d] = data.axes(0, d) * data.dsol[0](0, 0) + data.axes(1, d) * data.dsol[0](1, 0);
    gradu[1] += gw;
}

template <class T>
REAL TCoupledDrawdown<T>::PorePressureAt(const TPZVec<REAL> &x) const {
    for (int64_t ig = 0; ig < (int64_t)fCelP.size(); ig++) {
        if (!fCelP[ig]) continue;
        TPZGeoEl *gel = fGMesh->Element(ig);
        REAL mn[2] = {1.e300, 1.e300}, mx[2] = {-1.e300, -1.e300};
        for (int i = 0; i < gel->NCornerNodes(); i++)
            for (int d = 0; d < 2; d++) {
                mn[d] = std::min(mn[d], gel->NodePtr(i)->Coord(d));
                mx[d] = std::max(mx[d], gel->NodePtr(i)->Coord(d));
            }
        if (x[0] < mn[0] - 1.e-9 || x[0] > mx[0] + 1.e-9 || x[1] < mn[1] - 1.e-9 || x[1] > mx[1] + 1.e-9) continue;
        TPZManVector<REAL, 3> qsi(2, 0.), xx(3, 0.);
        xx[0] = x[0];
        xx[1] = x[1];
        if (!gel->ComputeXInverse(xx, qsi, 1.e-10)) continue;
        REAL u;
        TPZManVector<REAL, 2> g;
        ExcessPorePressure(ig, qsi, u, g);
        return u + fSoil.gammaW * (fGeo.D() - x[1]);
    }
    return std::nan("");
}

template <class T>
REAL TCoupledDrawdown<T>::MaxDisplacementIncrement() const {
    const auto &memory = *fMatU->GetMemory();
    REAL m = 0.;
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++) {
        if (fPoints[i].gel < 0 || i >= (int64_t)fU0.size()) continue;
        const TPZVec<REAL> &u = memory[i].m_u;
        if (u.size() < 2) continue;
        m = std::max(m, std::hypot(u[0] - fU0[i][0], u[1] - fU0[i][1]));
    }
    return m;
}

template <class T>
int64_t TCoupledDrawdown<T>::NPlasticPoints() const {
    const auto &memory = *fMatU->GetMemory();
    int64_t n = 0;
    for (int64_t i = 0; i < (int64_t)fPoints.size(); i++)
        if (fPoints[i].gel >= 0 && memory[i].m_elastoplastic_state.m_m_type != 0) n++;
    return n;
}

template <class T>
void TCoupledDrawdown<T>::DefineVTK(const std::string &base) {
    TPZStack<std::string> scal, vec;
    scal.Push("Pressure");
    scal.Push("ExcessPressure");
    vec.Push("Displacement");
    vec.Push("Flux");
    fAn->DefineGraphMesh(2, scal, vec, base + "_up.vtk");
    delete fPost;
    fPost = new TPZPostProcAnalysis();
    fPost->SetCompMesh(fCMeshU);
    TPZFStructMatrix<STATE> structmatrix(fPost->Mesh());
    structmatrix.SetNumThreads(0);
    fPost->SetStructuralMatrix(structmatrix);
    TPZStack<std::string> pscal, pvec, vars;
    for (const char *nm : {"PlasticStrainNorm", "StressXX", "StressYY", "SqrtStressJ2"}) pscal.Push(nm);
    for (auto &s : pscal) vars.Push(s);
    TPZVec<int> matids(1, ESoil);
    fPost->SetPostProcessVariables(matids, vars);
    fPost->DefineGraphMesh(2, pscal, pvec, base + "_tensoes.vtk");
}

template <class T>
void TCoupledDrawdown<T>::WriteVTK(int step) {
    fAn->SetStep(step);
    fAn->PostProcess(1);
    if (fPost) {
        fPost->TransferSolution();
        fPost->SetStep(step);
        fPost->PostProcess(1);
    }
}

template class TCoupledDrawdown<TMohrCoulomb>;
template class TCoupledDrawdown<TPZModifiedCamClay>;
