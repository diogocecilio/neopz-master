// PoroCamClayNativo.cpp — ver PoroCamClayNativo.h

#include "PoroCamClayNativo.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>

#include "TPZBndCondT.h"
#include "TPZNullMaterial.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "malhas.h"

// =====================================================================================================
// TPZMatCamClayPostProc
// =====================================================================================================
int TPZMatCamClayPostProc::IntegrationRuleOrder(const int elPMaxOrder) const {
    if (fIntegrationOrder > 0) return fIntegrationOrder;
    return TBase::IntegrationRuleOrder(elPMaxOrder);
}

int TPZMatCamClayPostProc::VariableIndex(const std::string &name) const {
    if (name == "PorePressure") return EPorePressure;
    if (name == "ExcessPorePressure") return EExcessPorePressure;
    if (name == "MeanEffectiveStress") return EMeanEffectiveStress;
    if (name == "DeviatoricStress") return EDeviatoricStress;
    if (name == "PreconsolidationPressure") return EPreconsolidation;
    if (name == "PlasticPoint") return EPlasticPoint;
    if (name == "EffectiveStress") return EEffectiveStress;
    if (name == "TotalStress") return ETotalStress;
    return TBase::VariableIndex(name);
}

int TPZMatCamClayPostProc::NSolutionVariables(int var) const {
    if (var == EEffectiveStress || var == ETotalStress) return 9;
    if (var >= EPorePressure && var <= EPlasticPoint) return 1;
    if (var == TBase::EVolHardening) return 1;
    return TBase::NSolutionVariables(var);
}

void TPZMatCamClayPostProc::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) {
    if (var == TBase::EVolHardening) {  // α = -ε_v^p (não tratado pelo TPZMatElastoPlastic::Solution)
        Solout.Resize(1);
        Solout[0] = data.intGlobPtIndex < 0 ? 0. : this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_hardening;
        return;
    }
    if (var < EPorePressure) {
        TBase::Solution(data, var, Solout);
        return;
    }
    Solout.Resize(NSolutionVariables(var));
    Solout.Fill(0.);
    const int64_t ip = data.intGlobPtIndex;
    if (ip < 0) return;
    const TPZElastoPlasticMem &mem = this->MemItem(ip);
    const TPZTensor<REAL> &s = mem.m_sigma;
    auto tensor = [&Solout](const TPZTensor<REAL> &t) {
        const REAL v[9] = {t[_XX_], t[_XY_], t[_XZ_], t[_XY_], t[_YY_], t[_YZ_], t[_XZ_], t[_YZ_], t[_ZZ_]};
        for (int i = 0; i < 9; i++) Solout[i] = v[i];
    };
    REAL p, q;
    switch (var) {
        case EPorePressure:
            Solout[0] = mem.fPorePressure;
            break;
        case EExcessPorePressure:
            Solout[0] = mem.fPorePressure - (fHydro ? fHydro(data.x) : 0.);
            break;
        case EMeanEffectiveStress:
            TPZModifiedCamClay::Invariants(s, p, q);
            Solout[0] = -p;
            break;
        case EDeviatoricStress:
            TPZModifiedCamClay::Invariants(s, p, q);
            Solout[0] = q;
            break;
        case EPreconsolidation: {
            TPZModifiedCamClay m(this->GetPlasticModel());
            if (fModelUpdate) fModelUpdate(data.x, m);
            Solout[0] = m.Pc(mem.m_elastoplastic_state.m_hardening);
        } break;
        case EPlasticPoint:
            Solout[0] = mem.m_elastoplastic_state.m_m_type;
            break;
        case EEffectiveStress:
            tensor(s);
            break;
        case ETotalStress: {
            TPZTensor<REAL> st(s);
            for (int c : {_XX_, _YY_, _ZZ_}) st[c] -= fAlpha * mem.fPorePressure;
            tensor(st);
        } break;
        default:
            break;
    }
}

// =====================================================================================================
// TPoroCamClayNativo
// =====================================================================================================
TPoroCamClayNativo::TPoroCamClayNativo(TPZGeoMesh *gmesh, const TParams &par) : fPar(par), fGMesh(gmesh) {
    // ---- malha de u (H1 vetorial com memória nos pontos de integração)
    fCMeshU = new TPZCompMesh(gmesh);
    fCMeshU->SetDimModel(3);
    fCMeshU->SetDefaultOrder(par.orderU);
    fCMeshU->SetAllCreateFunctionsContinuousWithMem();
    fMatU = new TPZMatCamClayPostProc(kMatVolume);
    TPZModifiedCamClay model(par.model);
    fMatU->SetPlasticityModel(model);
    fMatU->SetModelUpdate(par.modelUpdate);
    fMatU->SetHydrostatic(par.hydro);
    fMatU->SetAlpha(par.alpha);
    fMatU->SetIntegrationOrder(par.integrationOrder);
    fCMeshU->InsertMaterialObject(fMatU);
    {
        TPZFNMatrix<9, STATE> v1(3, 3, 0.);
        TPZManVector<STATE, 3> v2(3, 0.);
        for (auto &b : par.bcs) fCMeshU->InsertMaterialObject(fMatU->CreateBC(fMatU, b.id, 0, v1, v2));
    }
    fCMeshU->AutoBuild();

    // ---- malha de p (H1 escalar linear)
    fCMeshP = new TPZCompMesh(gmesh);
    fCMeshP->SetDimModel(3);
    fCMeshP->SetDefaultOrder(1);
    fCMeshP->SetAllCreateFunctionsContinuous();
    auto *nullmat = new TPZNullMaterial<STATE>(kMatVolume, 3, 1);
    fCMeshP->InsertMaterialObject(nullmat);
    {
        TPZFNMatrix<1, STATE> v1(1, 1, 0.);
        TPZManVector<STATE, 1> v2(1, 0.);
        for (auto &b : par.bcs) fCMeshP->InsertMaterialObject(nullmat->CreateBC(nullmat, b.id, 0, v1, v2));
    }
    fCMeshP->AutoBuild();

    // ---- malha multifísica u-p
    fMPhys = new TPZMultiphysicsCompMesh(gmesh);
    fMPhys->SetDimModel(3);
    fMPhys->SetAllCreateFunctionsMultiphysicElem();
    fMat = new TMatUP(kMatVolume);
    fMat->SetPlasticModel(par.model);
    fMat->SetModelUpdate(par.modelUpdate);
    fMat->SetHydrostatic(par.hydro);
    fMat->SetAlpha(par.alpha);
    fMat->SetSe(par.Se);
    fMat->SetPermeability(par.k);
    fMat->SetViscosity(par.mu);
    fMat->SetRhoF(par.rhof);
    fMat->SetGravity(par.g[0], par.g[1], par.g[2]);
    fMat->SetBodyForce(par.body[0], par.body[1], par.body[2]);
    fMat->SetIntegrationOrder(par.integrationOrder);
    fMPhys->InsertMaterialObject(fMat);
    for (auto &b : par.bcs) {
        TPZFNMatrix<9, STATE> v1(3, 3, 0.);
        for (int i = 0; i < 3; i++) v1(i, i) = b.v1diag[i];
        TPZManVector<STATE, 4> v2 = {b.v2[0], b.v2[1], b.v2[2], b.v2[3]};
        fMPhys->InsertMaterialObject(fMat->CreateBC(fMat, b.id, b.type, v1, v2));
    }
    TPZManVector<int, 2> active = {1, 1};
    TPZManVector<TPZCompMesh *, 2> meshvec = {fCMeshU, fCMeshP};
    fMPhys->BuildMultiphysicsSpace(active, meshvec);
    // o material u-p usa os índices de memória dos elementos de u (datavec[0].intGlobPtIndex):
    // a memória é a mesma do material da malha de u, lida no pós-processamento
    fMat->GetMemory() = fMatU->GetMemory();

    BuildVertexMaps();
    InitializeMemory();

    // ---- poropressão inicial nos vértices
    if (par.p0) {
        for (auto &kv : fVertexP) {
            TPZCompEl *cel = kv.second.first;
            const TPZConnect &con = cel->Connect(kv.second.second);
            TPZManVector<REAL, 3> x(3);
            fGMesh->NodeVec()[kv.first].GetCoordinates(x);
            TPZFMatrix<STATE> &solP = fCMeshP->Solution();
            solP(fCMeshP->Block().Position(con.SequenceNumber()), 0) = par.p0(x);
        }
        fMPhys->LoadSolutionFromMeshes();
    }

    // ---- análise: renumeração de banda, skyline não simétrica, LU
    fAn = new TPZLinearAnalysis(fMPhys, true);
    TPZSkylineNSymStructMatrix<STATE> str(fMPhys);
    fAn->SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    fAn->SetSolver(step);

    // equações de u (conectores da malha de u vêm primeiro na malha multifísica)
    const int64_t neq = fMPhys->NEquations();
    fIsU.assign(neq, false);
    for (int64_t ic = 0; ic < fCMeshU->NConnects(); ic++) {
        const TPZConnect &c = fMPhys->ConnectVec()[ic];
        if (c.SequenceNumber() < 0) continue;
        const int64_t pos = fMPhys->Block().Position(c.SequenceNumber());
        for (int k = 0; k < fMPhys->Block().Size(c.SequenceNumber()); k++) fIsU[pos + k] = true;
    }
}

TPoroCamClayNativo::~TPoroCamClayNativo() {
    delete fPost;
    delete fAn;
    delete fMPhys;
    delete fCMeshP;
    delete fCMeshU;
}

void TPoroCamClayNativo::BuildVertexMaps() {
    for (TPZCompEl *cel : fCMeshU->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 3) continue;
        TPZGeoEl *gel = cel->Reference();
        fGelToU[gel] = cel;
        for (int i = 0; i < gel->NCornerNodes(); i++) fVertexU[gel->NodeIndex(i)] = {cel, i};
    }
    for (TPZCompEl *cel : fCMeshP->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != 3) continue;
        TPZGeoEl *gel = cel->Reference();
        for (int i = 0; i < gel->NCornerNodes(); i++) fVertexP[gel->NodeIndex(i)] = {cel, i};
    }
}

void TPoroCamClayNativo::InitializeMemory() {
    auto &mem = *fMatU->GetMemory();
    for (auto &kv : fGelToU) {
        TPZGeoEl *gel = kv.first;
        TPZCompEl *cel = kv.second;
        TPZManVector<int64_t> idx;
        cel->GetMemoryIndices(idx);
        const TPZIntPoints &rule = cel->GetIntegrationRule();
        if (rule.NPoints() != idx.size()) DebugStop();
        TPZManVector<REAL, 3> qsi(3), x(3);
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            REAL w;
            rule.Point(ip, qsi, w);
            gel->X(qsi, x);
            TPZModifiedCamClay m(fPar.model);
            if (fPar.modelUpdate) fPar.modelUpdate(x, m);
            TPZElastoPlasticMem &it = mem[idx[ip]];
            it.m_elastoplastic_state.CleanUp();
            it.m_sigma = m.InitialStress();
            it.fPorePressure = fPar.p0 ? fPar.p0(x) : 0.;
            it.m_ER = m.GetElasticResponse();
            it.m_u.Resize(3);
            it.m_u.Fill(0.);
        }
    }
}

int64_t TPoroCamClayNativo::MPhysEquation(int64_t uconnect, int comp) const {
    const TPZConnect &c = fMPhys->ConnectVec()[uconnect];
    return fMPhys->Block().Position(c.SequenceNumber()) + comp;
}

void TPoroCamClayNativo::SetBCVal2(int id, const std::array<REAL, 4> &v2) {
    auto *bc = dynamic_cast<TPZBndCondT<STATE> *>(fMPhys->FindMaterial(id));
    if (!bc) throw std::runtime_error("condição de contorno não encontrada");
    TPZManVector<STATE, 4> v = {v2[0], v2[1], v2[2], v2[3]};
    bc->SetVal2(v);
}

REAL TPoroCamClayNativo::ExternalForceNorm() {
    fMat->SetAssembleMode(TMatUP::EExternalForces);
    fAn->AssembleResidual();
    fMat->SetAssembleMode(TMatUP::EFull);
    TPZFMatrix<STATE> rhs = fAn->Rhs();
    return Norm(rhs);
}

REAL TPoroCamClayNativo::FreeResidualU() {
    fAn->Assemble();
    const int64_t neq = fMPhys->NEquations();
    if (fFree.empty()) {
        auto K = fAn->MatrixSolver<STATE>().Matrix();
        fFree.resize(neq);
        for (int64_t i = 0; i < neq; i++) fFree[i] = std::fabs(K->GetVal(i, i)) < 1.e-6 * TMatUP::BigNumber();
    }
    const TPZFMatrix<STATE> rhs = fAn->Rhs();
    REAL m = 0.;
    for (int64_t i = 0; i < neq; i++)
        if (fFree[i] && fIsU[i]) m = std::max(m, std::fabs(rhs.GetVal(i, 0)));
    return m;
}

int TPoroCamClayNativo::Step(REAL dt, bool flow, bool dirichletChanged, const TPZFMatrix<STATE> *guess, REAL tol,
                             int maxit) {
    fMat->SetTimeStep(dt);
    fMat->SetFlow(flow);
    const REAL ref = std::max(ExternalForceNorm(), REAL(1.));
    const TPZFMatrix<STATE> xn = fAn->Solution();
    TPZFMatrix<STATE> x = guess ? *guess : xn;
    auto load = [&](const TPZFMatrix<STATE> &v) {
        fAn->LoadSolution(v);
        fMPhys->LoadSolutionFromMultiPhysics();
    };
    load(x);
    REAL nr = 0.;
    try {
        for (int it = 1; it <= maxit; it++) {
            fAn->Assemble();
            const int64_t neq = fMPhys->NEquations();
            if (fFree.empty()) {  // equações sem penalidade de Dirichlet
                auto K = fAn->MatrixSolver<STATE>().Matrix();
                fFree.resize(neq);
                for (int64_t i = 0; i < neq; i++)
                    fFree[i] = std::fabs(K->GetVal(i, i)) < 1.e-6 * TMatUP::BigNumber();
            }
            const TPZFMatrix<STATE> rhs = fAn->Rhs();
            REAL s = 0.;
            for (int64_t i = 0; i < neq; i++)
                if (fFree[i]) s += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
            nr = std::sqrt(s) / ref;
            if (nr < tol && (!dirichletChanged || it > 1)) {
                fMat->SetUpdateMem(true);  // aceita o passo: atualiza a memória dos pontos de integração
                fAn->AssembleResidual();
                fMat->SetUpdateMem(false);
                return it;
            }
            fAn->Solve();
            const TPZFMatrix<STATE> dx = fAn->Solution();
            x += dx;
            load(x);
        }
    } catch (...) {
        load(xn);
        throw;
    }
    load(xn);
    std::stringstream sout;
    sout << "TPoroCamClayNativo::Step: Newton não convergiu (|R|/|F| = " << nr << ")";
    throw std::runtime_error(sout.str());
}

void TPoroCamClayNativo::SetVertexU(TPZFMatrix<STATE> &vec,
                                    const std::function<std::array<REAL, 3>(const TPZVec<REAL> &)> &f) const {
    for (auto &kv : fVertexU) {
        TPZManVector<REAL, 3> x(3);
        fGMesh->NodeVec()[kv.first].GetCoordinates(x);
        const auto v = f(x);
        const int64_t ic = kv.second.first->ConnectIndex(kv.second.second);
        for (int c = 0; c < 3; c++) vec(MPhysEquation(ic, c), 0) = v[c];
    }
}

REAL TPoroCamClayNativo::NodalU(int64_t node, int comp) const {
    auto it = fVertexU.find(node);
    if (it == fVertexU.end()) throw std::runtime_error("nó sem grau de liberdade de vértice");
    const TPZConnect &con = it->second.first->Connect(it->second.second);
    const TPZFMatrix<STATE> &sol = fCMeshU->Solution();
    return sol.GetVal(fCMeshU->Block().Position(con.SequenceNumber()) + comp, 0);
}

REAL TPoroCamClayNativo::NodalP(int64_t node) const {
    auto it = fVertexP.find(node);
    if (it == fVertexP.end()) throw std::runtime_error("nó sem grau de liberdade de pressão");
    const TPZConnect &con = it->second.first->Connect(it->second.second);
    const TPZFMatrix<STATE> &sol = fCMeshP->Solution();
    return sol.GetVal(fCMeshP->Block().Position(con.SequenceNumber()), 0);
}

REAL TPoroCamClayNativo::ZonePressure(TPZGeoEl *gel) const {
    REAL s = 0.;
    for (int i = 0; i < 8; i++) s += NodalP(gel->NodeIndex(i));
    return s / 8.;
}

REAL TPoroCamClayNativo::MaxAbsPressure() const {
    REAL m = 0.;
    for (auto &kv : fVertexP) m = std::max(m, std::fabs(NodalP(kv.first)));
    return m;
}

TPZGeoEl *TPoroCamClayNativo::FindVolumeElement(const TPZVec<REAL> &x, TPZVec<REAL> &qsi) const {
    for (auto &kv : fGelToU) {
        TPZGeoEl *gel = kv.first;
        bool fora = false;
        for (int d = 0; d < 3 && !fora; d++) {
            REAL mn = 1.e300, mx = -1.e300;
            for (int i = 0; i < gel->NNodes(); i++) {
                TPZManVector<REAL, 3> c(3);
                gel->NodePtr(i)->GetCoordinates(c);
                mn = std::min(mn, c[d]);
                mx = std::max(mx, c[d]);
            }
            if (x[d] < mn - 1.e-12 || x[d] > mx + 1.e-12) fora = true;
        }
        if (fora) continue;
        qsi.Resize(3);
        qsi.Fill(0.);
        TPZManVector<REAL, 3> xx(x);
        if (gel->ComputeXInverse(xx, qsi, 1.e-12)) return gel;
    }
    return nullptr;
}

TPZTensor<REAL> TPoroCamClayNativo::StressAt(TPZGeoEl *gel, const TPZVec<REAL> &xi) const {
    TPZCompEl *cel = fGelToU.at(gel);
    TPZManVector<int64_t> idx;
    cel->GetMemoryIndices(idx);
    const TPZIntPoints &rule = cel->GetIntegrationRule();
    std::vector<std::array<REAL, 3>> pts(rule.NPoints());
    std::vector<REAL> g;
    for (int ip = 0; ip < rule.NPoints(); ip++) {
        TPZManVector<REAL, 3> q(3);
        REAL w;
        rule.Point(ip, q, w);
        pts[ip] = {q[0], q[1], q[2]};
        if (std::none_of(g.begin(), g.end(), [&](REAL v) { return std::fabs(v - q[0]) < 1.e-12; })) g.push_back(q[0]);
    }
    std::sort(g.begin(), g.end());
    const int n = int(g.size());
    auto lag = [&](REAL t, int l) {
        REAL v = 1.;
        for (int m = 0; m < n; m++)
            if (m != l) v *= (t - g[m]) / (g[l] - g[m]);
        return v;
    };
    auto near = [&](REAL t) {
        int b = 0;
        for (int m = 1; m < n; m++)
            if (std::fabs(g[m] - t) < std::fabs(g[b] - t)) b = m;
        return b;
    };
    TPZTensor<REAL> s;
    const auto &mem = *fMatU->GetMemory();
    for (int ip = 0; ip < rule.NPoints(); ip++) {
        const REAL w = lag(xi[0], near(pts[ip][0])) * lag(xi[1], near(pts[ip][1])) * lag(xi[2], near(pts[ip][2]));
        for (int i = 0; i < 6; i++) s[i] += w * mem[idx[ip]].m_sigma[i];
    }
    return s;
}

REAL TPoroCamClayNativo::SumReactions(const std::vector<int64_t> &nodes, int comp) {
    fMat->SetAssembleMode(TMatUP::ENoPenalty);
    fAn->AssembleResidual();
    fMat->SetAssembleMode(TMatUP::EFull);
    const TPZFMatrix<STATE> rhs = fAn->Rhs();  // = F_ext - F_int + Q p = -R
    REAL s = 0.;
    for (int64_t nd : nodes) {
        auto it = fVertexU.find(nd);
        if (it == fVertexU.end()) continue;
        s -= rhs.GetVal(MPhysEquation(it->second.first->ConnectIndex(it->second.second), comp), 0);
    }
    return s;
}

void TPoroCamClayNativo::DefineVTK(const std::string &base) {
    // campos contínuos (u, p) na malha multifísica
    TPZStack<std::string> scal, vec;
    scal.Push("Pressure");
    scal.Push("ExcessPressure");
    vec.Push("Displacement");
    vec.Push("Flux");
    fAn->DefineGraphMesh(3, scal, vec, base + "_up.vtk");
    // variáveis guardadas nos pontos de integração: TPZPostProcAnalysis sobre a malha de u
    fPost = new TPZPostProcAnalysis();
    fPost->SetCompMesh(fCMeshU);
    TPZFStructMatrix<STATE> structmatrix(fPost->Mesh());
    fPost->SetStructuralMatrix(structmatrix);
    TPZStack<std::string> pscal, pvec, ptens, vars;
    for (const char *nm : {"PorePressure", "ExcessPorePressure", "MeanEffectiveStress", "DeviatoricStress",
                           "VolHardening", "PreconsolidationPressure", "PlasticPoint"})
        pscal.Push(nm);
    pvec.Push("Displacement");
    ptens.Push("EffectiveStress");
    ptens.Push("TotalStress");
    for (auto &s : pscal) vars.Push(s);
    for (auto &s : pvec) vars.Push(s);
    for (auto &s : ptens) vars.Push(s);
    TPZVec<int> matids(1, kMatVolume);
    fPost->SetPostProcessVariables(matids, vars);
    fPost->DefineGraphMesh(3, pscal, pvec, ptens, base + "_tensoes.vtk");
}

void TPoroCamClayNativo::WriteVTK(int step) {
    fAn->SetStep(step);
    fAn->PostProcess(fPar.vtkResolution);
    fPost->TransferSolution();
    fPost->SetStep(step);
    fPost->PostProcess(fPar.vtkResolution);
}
