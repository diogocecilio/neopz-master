// SeepageProblem.cpp — ver SeepageProblem.h

#include "SeepageProblem.h"

#include <algorithm>
#include <cmath>

#include "TPZBndCondT.h"
#include "pzskylstrmatrix.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "pzstack.h"
#include "pzstepsolver.h"
#include "pzquad.h"

void TPZDarcyFlowAnisotropic::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                                         TPZFMatrix<STATE> &ef) {
    // gradientes nas direções globais: dphi_global = axesᵀ dphix
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &axes = data.axes;
    const int nphi = dphi.Cols();
    REAL kv = 1.;
    if (!fKv.empty() && data.gelElId >= 0 && data.gelElId < (int)fKv.size()) kv = fKv[data.gelElId];
    const REAL kxx = fAlpha * kv, kyy = kv;
    TPZFNMatrix<60, REAL> g(2, nphi, 0.);
    for (int in = 0; in < nphi; in++)
        for (int d = 0; d < 2; d++) g(d, in) = axes(0, d) * dphi(0, in) + axes(1, d) * dphi(1, in);
    for (int in = 0; in < nphi; in++)
        for (int jn = 0; jn < nphi; jn++)
            ek(in, jn) += weight * (kxx * g(0, in) * g(0, jn) + kyy * g(1, in) * g(1, jn));
    (void)ef;
}

TSeepageProblem::TSeepageProblem(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TParams &par)
    : fGMesh(gmesh), fGeo(geo), fPar(par) {
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(2);
    fCMesh->SetDefaultOrder(par.porder);
    fCMesh->SetAllCreateFunctionsContinuous();
    fMat = new TPZDarcyFlowAnisotropic(ESoil, 2);
    fMat->SetAnisotropy(par.alpha);
    fCMesh->InsertMaterialObject(fMat);
    // Dirichlet: u = -γw min(y', hw) (crista: y' = 0; face; pé: y' = H >= hw)
    const REAL D = geo.D(), gw = par.gammaW, hw = par.hw;
    ForcingFunctionBCType<STATE> ud = [D, gw, hw](const TPZVec<REAL> &x, TPZVec<STATE> &rhs, TPZFMatrix<STATE> &) {
        const REAL depth = std::max(D - x[1], REAL(0.));
        rhs[0] = -gw * std::min(depth, hw);
    };
    TPZFNMatrix<1, STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 1> val2(1, 0.);
    for (int id : {ECrest, EFace, EToe}) {
        auto *bc = fMat->CreateBC(fMat, id, 0, val1, val2);
        bc->SetForcingFunctionBC(ud);
        fCMesh->InsertMaterialObject(bc);
    }
    std::set<int> matids = {ESoil, ECrest, EFace, EToe};
    fCMesh->AutoBuild(matids);
    fCelOfGel.assign(gmesh->NElements(), nullptr);
    for (TPZCompEl *cel : fCMesh->ElementVec())
        if (cel && cel->Reference() && cel->Reference()->Dimension() == 2) fCelOfGel[cel->Reference()->Index()] = cel;

    fAn = std::make_unique<TPZLinearAnalysis>(fCMesh, true);
    TPZSkylineStructMatrix<STATE> skl(fCMesh);
    skl.SetNumThreads(0);
    fAn->SetStructuralMatrix(skl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    fAn->SetSolver(step);
}

TSeepageProblem::~TSeepageProblem() {
    fAn.reset();
    delete fCMesh;
}

void TSeepageProblem::SetElementPermeability(const std::vector<REAL> &kvByIndex) {
    // o material lê kv pelo identificador do elemento geométrico (TPZMaterialData::gelElId)
    int64_t maxid = 0;
    for (TPZGeoEl *gel : fGMesh->ElementVec())
        if (gel) maxid = std::max<int64_t>(maxid, gel->Id());
    std::vector<REAL> kvById(maxid + 1, 1.);
    for (TPZGeoEl *gel : fGMesh->ElementVec())
        if (gel && gel->Index() < (int64_t)kvByIndex.size()) kvById[gel->Id()] = kvByIndex[gel->Index()];
    fMat->SetElementPermeability(kvById);
}

void TSeepageProblem::Solve() {
    fAn->Assemble();
    fAn->Solve();
}

void TSeepageProblem::Evaluate(int64_t gelIndex, const TPZVec<REAL> &qsi, REAL &u, TPZManVector<REAL, 2> &gradu) {
    auto *intel = dynamic_cast<TPZInterpolationSpace *>(fCelOfGel[gelIndex]);
    TPZMaterialDataT<STATE> data;
    TPZManVector<REAL, 3> q(qsi);
    intel->ComputeSolution(q, data, false);
    u = data.sol[0][0];
    // dsol nos eixos locais do elemento -> componentes globais
    gradu.Resize(2);
    for (int d = 0; d < 2; d++) gradu[d] = data.axes(0, d) * data.dsol[0](0, 0) + data.axes(1, d) * data.dsol[0](1, 0);
}

REAL TSeepageProblem::Functional() {
    REAL J = 0.;
    for (TPZCompEl *cel : fCMesh->ElementVec()) {
        auto *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel || !intel->Reference() || intel->Reference()->Dimension() != 2) continue;
        TPZGeoEl *gel = intel->Reference();
        std::unique_ptr<TPZIntPoints> rule(gel->CreateSideIntegrationRule(gel->NSides() - 1, 2 * fPar.porder + 2));
        REAL kv = 1.;
        // a mesma permeabilidade de TPZDarcyFlowAnisotropic::Contribute (por identificador do elemento)
        if (const std::vector<REAL> &kvs = fMat->ElementPermeability(); !kvs.empty() && gel->Id() < (int64_t)kvs.size())
            kv = kvs[gel->Id()];
        TPZManVector<REAL, 3> qsi(2, 0.);
        for (int ip = 0; ip < rule->NPoints(); ip++) {
            REAL w = 0.;
            rule->Point(ip, qsi, w);
            TPZFNMatrix<9, REAL> jac(2, 2), axes(2, 3), jacinv(2, 2);
            REAL detjac = 0.;
            gel->Jacobian(qsi, jac, axes, detjac, jacinv);
            REAL u = 0.;
            TPZManVector<REAL, 2> g(2, 0.);
            Evaluate(gel->Index(), qsi, u, g);
            J += 0.5 * w * std::fabs(detjac) * kv * (fPar.alpha * g[0] * g[0] + g[1] * g[1]);
        }
    }
    return J;
}

void TSeepageProblem::DefineVTK(const std::string &file) {
    TPZStack<std::string> scal, vec;
    scal.Push("Pressure");
    vec.Push("Derivative");
    fAn->DefineGraphMesh(2, scal, vec, file);
    fVTKDefined = true;
}

void TSeepageProblem::WriteVTK(int step) {
    if (!fVTKDefined) return;
    fAn->SetStep(step);
    fAn->PostProcess(1);
}
