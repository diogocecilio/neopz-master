#include "SlopeModel.h"

#include "TPZBndCondT.h"
#include "TPZNullMaterial.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzmultiphysicselement.h"
#include "pzquad.h"

namespace slope {

REAL ContourPressure(int id, REAL y, const Data &d)
{
    switch (id) {
    case bcSuperficie: return d.gammaw * std::max(REAL(0), d.Nivel() - y);
    case bcEsquerda:
    case bcEsquerdaP: return d.gammaw * (d.yCrista - y);
    case bcDireita:
    case bcDireitaP: return d.gammaw * (d.Nivel() - y);
    default: DebugStop();
    }
    return 0.;
}

REAL ExcessPorePressure(REAL p, REAL y, const Data &d) { return p - d.gammaw * (d.yCrista - y); }

TMC MohrCoulomb(const Data &d)
{
    TPZElasticResponse ER;
    ER.SetEngineeringData(d.E, d.nu);
    return TMC(TPZYCMohrCoulombPV2(d.phi, d.phi, d.c, ER), ER);
}

TPZCompMesh *CreateAtomicCompMesh(TPZGeoMesh *gmesh, int nstate, int ordem, const std::set<int> &bcids,
                                  const std::vector<int> &ordens)
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(ordem);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(matSolo, 2, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    if (!ordens.empty()) {
        for (TPZCompEl *cel : cmesh->ElementVec()) {
            if (!cel || cel->Reference()->MaterialId() != matSolo) continue;
            const int o = ordens[cel->Reference()->Index()];
            if (o != ordem) dynamic_cast<TPZInterpolatedElement *>(cel)->PRefine(o); // ajusta os lados com os vizinhos
        }
        cmesh->ExpandSolution();
    }
    return cmesh;
}

TPZMultiphysicsCompMesh *CreateUPCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d, const std::vector<int> &ordens)
{
    using B = TPZMatPoroElastoPlasticUPBase;
    const std::set<int> bcids = {bcBase, bcDireita, bcSuperficie, bcEsquerda, bcDireitaP, bcEsquerdaP};
    TPZManVector<TPZCompMesh *, 2> meshvec = {CreateAtomicCompMesh(gmesh, 2, porder, bcids, ordens), // u
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids)};    // p_w
    auto *mat = new MatUP(matSolo, B::EPlaneStrain);
    mat->SetPlasticModel(MohrCoulomb(d));
    mat->SetBiot(1., 0.);       // graos e agua incompressiveis
    mat->SetPermeability(d.kv); // mobilidade k / gamma_w
    mat->SetBodyForce(TPZManVector<REAL, 3>({0., -d.gamma, 0.}));
    mat->SetFluidWeight(TPZManVector<REAL, 3>({0., -d.gammaw, 0.}));

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2, 2, 0.), sel(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    sel(0, 0) = 1.; // so u_x
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcBase, B::EDirichletU, val1, val2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcDireita, B::EDirichletUDirectional, sel, val2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcEsquerda, B::EDirichletUDirectional, sel, val2));
    for (int id : {bcSuperficie, bcDireitaP, bcEsquerdaP}) { // a tracao nula da superficie e natural
        auto *bc = mat->CreateBC(mat, id, B::EDirichletP, val1, val2);
        bc->SetForcingFunctionBC([id, d](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
            v2[0] = ContourPressure(id, x[1], d);
        });
        mphys->InsertMaterialObject(bc);
    }
    TPZManVector<int, 2> active(2, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matSolo}, bcids);
    return mphys;
}

TPZMultiphysicsCompMesh *CreateUnsatCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d, REAL unitWeight,
                                             const std::vector<int> &ordens)
{
    const std::set<int> bcids = {bcBase, bcDireita, bcSuperficie, bcEsquerda};
    TPZManVector<TPZCompMesh *, 3> meshvec = {CreateAtomicCompMesh(gmesh, 2, porder, bcids, ordens), // u
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids),     // p_w
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids)};    // p_a
    auto *mat = new MatUnsat(matSolo);
    mat->SetPlasticModel(MohrCoulomb(d));
    TPZUnsaturatedHydraulics h; // Gardner com a_G = 0: S_w = 1 e k_rw = 1 em qualquer pressao
    h.fRetention = TPZUnsaturatedHydraulics::EGardner;
    h.fAG = 0.;
    h.fSwr = 0.;
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, 1.); // Terzaghi
    mat->SetPorosity(0.45);
    mat->SetBiot(1., 0.);
    mat->SetWaterBulkModulus(1.e30);          // sem armazenamento: percolacao permanente
    mat->SetConductivity(d.kv * d.gammaw, d.gammaw, 1.e-3 / 1.8e-5);
    mat->SetAir(1.2 * 9.81e-3, 101.325);
    mat->SetUnitWeight(unitWeight);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    for (int id : {bcBase, bcDireita, bcEsquerda}) // os apoios sao impostos pela analise; so marcam o contorno
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, val1, val2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcSuperficie, TUnsat::ENeumannU, val1, val2));
    TPZManVector<int, 3> active(3, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matSolo}, bcids);
    return mphys;
}

void EraseMesh(TPZMultiphysicsCompMesh *mphys)
{
    TPZVec<TPZCompMesh *> meshvec = mphys->MeshVector();
    delete mphys;
    for (TPZCompMesh *m : meshvec) delete m;
}

std::map<int64_t, REAL> BoundaryPressures(TPZGeoMesh *gmesh, const Data &d)
{
    std::map<int64_t, REAL> pb;
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel) continue;
        const int id = gel->MaterialId();
        if (id != bcSuperficie && id != bcEsquerda && id != bcDireita) continue;
        for (int i = 0; i < gel->NCornerNodes(); i++)
            pb[gel->NodeIndex(i)] = ContourPressure(id, gel->NodePtr(i)->Coord(1), d);
    }
    return pb;
}

int RuleOrder(TPZCompEl *cel)
{
    auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
    auto *mat = dynamic_cast<TPZMatCombinedSpaces *>(cel->Material());
    if (!mf || !mat) DebugStop();
    TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
    for (int64_t i = 0; i < mf->NMeshes(); i++)
        if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
    return mat->IntegrationRuleOrder(ordens);
}

void ForEachPoint(TPZMultiphysicsCompMesh *cmesh,
                  const std::function<void(TPZCompEl *, const TPZVec<REAL> &, int64_t)> &f)
{
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZGeoEl *gel = cel->Reference();
        const int o = RuleOrder(cel); // a regra de TPZMultiphysicsCompEl::CalcStiff
        std::unique_ptr<TPZIntPoints> rule(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
        TPZManVector<int, 2> ordem(2, o);
        rule->SetOrder(ordem);
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        if (mem.size() != rule->NPoints()) DebugStop(); // a memoria segue a ordem dos pontos da regra
        for (int ip = 0; ip < rule->NPoints(); ip++) {
            TPZManVector<REAL, 3> qsi(2), x(3);
            REAL w;
            rule->Point(ip, qsi, w);
            gel->X(qsi, x);
            f(cel, x, mem[ip]);
        }
    }
}

REAL SqrtJ2Plastic(const TPZElastoPlasticMem &mem)
{
    TPZTensor<REAL> ep = mem.m_elastoplastic_state.m_eps_p;
    ep.XY() *= 0.5; // cisalhamento de engenharia -> tensor
    ep.XZ() *= 0.5;
    ep.YZ() *= 0.5;
    return std::sqrt(std::max<REAL>(ep.J2(), 0.));
}

} // namespace slope
