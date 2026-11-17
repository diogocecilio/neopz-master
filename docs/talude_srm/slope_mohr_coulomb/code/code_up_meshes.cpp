/// Toe ground (-3) and face (-6) need a pore pressure and a water load: twin boundary elements -13 and -16
void AddWaterLoadElements(TPZGeoMesh *gmesh) {
    const int64_t nel = gmesh->NElements();
    for (int64_t i = 0; i < nel; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (gel && !gel->HasSubElement() && (gel->MaterialId() == -3 || gel->MaterialId() == -6))
            TPZGeoElBC(gel, gel->NSides() - 1, gel->MaterialId() - 10);
    }
}

/// Atomic H1 space of the u-p mesh with a null material: nstate = 2 displacement, 1 pore pressure
TPZCompMesh *CreateAtomicMesh(TPZGeoMesh *gmesh, int nstate, int order, const std::set<int> &bcids) {
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(1, 2, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

/// Biot u-p mesh of the slope (plane strain, Taylor-Hood: u P2 and p P1, porder = 2).
/// Base fixed, sides on rollers, base and sides impermeable. Reservoir level H = H0 - lambda (H0 - H1), lambda =
/// load factor of the step (0: crest, 1: drawn down; bisection interpolates it). Ground surface (toe -3, crest -4,
/// face -6): p = gamma_w (H - y)+, i.e. drained above the water (the phreatic surface stays at the crest, as in
/// the rapid drawdown of Ceron et al. 2025); water load t = -gamma_w (H - y)+ n on the twins -13 and -16.
template <class TPlastic>
TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, int porder, const TPlastic &model, const Soil &s,
                                        const Water &w) {
    using B = TPZMatPoroElastoPlasticUPBase;
    if (porder != 2) DebugStop(); // PoreField and the nodal Dirichlet values assume a P1 pressure
    const std::set<int> bcids = {-1, -2, -3, -4, -5, -6, -13, -16};
    TPZManVector<TPZCompMesh *, 2> meshvec = {CreateAtomicMesh(gmesh, 2, porder, bcids),
                                              CreateAtomicMesh(gmesh, 1, porder - 1, bcids)};
    auto *mat = new TPZMatPoroElastoPlasticUP<TPlastic, TPZElastoPlasticMem>(1, B::EPlaneStrain);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., 0.);                  // incompressible grains and water
    mat->SetPermeability(w.kh / w.gammaw); // mobility
    TPZManVector<REAL, 3> b = {0., -s.gamma, 0.}, gw = {0., -w.gammaw, 0.};
    mat->SetBodyForce(b); // saturated unit weight
    mat->SetFluidWeight(gw);
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);
    auto pw = [mat, w](const TPZVec<REAL> &x) {
        const REAL H = w.H0 - mat->LoadFactor() * (w.H0 - w.H1);
        return w.gammaw * std::max<REAL>(H - x[1], 0.);
    };
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    mphys->InsertMaterialObject(mat->CreateBC(mat, -1, B::EDirichletU, val1, val2));
    for (int id : {-3, -4, -6}) {
        auto *bc = mat->CreateBC(mat, id, B::EDirichletP, val1, val2);
        bc->SetForcingFunctionBC([pw](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) { v[0] = pw(x); });
        mphys->InsertMaterialObject(bc);
    }
    const REAL normal[2][2] = {{0., 1.}, {M_SQRT1_2, M_SQRT1_2}}; // outward normals of -13 (toe) and -16 (face)
    for (int i = 0; i < 2; i++) {
        auto *bc = mat->CreateBC(mat, i ? -16 : -13, B::ENeumannUFixed, val1, val2);
        const REAL nx = normal[i][0], ny = normal[i][1];
        bc->SetForcingFunctionBC([pw, nx, ny](const TPZVec<REAL> &x, TPZVec<STATE> &v, TPZFMatrix<STATE> &) {
            v[0] = -pw(x) * nx;
            v[1] = -pw(x) * ny;
        });
        mphys->InsertMaterialObject(bc);
    }
    val1(0, 0) = 1.; // u_x = 0
    for (int id : {-2, -5}) mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    TPZManVector<int, 2> active(2, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {1}, bcids);
    return mphys;
}
