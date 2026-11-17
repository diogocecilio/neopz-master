/// Seepage body force of the stability analysis, b = lambda (gamma_sat g - grad p+): lambda is the gravity factor
/// that SlopeAnalysis sets through the body force of the material (0, -lambda gamma_sat, 0)
template <class TPlastic>
void SetSeepageForce(TPZCompMesh *cmesh, const PoreField &pf, REAL gamma) {
    auto *mat = dynamic_cast<TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem> *>(cmesh->FindMaterial(1));
    if (!mat) DebugStop();
    mat->SetForcingFunction(
        [mat, &pf, gamma](const TPZVec<REAL> &x, TPZVec<STATE> &f) {
            const REAL lambda = -mat->GetBodyForce()[1] / gamma;
            REAL g[2];
            pf.GradPositive(x, g);
            f[0] = -lambda * g[0];
            f[1] = -lambda * (gamma + g[1]);
            f[2] = 0.;
        },
        0);
}
