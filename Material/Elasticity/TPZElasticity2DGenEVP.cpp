#include "TPZElasticity2DGenEVP.h"

TPZElasticity2DGenEVP::TPZElasticity2DGenEVP(int id, REAL E, REAL nu,
                                             REAL fx, REAL fy,
                                             REAL rho, REAL thickness)
    : TPZElasticity2D(id, E, nu, fx, fy),
      fRho(rho),
      fThick(thickness),
      fWhich(EMatrixA_K)
{
}

int TPZElasticity2DGenEVP::ClassId() const
{
    return Hash("TPZElasticity2DGenEVP") ^ (TPZElasticity2D::ClassId() << 1);
}

void TPZElasticity2DGenEVP::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity2D::Write(buf, withclassid);
    buf.Write(fRho);
    buf.Write(fThick);
    int which = static_cast<int>(fWhich);
    buf.Write(&which, 1);
}

void TPZElasticity2DGenEVP::Read(TPZStream &buf, void *context)
{
    TPZElasticity2D::Read(buf, context);
    buf.Read(&fRho);
    buf.Read(&fThick);
    int which = 0;
    buf.Read(&which, 1);
    fWhich = static_cast<EWhichMatrix>(which);
}

void TPZElasticity2DGenEVP::Contribute(const TPZMaterialDataT<STATE> &data,
                                       REAL weight,
                                       TPZFMatrix<STATE> &ek,
                                       TPZFMatrix<STATE> &ef)
{
    if (fWhich == EMatrixA_K) {
        TPZElasticity2D::Contribute(data, weight, ek, ef); // monta K
    } else {
        ContributeMass(data, weight, ek, ef);               // monta M
    }
}

// void TPZElasticity2DGenEVP::ContributeBC(const TPZMaterialDataT<STATE> &data,
//                                          REAL weight,
//                                          TPZFMatrix<STATE> &ek,
//                                          TPZFMatrix<STATE> &ef,
//                                          TPZBndCondT<STATE> &bc)
// {
//     if (fWhich == EMatrixA_K) {
//         TPZElasticity2D::ContributeBC(data, weight, ek, ef, bc);
//     }
// }
void TPZElasticity2DGenEVP::ContributeBCMass(const TPZMaterialDataT<STATE> &data,
                                             REAL weight,
                                             TPZFMatrix<STATE> &ek,
                                             TPZFMatrix<STATE> &ef,
                                             TPZBndCondT<STATE> &bc)
{
    // mesmo padrão do TPZElasticity2D::ContributeBC
    const TPZFMatrix<REAL> &phi = data.fPhi;
    const int phr = phi.Rows();
    const int nstate = 2; // u,v
    const auto &BIGNUMBER = TPZMaterial::fBigNumber;

    switch (bc.Type()) {

        // ---------------- Dirichlet (penalty) ----------------
        case 0: // Dirichlet total: impõe u = v2[0], v = v2[1]
        {
            // Val2 tem o deslocamento prescrito
            TPZManVector<STATE,2> v2(2,0.);
            if (bc.Val2().NElements() >= 2) {
                v2[0] = bc.Val2()[0];
                v2[1] = bc.Val2()[1];
            }
            for (int in = 0; in < phr; in++) {
                // RHS
                ef(nstate*in+0,0) += BIGNUMBER * v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * v2[1] * phi(in,0) * weight;
                // K (na massa B)
                for (int jn = 0; jn < phr; jn++) {
                    const STATE kpen = BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    ek(nstate*in+0, nstate*jn+0) += kpen;
                    ek(nstate*in+1, nstate*jn+1) += kpen;
                }
            }
        }
        break;

        case 3: // Dirichlet direcional (penaliza apenas as comps ligadas)
        {
            TPZManVector<STATE,2> mask(2,1.);
            if (bc.Val2().NElements() >= 2) {
                mask[0] = bc.Val2()[0];
                mask[1] = bc.Val2()[1];
            }
            for (int in = 0; in < phr; in++) {
                for (int jn = 0; jn < phr; jn++) {
                    const STATE kpen = BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    if (std::abs(mask[0]) > 0) ek(nstate*in+0, nstate*jn+0) += kpen;
                    if (std::abs(mask[1]) > 0) ek(nstate*in+1, nstate*jn+1) += kpen;
                }
            }
        }
        break;

    }
}


void TPZElasticity2DGenEVP::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                         REAL weight,
                                         TPZFMatrix<STATE> &ek,
                                         TPZFMatrix<STATE> &ef,
                                         TPZBndCondT<STATE> &bc)
{
    if (fWhich == EMatrixA_K) {
        // BCs usuais (Dirichlet/Neumann) afetam K
        TPZElasticity2D::ContributeBC(data, weight, ek, ef, bc);
    } else {
        // Montando B: só faça algo se houver massa na fronteira
        ContributeBCMass(data, weight, ek, ef, bc);
    }
}

void TPZElasticity2DGenEVP::ContributeMass(const TPZMaterialDataT<STATE> &data,
                                           REAL weight,
                                           TPZFMatrix<STATE> &ek,
                                           TPZFMatrix<STATE> & /*ef*/)
{
    const int nshape = data.phi.Rows();
    const int dim = 2; // 2D
    const int ndof = dim * nshape;

    TPZFMatrix<STATE> N(dim, ndof, 0.0);
    for (int i = 0; i < nshape; i++) {
        STATE phi = data.phi(i,0);
        N(0, 2*i  ) = phi;
        N(1, 2*i+1) = phi;
    }

    TPZFMatrix<STATE> NtN(ndof, ndof, 0.0);
    for (int r = 0; r < ndof; r++) {
        for (int c = 0; c < ndof; c++) {
            NtN(r,c) = N(0,r)*N(0,c) + N(1,r)*N(1,c);
        }
    }

    const STATE coeff = fRho * fThick * weight;
    for (int r = 0; r < ndof; r++) {
        for (int c = 0; c < ndof; c++) {
            ek(r,c) += coeff * NtN(r,c);
        }
    }
}

int TPZElasticity2DGenEVP::VariableIndex(const std::string &name) const {
    if (name == "Ux" || name == "EVP_Ux") return EVP_UX;
    if (name == "Uy" || name == "EVP_Uy") return EVP_UY;
    if (name == "U"  || name == "EVP_U")  return EVP_UMAG;
    if (name == "displacement")           return EVP_UVEC;   // <—
    if (name == "POrder" || name == "p")  return EVP_PORDER;
    DebugStop(); return -1;
}

int TPZElasticity2DGenEVP::NSolutionVariables(int var) const {
    switch (var) {
        case EVP_UX:    return 1;
        case EVP_UY:    return 1;
        case EVP_UMAG:  return 1;
        case EVP_UVEC:  return 2; // <— (ux,uy)
        case EVP_PORDER:return 1;
        default: DebugStop(); return 0;
    }
}

void TPZElasticity2DGenEVP::Solution(const TPZMaterialDataT<STATE> &data,
                                     int var, TPZVec<STATE> &sol)
{
    const auto &u = data.sol[0]; // (ux,uy)
    const STATE ux = (u.size() > 0 ? u[0] : 0.);
    const STATE uy = (u.size() > 1 ? u[1] : 0.);

    switch (var) {
        case EVP_UX:   sol.Resize(1); sol[0] = ux; break;
        case EVP_UY:   sol.Resize(1); sol[0] = uy; break;
        case EVP_UMAG: sol.Resize(1); sol[0] = std::sqrt(ux*ux + uy*uy); break;
        case EVP_UVEC: sol.Resize(2); sol[0] = ux; sol[1] = uy; break; // <—
        case EVP_PORDER: sol.Resize(1); sol[0] = data.p; break;
        default: DebugStop(); break;
    }
}

// multiphysics wrapper
void TPZElasticity2DGenEVP::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                     int var, TPZVec<STATE> &sol)
{
    Solution(datavec[0], var, sol);
}
void TPZElasticity2DGenEVP::FillDataRequirements(TPZMaterialData &data)
{
    TPZElasticity2D::FillDataRequirements(data);
    data.SetAllRequirements(false);
    data.fNeedsSol = true;   // preciso de data.sol[0]
}
