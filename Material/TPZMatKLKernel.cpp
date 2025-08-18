#include "TPZMatKLKernel.h"

#include "pzcompel.h"
#include "pzgeoel.h"
#include "tpzintpoints.h"
#include "TPZElementMatrixT.h"   // <<-- cabeçalho certo do ElementMatrix TEMPLATED
#include <cmath>

TPZMatKLKernel::TPZMatKLKernel(int matid, int dim, KernelFn ker/*, int qx, int qy*/)
: TBase(matid), fKernel(std::move(ker)), fDim(dim)/*, fQx(qx), fQy(qy) */ {}

void TPZMatKLKernel::Print(std::ostream &out) const {
    out << " dim=" << fDim
        //<< " qx=" << fQx << " qy=" << fQy
        << " target=" << (fTarget==Target::A ? "A(C)" : "B(M)") << "\n";
}

/** C_ij = ∬ phi_i(x) K(x,y) phi_j(y) dx dy  (Nyström) */
void TPZMatKLKernel::CalcStiffNystrom(TPZInterpolationSpace* elx,
                                      TPZInterpolationSpace* ely,
                                      TPZElementMatrixT<STATE> &ce) const
{
    TPZMaterialDataT<STATE> dataX, dataY;
    elx->InitMaterialData(dataX);
    ely->InitMaterialData(dataY);

    auto ruleX = elx->GetIntegrationRule().Clone();
    auto ruleY = ely->GetIntegrationRule().Clone();

    // SetOrder usa TPZVec<int> nas versões recentes
    //TPZManVector<int,3> ordx(elx->Dimension(), fQx);
    //TPZManVector<int,3> ordy(ely->Dimension(), fQy);
    //ruleX->SetOrder(ordx);
    //ruleY->SetOrder(ordy);

    const int nx = elx->NShapeF();
    const int ny = ely->NShapeF();
    ce.fMat.Redim(nx, ny);
    ce.fMat.Zero();

    TPZManVector<REAL,3> qsiX(elx->Dimension()), qsiY(ely->Dimension());
    REAL wx = 0., wy = 0.;
    const int npx = ruleX->NPoints();
    const int npy = ruleY->NPoints();

    for (int ipx=0; ipx<npx; ++ipx) {
        ruleX->Point(ipx, qsiX, wx);

        elx->ComputeRequiredData(dataX, qsiX);

        const REAL wJx = wx * dataX.detjac;
        const TPZVec<REAL> &X = dataX.x;

        for (int ipy=0; ipy<npy; ++ipy) {
            ruleY->Point(ipy, qsiY, wy);

            ely->ComputeRequiredData(dataY, qsiY);

            const REAL wJy = wy * dataY.detjac;
            const TPZVec<REAL> &Y = dataY.x;

            const STATE Kxy = fKernel ? fKernel(X,Y) : (STATE)0;
            const STATE W   = (STATE)(wJx*wJy) * Kxy;

            const TPZFMatrix<STATE> &phiX = dataX.phi;
            const TPZFMatrix<STATE> &phiY = dataY.phi;

            for (int i=0;i<nx;i++){
                const STATE pix = phiX(i,0);
                for (int j=0;j<ny;j++){
                    ce.fMat(i,j) += W * pix * (STATE)phiY(j,0);
                }
            }
        }
    }
}

/** B_ij = ∫ phi_i phi_j dx (massa consistente) */
void TPZMatKLKernel::CalcStiffMass(TPZInterpolationSpace* el,
                                   TPZElementMatrixT<STATE> &be,
                                   int qmass) const
{
    TPZMaterialDataT<STATE> data;
    el->InitMaterialData(data);

    auto rule = el->GetIntegrationRule().Clone();
    TPZManVector<int,3> ord(el->Dimension(), qmass);
    rule->SetOrder(ord);

    const int n = el->NShapeF();
    be.fMat.Redim(n,n);
    be.fMat.Zero();

    TPZManVector<REAL,3> qsi(el->Dimension());
    REAL w = 0.;
    const int np = rule->NPoints();

    for (int ip=0; ip<np; ++ip) {
        rule->Point(ip, qsi, w);
        el->ComputeRequiredData(data, qsi);

        const REAL weight = w * data.detjac;

        const TPZFMatrix<STATE> &phi = data.phi;
        for (int i=0;i<n;i++){
            const STATE pi = phi(i,0);
            for (int j=0;j<n;j++){
                be.fMat(i,j) += weight * pi * (STATE)phi(j,0);
            }
        }
    }
}

int TPZMatKLKernel::ClassId() const {
    return Hash("TPZMatKLKernel") ^ (TPZMaterial::ClassId() << 1);
}

void TPZMatKLKernel::Write(TPZStream &buf, int withclassid) const {
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fDim, 1);
    //buf.Write(&fQx, 1);
    //buf.Write(&fQy, 1);
    int targ = (fTarget==Target::A ? 0 : 1);
    buf.Write(&targ, 1);
    // fKernel não é serializado (functor).
}

void TPZMatKLKernel::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fDim, 1);
    //buf.Read(&fQx, 1);
    //buf.Read(&fQy, 1);
    int targ=0; buf.Read(&targ, 1);
    fTarget = (targ==0 ? Target::A : Target::B);
}
bool TPZMatKLKernel::HasForcingFunction() const
{
    return false;
}

void TPZMatKLKernel::FillDataRequirements(TPZMaterialData&data)const
{
    data.SetAllRequirements(false);
    //data.fNeedsSol=false;
    //data.fNeedsNormal=false;
    //data.fNeedsNeighborSol=false;
}

