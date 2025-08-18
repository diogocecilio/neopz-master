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
    // TPZManVector<int,3> ord(el->Dimension(), qmass);
    // rule->SetOrder(ord);

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
    data.fNeedsSol=true;
    //data.fNeedsNormal=false;
    //data.fNeedsNeighborSol=false;
}
int TPZMatKLKernel::VariableIndex(const std::string& name) const {
  if (name=="Solution"||name=="u") return ESolution;
  if (name=="Gradient"||name=="grad") return EGradient;
  if (name=="ExactSolution") return EExact;
  if (name=="ExactGradient") return EExactGrad;
  if (name=="Error") return EError;
  if (name=="ErrorGrad") return EErrorGrad;
  return -1;
}

int TPZMatKLKernel::NSolutionVariables(int var) const {
  switch (var) {
    case ESolution: case EExact: case EError: return 1;
    case EGradient: case EExactGrad: case EErrorGrad: return Dimension();
    default: return 0;
  }
}

void TPZMatKLKernel::Solution(const TPZMaterialDataT<STATE>& data, int var,
                              TPZVec<STATE>& Solout) {
  const int dim = Dimension();
  switch (var) {
    case ESolution: {
      Solout.Resize(1); Solout[0] = data.sol[0][0]; break;
    }
    case EGradient: {
      Solout.Resize(dim);
      for (int i=0;i<dim;i++) Solout[i] = data.dsol[0](i,0);
      break;
    }
    case EExact:
    case EExactGrad:
    case EError:
    case EErrorGrad: {
      STATE ue = 0; TPZFMatrix<STATE> due(dim,1,0.);
      if (fExact) fExact(data.x, ue, due);
      if (var==EExact) { Solout.Resize(1); Solout[0] = ue; break; }
      if (var==EExactGrad) {
        Solout.Resize(dim); for (int i=0;i<dim;i++) Solout[i]=due(i,0); break;
      }
      if (var==EError) { // u_h - u_ex
        Solout.Resize(1); Solout[0] = data.sol[0][0] - ue; break;
      }
      if (var==EErrorGrad) {
        Solout.Resize(dim);
        for (int i=0;i<dim;i++) Solout[i] = data.dsol[0](i,0) - due(i,0);
        break;
      }
      break;
    }
    default: Solout.Resize(0); break;
  }
}

