#include "TPZMatKLKernel.h"

#include "TPZMatKLKernel.h"

#include "pzcompel.h"
#include "pzgeoel.h"
#include "tpzintpoints.h"
#include "TPZElementMatrixT.h"   // <<-- cabeçalho certo do ElementMatrix TEMPLATED
#include <cmath>
#include <memory>
#include <vector>

// ctor default seguro (id -1), dim=2, kernel exponencial com Lx=1, Ly=0.1
TPZMatKLKernel::TPZMatKLKernel()
: TBase(-1), fDim(2), fTarget(Target::A), fKind(EKernelKind::ExpSeparable)
{
    fParams.Resize(3);
    fParams[0] = 1.0;   // Lx
    fParams[1] = 0.1;   // Ly
    fParams[2] = 0.0;
    RebuildKernel();
}

// ctor com id e dimensão, mantém kernel padrão Exp(Lx=1, Ly=0.1)
TPZMatKLKernel::TPZMatKLKernel(int matid, int dim)
: TBase(matid), fDim(dim), fTarget(Target::A), fKind(EKernelKind::ExpSeparable)
{
    fParams.Resize(3);
    fParams[0] = 1.0;
    fParams[1] = 0.1;
    fParams[2] = 0.0;
    RebuildKernel();
}

// ctor com id, dimensão e parâmetros do kernel exponencial
TPZMatKLKernel::TPZMatKLKernel(int matid, int dim, REAL Lx, REAL Ly)
: TBase(matid), fDim(dim), fTarget(Target::A), fKind(EKernelKind::ExpSeparable)
{
    fParams.Resize(3);
    // garante positividade mínima
    fParams[0] = (Lx > 0 ? Lx : 1.0);
    fParams[1] = (Ly > 0 ? Ly : 1.0);
    fParams[2] = 0.0;
    RebuildKernel();         // monta fKernel a partir de fKind+fParams
}

void TPZMatKLKernel::Print(std::ostream &out) const {
    out << " dim=" << fDim
        //<< " qx=" << fQx << " qy=" << fQy
        << " target=" << (fTarget==Target::A ? "A(C)" : "B(M)") << "\n";
}

void TPZMatKLKernel::SetExpKernel(REAL Lx, REAL Ly)
{
    fKind = EKernelKind::ExpSeparable;
    if (fParams.size() < 3) fParams.Resize(3, 0.);
    fParams[0] = (Lx > 0 ? Lx : 1.);
    fParams[1] = (Ly > 0 ? Ly : 1.);
    fParams[2] = 0.;
    RebuildKernel();
}

void TPZMatKLKernel::SetGaussKernel(REAL Lx, REAL Ly)
{
    fKind = EKernelKind::GaussSeparable;
    if (fParams.size() < 3) fParams.Resize(3, 0.);
    fParams[0] = (Lx > 0 ? Lx : 1.);
    fParams[1] = (Ly > 0 ? Ly : 1.);
    fParams[2] = 0.;
    RebuildKernel();
}

void TPZMatKLKernel::RebuildKernel()
{
    switch (fKind) {

    case EKernelKind::ExpSeparable: {
        const REAL Lx = (fParams.size() > 0 && fParams[0] > 0) ? fParams[0] : 1.;
        const REAL Ly = (fParams.size() > 1 && fParams[1] > 0) ? fParams[1] : 1.;
        fKernel = [Lx, Ly](const TPZVec<REAL>& x, const TPZVec<REAL>& y) -> STATE {
            return (STATE)std::exp(-std::fabs(x[0]-y[0]) / Lx
                                   -std::fabs(x[1]-y[1]) / Ly);
        };
        break;
    }

    case EKernelKind::GaussSeparable: {
        const REAL Lx = (fParams.size() > 0 && fParams[0] > 0) ? fParams[0] : 1.;
        const REAL Ly = (fParams.size() > 1 && fParams[1] > 0) ? fParams[1] : 1.;
        fKernel = [Lx, Ly](const TPZVec<REAL>& x, const TPZVec<REAL>& y) -> STATE {
            const REAL dx = (x[0]-y[0]) / Lx;
            const REAL dy = (x[1]-y[1]) / Ly;
            return (STATE)std::exp(-(dx*dx + dy*dy));
        };
        break;
    }

    case EKernelKind::User: {
        // nada a fazer; assume que fKernel foi setado por SetKernel(...)
        if (!fKernel) {
            // fallback seguro
            fKernel = [](const TPZVec<REAL>&, const TPZVec<REAL>&)->STATE { return (STATE)0; };
        }
        break;
    }

    default: {
        // fallback seguro
        fKernel = [](const TPZVec<REAL>&, const TPZVec<REAL>&)->STATE { return (STATE)0; };
        break;
    }
    }
}
namespace {
/// Dados de quadratura de um elemento: pontos x, pesos w·detJ e funções de forma (n x np)
struct TKLQuadrature {
    std::vector<TPZManVector<REAL,3>> x;
    std::vector<REAL> wJ;
    TPZFNMatrix<200, STATE> phi;
};

/// Calcula os dados de quadratura do elemento com a regra padrão acrescida de extraOrder.
void KLQuadrature(TPZInterpolationSpace *el, int extraOrder, TKLQuadrature &q)
{
    TPZMaterialDataT<STATE> data;
    el->InitMaterialData(data);
    std::unique_ptr<TPZIntPoints> rule(el->GetIntegrationRule().Clone());
    if (extraOrder > 0) {
        TPZManVector<int,3> ord(el->Dimension(), 0);
        rule->GetOrder(ord);
        for (int d = 0; d < ord.size(); d++) ord[d] += extraOrder;
        rule->SetOrder(ord);
    }
    const int np = rule->NPoints();
    const int n = el->NShapeF();
    q.x.resize(np);
    q.wJ.resize(np);
    q.phi.Redim(n, np);
    TPZManVector<REAL,3> qsi(el->Dimension(), 0.);
    REAL w = 0.;
    for (int ip = 0; ip < np; ++ip) {
        rule->Point(ip, qsi, w);
        el->ComputeRequiredData(data, qsi);
        q.x[ip] = data.x;
        q.wJ[ip] = w * data.detjac;
        for (int i = 0; i < n; ++i) q.phi(i, ip) = data.phi(i, 0);
    }
}
} // namespace

/** C_ij = ∬ phi_i(x) K(x,y) phi_j(y) dx dy, integrada com as regras dos dois elementos (produto tensorial
 *  das quadraturas). Os dados de quadratura de cada elemento são calculados uma única vez por par (antes,
 *  ComputeRequiredData de ely era chamado npx*npy vezes) e, no bloco diagonal (elx == ely), onde o kernel
 *  exponencial tem a quina em x = y, a ordem da regra é aumentada de fDiagonalExtraOrder. */
void TPZMatKLKernel::CalcStiffNystrom(TPZInterpolationSpace* elx,
                                      TPZInterpolationSpace* ely,
                                      TPZElementMatrixT<STATE> &ce) const
{
    CalcStiffGalerkin(elx, ely, ce);
}

void TPZMatKLKernel::CalcStiffGalerkin(TPZInterpolationSpace* elx,
                                       TPZInterpolationSpace* ely,
                                       TPZElementMatrixT<STATE> &ce) const
{
    const int extra = (elx == ely) ? fDiagonalExtraOrder : 0;
    TKLQuadrature qx, qy;
    KLQuadrature(elx, extra, qx);
    if (elx == ely) {
        qy.x = qx.x;
        qy.wJ = qx.wJ;
        qy.phi = qx.phi;
    } else {
        KLQuadrature(ely, extra, qy);
    }
    const int nx = elx->NShapeF();
    const int ny = ely->NShapeF();
    const int npx = (int)qx.wJ.size();
    const int npy = (int)qy.wJ.size();
    ce.fMat.Redim(nx, ny);
    ce.fMat.Zero();
    // W(ipx, ipy) = wJx wJy K(x, y);  C = PhiX W PhiY^T
    TPZFNMatrix<400, STATE> W(npx, npy, 0.), PW;
    for (int ipx = 0; ipx < npx; ++ipx)
        for (int ipy = 0; ipy < npy; ++ipy)
            W(ipx, ipy) = (STATE)(qx.wJ[ipx] * qy.wJ[ipy]) * (fKernel ? fKernel(qx.x[ipx], qy.x[ipy]) : (STATE)0);
    qx.phi.Multiply(W, PW);  // PW = PhiX W (nx x npy)
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            STATE v = 0.;
            for (int ipy = 0; ipy < npy; ++ipy) v += PW(i, ipy) * qy.phi(j, ipy);
            ce.fMat(i, j) = v;
        }
}

/** B_ij = ∫ phi_i phi_j dx (massa consistente) */
void TPZMatKLKernel::CalcStiffMass(TPZInterpolationSpace* el,
                                   TPZElementMatrixT<STATE> &be,
                                   int qmass) const
{
    TPZMaterialDataT<STATE> data;
    el->InitMaterialData(data);

    std::unique_ptr<TPZIntPoints> rule(el->GetIntegrationRule().Clone());

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
    // buf.Write(&fDim, 1);
    // //buf.Write(&fQx, 1);
    // //buf.Write(&fQy, 1);
    // int targ = (fTarget==Target::A ? 0 : 1);
    // buf.Write(&targ, 1);
    // // fKernel não é serializado (functor).
}

void TPZMatKLKernel::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    // buf.Read(&fDim, 1);
    // //buf.Read(&fQx, 1);
    // //buf.Read(&fQy, 1);
    // int targ=0; buf.Read(&targ, 1);
    // fTarget = (targ==0 ? Target::A : Target::B);
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
// --- .cpp ---
// acrescente os novos ids (mantive ordem e adicionei SolutionSquared sem quebrar os existentes)
enum EVarIds {
    EVar_Solution = 1, EVar_ExactSolution, EVar_Error, EVar_ErrorSquared,
    EVar_Gradient, EVar_ExactGradient, EVar_ErrorGrad, EVar_GradErrorSquared,

    // novo para normalização por massa
    EVar_SolutionSquared,

    // pós-processo do campo estocástico carregado como Solution
    EVar_KLField, EVar_KLFieldGrad, EVar_KLFieldLog
};

int TPZMatKLKernel::VariableIndex(const std::string& name) const {
    if (name=="Solution")           return EVar_Solution;
    if (name=="ExactSolution")      return EVar_ExactSolution;
    if (name=="Error")              return EVar_Error;
    if (name=="ErrorSquared")       return EVar_ErrorSquared;
    if (name=="Gradient")           return EVar_Gradient;
    if (name=="ExactGradient")      return EVar_ExactGradient;
    if (name=="ErrorGrad")          return EVar_ErrorGrad;
    if (name=="GradErrorSquared")   return EVar_GradErrorSquared;

    // **necessário para BuildPhiSqrtLambdaNodal**
    if (name=="SolutionSquared")    return EVar_SolutionSquared;

    // pós-processo do campo
    if (name=="KLField")            return EVar_KLField;
    if (name=="KLFieldGrad")        return EVar_KLFieldGrad;
    if (name=="KLFieldLog")         return EVar_KLFieldLog;

    return -1;
}

int TPZMatKLKernel::NSolutionVariables(int var) const {
    switch (var) {
        case EVar_Solution:
        case EVar_ExactSolution:
        case EVar_Error:
        case EVar_ErrorSquared:
        case EVar_GradErrorSquared:
        case EVar_SolutionSquared:   // <- novo
        case EVar_KLField:
        case EVar_KLFieldLog:
            return 1;

        case EVar_Gradient:
        case EVar_ExactGradient:
        case EVar_ErrorGrad:
        case EVar_KLFieldGrad:
            return this->Dimension(); // tipicamente 2
    }
    return 0;
}

void TPZMatKLKernel::Solution(const TPZMaterialDataT<STATE>& data, int var,
                              TPZVec<STATE>& out)
{
    // solução primal no ponto
    const STATE uh = data.sol[0][0];

    // solução exata (se fornecida)
    STATE uex = 0.;
    TPZFMatrix<STATE> duex; duex.Redim(2,1); duex.Zero();
    if (fExact) fExact(data.x, uex, duex);

    const STATE e  = uh - uex;
    const STATE e2 = e*e;

    // gradiente numérico da solução
    const int dim = this->Dimension();
    TPZManVector<STATE,3> gradh(3,0.);
    if (data.dsol.size()) {
        for (int i=0; i<dim; i++) gradh[i] = data.dsol[0](i,0);
    }
    const STATE ex = gradh[0] - duex(0,0);
    const STATE ey = (dim>1 ? gradh[1] - duex(1,0) : 0.);
    const STATE gradE2 = ex*ex + ey*ey;

    switch (var) {
        case EVar_Solution:         out[0] = uh;                break;
        case EVar_ExactSolution:    out[0] = uex;               break;
        case EVar_Error:            out[0] = e;                 break;
        case EVar_ErrorSquared:     out[0] = e2;                break;

        case EVar_Gradient:
            for (int i=0;i<dim;i++) out[i] = gradh[i];
            break;

        case EVar_ExactGradient:
            out[0] = duex(0,0); if (dim>1) out[1] = duex(1,0);
            break;

        case EVar_ErrorGrad:
            out[0] = ex; if (dim>1) out[1] = ey;
            break;

        case EVar_GradErrorSquared:
            out[0] = gradE2;
            break;

            // **novo**: usado na normalização por massa
        case EVar_SolutionSquared:
            out[0] = uh*uh;
            break;

            // pós-processo do campo (quando você carregar E(x) como Solution)
        case EVar_KLField:
            out[0] = uh;
            break;

        case EVar_KLFieldGrad:
            for (int i=0;i<dim;i++) out[i] = gradh[i];
            break;

        case EVar_KLFieldLog: {
            constexpr STATE eps = (STATE)1e-30;
            const STATE val = uh > eps ? uh : eps;
            out[0] = std::log(val);
        } break;

        default:
            out.Fill(0.);
    }
}

// #include "TPZSavable.h"          // já deve existir em seu include path
//
// // --- registra a classe para o PersistenceManager ---
// template class TPZRestoreClass<TPZMatKLKernel>;
// static TPZRestoreClass<TPZMatKLKernel> gTPZMatKLKernelRestore;
// namespace
template class TPZRestoreClass<TPZMatKLKernel>;
