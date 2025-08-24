#include "TPZMatElastic2DMem.h"
#include "TPZMaterialDataT.h"
#include "pzerror.h"
#include <cmath>
#include <algorithm>

// ---------- ctor ----------
template <class TMEM>
TPZMatElastic2DMem<TMEM>::TPZMatElastic2DMem(int matid, EPlaneType plane)
: TBase(matid), fPlane(plane)
{
}

// ---------- requisitos ----------
template <class TMEM>
void TPZMatElastic2DMem<TMEM>::FillDataRequirements(TPZMaterialData &data) const {
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fDeformedDirections = true;
    data.fNeedsNormal = false;
    data.fNeedsHSize = false;
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::FillBoundaryConditionDataRequirements(int, TPZMaterialData &data) const {
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

// ---------- contribuição volumétrica ----------
template <class TMEM>
void TPZMatElastic2DMem<TMEM>::Contribute(const TPZMaterialDataT<STATE> &data,
                                          REAL weight,
                                          TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &phi = data.fPhi;
    const TPZFMatrix<REAL> &jac = data.jacobian;
    TPZFMatrix<REAL> dpsi;
    (void)ef;
    if (data.phi.Rows() == 0) return;

    // parâmetros elásticos do IP (da memória ou fallback)
    REAL E, nu;
    ERFromMem(data, E, nu);

    TPZFMatrix<REAL> D(3,3,0.0);
    BuildConstitutiveMatrix(E, nu, D);

    const int nshape = phi.Rows();
    const int ndofel = 2 * nshape;

    // B (3 x 2n)
    TPZFMatrix<REAL> B(3, ndofel, 0.0);
    BuildBMatrix(dphi, B);


    TPZFMatrix<REAL> Bt;
    B.Transpose(&Bt);

    TPZFMatrix<REAL> BtD;
    Bt.Multiply(D, BtD);

    TPZFMatrix<REAL> Kloc;
    BtD.Multiply(B, Kloc);

    // Kloc = (B^T D B) * detjac * weight
    ek +=  weight*Kloc;


}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::Contribute(const TPZMaterialDataT<STATE> &data,
                                          REAL weight, TPZFMatrix<REAL> &ef)
{
    (void)data; (void)weight; (void)ef;
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                                            TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
                                            TPZBndCondT<STATE> &bc)
{
    const int nshape = data.phi.Rows();
    const auto &phi = data.phi;
    const int type = bc.Type();
    constexpr STATE big = (STATE)1.e12;
    if (type == 0) { // Dirichlet (penalidade)

        const auto &val2 = bc.Val2(); // vetor [u_Dx, u_Dy]
        for (int i=0;i<nshape;i++){
            for (int j=0;j<nshape;j++){
                const STATE k = big * phi(i,0) * phi(j,0) * weight;
                ek(2*i+0,2*j+0) += k;
                ek(2*i+1,2*j+1) += k;
            }
            ef(2*i+0) += big * phi(i,0) * (val2.size()>0 ? val2[0] : (STATE)0.) * weight;
            ef(2*i+1) += big * phi(i,0) * (val2.size()>1 ? val2[1] : (STATE)0.) * weight;
        }
    } else if (type == 1) { // Neumann: tração
        const auto &val2 = bc.Val2(); // vetor [tx, ty]
        for (int i=0;i<nshape;i++){
            ef(2*i+0) += phi(i,0) * val2[0]  * weight;
            ef(2*i+1) += phi(i,0) * val2[1]  * weight;
        }
    } else if (type == 3) { // Directional Dirichlet
        for(int in = 0 ; in < nshape; in++) {
            const auto &v2 = bc.Val2();
            for (int jn = 0 ; jn < nshape; jn++) {
                ek(2*in+0,2*jn+0) += big * phi(in,0) * phi(jn,0) * weight * v2[0];
                ek(2*in+1,2*jn+1) += big * phi(in,0) * phi(jn,0) * weight * v2[1];
            }//jn
        }//in
    } else {
        // outros tipos: sem contribuição especial
    }
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                                            TPZFMatrix<REAL> &ef, TPZBndCondT<STATE> &bc)
{
    TPZFMatrix<REAL> dummy;
    ContributeBC(data, weight, dummy, ef, bc);
}

// ---------- pós-processo ----------
template <class TMEM>
int TPZMatElastic2DMem<TMEM>::VariableIndex(const std::string &name) const {
    if (name == "Displacement") return EDisplacement;
    if (name == "Strain")       return EStrain;
    if (name == "Stress")       return EStress;
    if (name == "Young")        return EYoung;
    if (name == "Poisson")      return EPoisson;
    if (name == "POrder")       return EPOrder;
    return -1;
}

template <class TMEM>
int TPZMatElastic2DMem<TMEM>::NSolutionVariables(int var) const {
    switch (var) {
        case EDisplacement: return 3; // (ux,uy,0)
        case EStrain:       return 3; // (exx, eyy, gxy)
        case EStress:       return 3; // (sxx, syy, sxy)
        case EYoung:        return 1;
        case EPoisson:      return 1;
        case EPOrder:       return 1;
    }
    return 0;
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout)
{
    //const int64_t memid = data.intGlobPtIndex;
    //TMEM &Memory = this->MemItem(memid);
    // deslocamentos
    const REAL ux = (data.sol.size() && data.sol[0].size()>0) ? (REAL)data.sol[0][0] : 0.0;
    const REAL uy = (data.sol.size() && data.sol[0].size()>1) ? (REAL)data.sol[0][1] : 0.0;

    // derivadas (dudx,dudy,dvdx,dvdy)
    REAL dudx=0, dudy=0, dvdx=0, dvdy=0;
    if (data.dsol.size()){
        if (data.dsol[0].Rows()>0 && data.dsol[0].Cols()>0) dudx = data.dsol[0](0,0);
        if (data.dsol[0].Rows()>1 && data.dsol[0].Cols()>0) dudy = data.dsol[0](1,0);
        if (data.dsol[0].Rows()>0 && data.dsol[0].Cols()>1) dvdx = data.dsol[0](0,1);
        if (data.dsol[0].Rows()>1 && data.dsol[0].Cols()>1) dvdy = data.dsol[0](1,1);
    }

    // deformações de engenharia
    const REAL exx = dudx;
    const REAL eyy = dvdy;
    const REAL gxy = dudy + dvdx;

    REAL E, nu;
    ERFromMem(data, E, nu);
    //REAL E=fE_fallback;
    //REAL nu = fNu_fallback;

    TPZFNMatrix<9,REAL> D(3,3,0.0);
    BuildConstitutiveMatrix(E, nu, D);

    const REAL sxx = D(0,0)*exx + D(0,1)*eyy + D(0,2)*gxy;
    const REAL syy = D(1,0)*exx + D(1,1)*eyy + D(1,2)*gxy;
    const REAL sxy = D(2,0)*exx + D(2,1)*eyy + D(2,2)*gxy;

    switch (var) {
        case EDisplacement:
            Solout.Resize(3); Solout[0]=ux;  Solout[1]=uy;  Solout[2]=0.0; break;
        case EStrain:
            Solout.Resize(3); Solout[0]=exx; Solout[1]=eyy; Solout[2]=gxy; break;
        case EStress:
            Solout.Resize(3); Solout[0]=sxx; Solout[1]=syy; Solout[2]=sxy; break;
        case EYoung:
            Solout.Resize(1); Solout[0]=E; break;
        case EPoisson:
            Solout.Resize(1); Solout[0]=nu; break;
        case EPOrder:
            Solout.Resize(1); Solout[0]=(REAL)data.p; break;
        default:
            Solout.Resize(0);
    }
}

// ---------- persistência ----------
template <class TMEM>
int TPZMatElastic2DMem<TMEM>::ClassId() const {
    return Hash("TPZMatElastic2DMem") ^ TWithMem::ClassId() << 1;
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    buf.Write(&fE_fallback, 1);
    buf.Write(&fNu_fallback, 1);
    int p = (int)fPlane; buf.Write(&p,1);
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::Read (TPZStream &buf, void *ctx) {
    TBase::Read(buf, ctx);
    buf.Read(&fE_fallback, 1);
    buf.Read(&fNu_fallback, 1);
    int p=0; buf.Read(&p,1); fPlane = (EPlaneType)p;
}

// ---------- utils ----------
template <class TMEM>
void TPZMatElastic2DMem<TMEM>::BuildConstitutiveMatrix(REAL E, REAL nu, TPZFMatrix<REAL> &D) const
{

    D.Zero();

    // Pequeno cuidado numérico

    if (E <= 0)
    {
        std::cout << "Módulo de elasticidade nulo" <<std::endl;
        DebugStop();
    } // evita matriz nula

    // Módulo de cisalhamento (serve para ambos os casos)
    const REAL mu = E/(2.0*(1.0+nu));

    if (fPlane == EPlaneType::PlaneStress) {
        // D = E/(1-ν²) * [ [1, ν, 0],
        //                  [ν, 1, 0],
        //                  [0, 0, (1-ν)/2 * (1-ν²)/(1+ν) ] ]
        // Tradicionalmente usa-se D(2,2) = G = E/(2(1+ν))
        const REAL denom = 1.0 - nu*nu;
        const REAL c = E/denom;

        D(0,0) = c;        D(0,1) = c*nu;     D(0,2) = 0.0;
        D(1,0) = c*nu;     D(1,1) = c;        D(1,2) = 0.0;
        D(2,0) = 0.0;      D(2,1) = 0.0;      D(2,2) = mu; // = E/(2(1+ν))
    }
    else { // PlaneStrain
        // D = E/((1+ν)(1-2ν)) * [ [1-ν, ν, 0],
        //                         [ν, 1-ν, 0],
        //                         [0, 0, (1-2ν)/2] ]
        const REAL denom = (1.0+nu)*(1.0-2.0*nu);
        const REAL c = E/denom;

        D(0,0) = c*(1.0 - nu);  D(0,1) = c*nu;           D(0,2) = 0.0;
        D(1,0) = c*nu;          D(1,1) = c*(1.0 - nu);   D(1,2) = 0.0;
        D(2,0) = 0.0;           D(2,1) = 0.0;            D(2,2) = c*(1.0 - 2.0*nu)/2.0; // = mu
        // (note que esse termo é exatamente μ = E/(2(1+ν)))
    }
    //D.Print(std::cout);
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::BuildBMatrix(const TPZFMatrix<REAL> &dphix, TPZFMatrix<REAL> &B) const
{
    // dphix(dim, nshape) — já nas coordenadas globais (x,y)
    const int nshape = dphix.Cols();
    B.Zero();
    for (int a=0; a<nshape; a++){
        const REAL dNdx = dphix(0,a);
        const REAL dNdy = dphix(1,a);
        const int  iu   = 2*a;     // coluna DOF ux
        const int  iv   = 2*a + 1; // coluna DOF uy

        // exx = dudx
        B(0, iu) = dNdx;

        // eyy = dvdy
        B(1, iv) = dNdy;

        // gxy = dudy + dvdx
        B(2, iu) = dNdy; // dudy
        B(2, iv) = dNdx; // dvdx
    }
    //B.Print(std::cout);
}

template <class TMEM>
void TPZMatElastic2DMem<TMEM>::ERFromMem(const TPZMaterialDataT<STATE> &data, REAL &E, REAL &nu) const
{
    E  = fE_fallback;
    nu = fNu_fallback;

    const int64_t memid = data.intGlobPtIndex;
    if (memid < 0) return;

    const TMEM &m = this->WithMem()->MemItem(memid);
    // Ajuste os nomes dos getters se necessário
    E  = m.m_ER.E();
    nu = m.m_ER.Poisson();

}
template <class TMEM>
void TPZMatElastic2DMem<TMEM>::SetElasticResponse(const TPZElasticResponse &ER)
{
    // 1) guarda como fallback (usado se algum IP não tiver memória escrita)
    fE_fallback  = ER.E();         // ou ER.YoungModulus()
    fNu_fallback = ER.Poisson();   // ou ER.PoissonRatio()

    // 2) grava o "default" da memória (o que novos IPs recebem ao serem criados)
    TMEM memory;
    // Se TMEM tiver outros campos, inicialize-os aqui também.
    memory.m_ER.SetEngineeringData(ER.E(), ER.Poisson());
    // Em muitas branches, SetDefaultMem é herdado de TPZMatWithMem<TMEM>,
    // exatamente como no seu exemplo plástico:
    this->SetDefaultMem(memory);
    // (se sua branch expõe via WithMem(): this->WithMem()->SetDefaultMem(memory); )
}
// ---------------------- instanciação explícita padrão ----------------------
template class TPZMatElastic2DMem<TPZElasticMem>;
