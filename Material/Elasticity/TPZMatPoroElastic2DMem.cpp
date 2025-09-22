#include "TPZMatPoroElastic2DMem.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "Elasticity/TPZElasticMem.h" // ajuste o path se necessário

// ---------- ctor ----------
template <class TMEM>
TPZMatPoroElastic2DMem<TMEM>::TPZMatPoroElastic2DMem(int matid, EPlaneType plane)
: TBase(matid), fPlane(plane) {}

// ---------- memória elástica ----------
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::SetElasticResponse(const TPZElasticResponse &ER)
{
    fE_fallback  = ER.E();
    fNu_fallback = ER.Poisson();
    TMEM m;
    m.m_ER.SetEngineeringData(ER.E(), ER.Poisson());
    this->SetDefaultMem(m);
}

// ---------- requisitos (CombinedSpaces) ----------
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (auto &d : datavec){
        d.SetAllRequirements(false);
        d.fNeedsSol = true;
        d.fDeformedDirections = true;
        d.fNeedsNormal = false;
        d.fNeedsHSize = false;
    }
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::FillBoundaryConditionDataRequirements(int, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    for (auto &d : datavec){
        d.SetAllRequirements(false);
        d.fNeedsSol = true;
        d.fNeedsNormal = true;
    }
}

// ---------- pós-processo ----------
template <class TMEM>
int TPZMatPoroElastic2DMem<TMEM>::VariableIndex(const std::string &name) const {
    if (name=="Displacement") return EDisplacement;
    if (name=="Pressure")     return EPressure;
    if (name=="Flux")         return EFlux;
    if (name=="Strain")       return EStrain;
    if (name=="Stress")       return EStress;
    if (name=="Young")        return EYoung;
    if (name=="Poisson")      return EPoisson;
    if (name=="POrder")       return EPOrder;
    if (name=="ExactPressureSolution")      return EExactPressure;
    return -1;
}

template <class TMEM>
int TPZMatPoroElastic2DMem<TMEM>::NSolutionVariables(int var) const {
    switch (var) {
        case EDisplacement: return 3;
        case EPressure:     return 1;
        case EFlux:         return 3;
        case EStrain:       return 3;
        case EStress:       return 3;
        case EYoung:        return 1;
        case EPoisson:      return 1;
        case EPOrder:       return 1;
        case EExactPressure:       return 1;
    }
    return 0;
}




template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                            int var, TPZVec<REAL> &Sol)
{
    const auto &dataU = datavec[0];
    const auto &dataP = datavec[1];
    // solução exata (se fornecida)
    STATE uex = 0.;
    TPZFMatrix<STATE> duex; duex.Redim(2,1); duex.Zero();
    if (fExact) fExact(dataP.x, uex, duex);
    Sol.Resize(0);
    if (var==EYoung)   { Sol.Resize(1); Sol[0]=fE_fallback;  return; }
    if (var==EPoisson) { Sol.Resize(1); Sol[0]=fNu_fallback; return; }
    if (var==EPressure){ Sol.Resize(1); Sol[0] = (dataP.sol.size()? (REAL)dataP.sol[0][0] : 0.0); return; }
    if (var==EExactPressure){ Sol.Resize(1); Sol[0] = uex; return; }
    if (var==EPOrder)  { Sol.Resize(1); Sol[0]=(REAL)dataU.p; return; }



    // gradientes
    REAL dudx=0, dudy=0, dvdx=0, dvdy=0, dpdx=0, dpdy=0;
    if (dataU.dsol.size()){
        dudx = dataU.dsol[0](0,0);
        dudy = dataU.dsol[0](1,0);
        dvdx = dataU.dsol[0](0,1);
        dvdy = dataU.dsol[0](1,1);
    }
    if (dataP.dsol.size()){
        dpdx = dataP.dsol[0](0,0);
        dpdy = dataP.dsol[0](1,0);
    }

    if (var==EDisplacement){
        const REAL ux = dataU.sol[0][0];
        const REAL uy = dataU.sol[0][1];
        Sol.Resize(3); Sol[0]=ux; Sol[1]=uy; Sol[2]=0.0;
        return;
    }
    if (var==EStrain){
        Sol.Resize(3); Sol[0]=dudx; Sol[1]=dvdy; Sol[2]=dudy+dvdx;
        return;
    }
    if (var==EStress){
        STATE E,nu; ERFromMem(dataU,E,nu);
        TPZFNMatrix<9,STATE> D(3,3,0.0); BuildConstitutiveMatrix(E,nu,D);
        const REAL exx=dudx, eyy=dvdy, gxy=dudy+dvdx;
        const REAL sxx = D(0,0)*exx + D(0,1)*eyy + D(0,2)*gxy;
        const REAL syy = D(1,0)*exx + D(1,1)*eyy + D(1,2)*gxy;
        const REAL sxy = D(2,0)*exx + D(2,1)*eyy + D(2,2)*gxy;
        Sol.Resize(3); Sol[0]=sxx; Sol[1]=syy; Sol[2]=sxy; return;
    }
    if (var==EFlux){
        const STATE k_over_mu = fk/fmu;
        Sol.Resize(3);
        Sol[0] = -k_over_mu * (dpdx - frhof*fG[0]);
        Sol[1] = -k_over_mu * (dpdy - frhof*fG[1]);
        Sol[2] = 0.0; return;
    }
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Contribute(
    const TPZVec<TPZMaterialDataT<STATE>> &datavec,
    REAL weight,
    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (datavec.size() != 2) {
        PZError << "TPZMatPoroElastic2DMem::Contribute needs datavec.size()==2 (u,p)\n";
        DebugStop();
    }

    // ---------------- aliases ----------------
    auto &dataU  = datavec[0];
    auto &dataP  = datavec[1];
    auto &phiU   = dataU.phi;      // H1 vetorial (u)
    auto &dphiU  = dataU.dphix;    // grad N_u  (linhas = {dx,dy})
    auto &phiP   = dataP.phi;      // H1 escalar (p)
    auto &dphiP  = dataP.dphix;    // grad N_p

    const int nU   = phiU.Rows();
    const int nP   = phiP.Rows();
    const int dim  = 2;
    const int off_u = 0;
    const int off_p = dim*nU;

    // ------------- parâmetros materiais -------------
    STATE E, nu; ERFromMem(dataU, E, nu);

    TPZFMatrix<STATE> D(3,3,0.);
    BuildConstitutiveMatrix(E, nu, D);

    // Matrizes locais
    TPZFMatrix<STATE> Bu, Bp, Ke, Qe, He, Se;

    // construir “B”s
    BuildBu(dphiU, Bu);            // Bu: (3 x 2*nU) em Voigt
    BuildBp(dphiP, Bp);            // Bp: (dim x nP)

    // --- Ke = Bu^T D Bu
    TPZFMatrix<STATE> DBu; D.Multiply(Bu, DBu);
    TPZFMatrix<STATE> But; Bu.Transpose(&But);
    But.Multiply(DBu, Ke);

    // --- He = (k/μ) * (Bp^T Bp)
    TPZFMatrix<STATE> Bpt; Bp.Transpose(&Bpt);
    Bpt.Multiply(Bp, He);
    He *= (fk/fmu);
    //std::cout << "He = "<<std::endl;
    //He.Print("He");
    // --- Se = Se * (phiP^T phiP)
    TPZFMatrix<STATE> phiPt; phiP.Transpose(&phiPt);
    phiP.Multiply(phiPt, Se);      // atenção: phiP * phiP^T
    Se *= fSe;
    //std::cout << "Se = "<<std::endl;
    //Se.Print("Se");

    // --- Qe = (Bu^T m_u) (phiP^T)  com m_u=[1,1,0]^T
    TPZFMatrix<STATE> mu(3,1,0.),Qet;  mu(0,0)=1.; mu(1,0)=1.;
    TPZFMatrix<STATE> g; But.Multiply(mu, g);   // g: (2*nU x 1)
    TPZFMatrix<STATE> phiPt_loc; phiP.Transpose(&phiPt_loc);
    g.Multiply(phiPt_loc, Qe);                  // (2*nU x nP)
    Qe *= falpha;
    //std::cout << "Qe = "<<std::endl;
    //Qe.Print("Qe");
    Qe.Transpose(&Qet);
    //Qet*=1./fTimeStep;
    //Se*=1./fTimeStep;
    // =================== Montagem em ek ===================
    const int nueqs = dim * nU;   // 2*nU
    const int npeqs = nP;


    for (int i = 0; i < nueqs; ++i)
    {
        for (int j = 0; j < nueqs; ++j)
        {
            ek(off_u + i, off_u + j) += Ke(i, j)* weight;
        }

    }

    for (int i = 0; i < nueqs; ++i)
    {
        for (int j = 0; j < npeqs; ++j)
        {
            ek(off_u + i, off_p + j) += -Qe(i, j)* weight;
        }

    }




    for (int i = 0; i < npeqs; ++i)
    {
        for (int j = 0; j < nueqs; ++j)
        {
            ek(off_p + i, off_u + j) +=Qet(i, j)* weight;
        }

    }

        for (int i = 0; i < npeqs; ++i)
        {
            for (int j = 0; j < npeqs; ++j)
            {
                ek(off_p + i, off_p + j) +=fTimeStep* He(i, j)* weight;
                ek(off_p + i, off_p + j) +=Se(i, j)* weight;
            }

        }


    // =================== Vetor de carga ef ===================
    auto m = this->WithMem();

    const int gp = dataU.intGlobPtIndex;

    bool update = m->GetUpdateMem();
    if(update)UpdateMemory(datavec);

    TPZFMatrix<STATE> gvec(dim,1,0.);
    gvec(0,0) = fG[0]; gvec(1,0) = fG[1];

    // q_H = (k/μ) Bp^T ∇p0
    // q_h = ρ_f (k/μ) Bp^T g
    //TPZFMatrix<STATE> gradp(dim,1,0.); for (int k=0;k<dim;k++) gradp(k,0)=m->MemItem(gp).fdPorePressure[k];

    TPZFMatrix<REAL> gradp=dataP.dsol[0];
    TPZFMatrix<STATE> qH; Bpt.Multiply(gradp, qH); // (nP x 1)
    TPZFMatrix<STATE> qh; Bpt.Multiply(gvec, qh);
    //gradp.Print("GradP");
    //qH.Print("qH");
    //qh.Print("qh");
    //qh-qH
    for (int i = 0; i < nP; i++)
    {
        //ef(off_p + i, 0) += weight *(fk/fmu)* (qh(i,0)*frhof  - qH(i,0) )*fTimeStep ;
    }





}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight, TPZFMatrix<STATE> &ef)
{

    const auto& dataU = datavec[0];
    const auto& dataP = datavec[1];
    const int nU = dataU.phi.Rows();
    const int nP = dataP.phi.Rows();
    const int dim = 2;
    const int off_p = dim*nU;

    auto m = this->WithMem();

    bool update = m->GetUpdateMem();
    if(update)UpdateMemory(datavec);

    const int gp = dataU.intGlobPtIndex;

    TPZFMatrix<REAL> gradu=dataU.dsol[0];
    REAL pressure=dataP.sol[0][0];


    // - Q^T u^n  -> - α div(u^n) * φ_p
    const STATE divu_prev = m->MemItem(gp).fGradSolU(0,0) + m->MemItem(gp).fGradSolU(1,1);
    for (int i=0;i<nP;i++)
        ef(off_p + i,0) += weight *   falpha * (gradu(0,0)+gradu(1,1)) * dataP.phi(i,0) ;
        //ef(off_p + i,0) += weight *   falpha * divu_prev * dataP.phi(i,0) ;

    //std::cout<<"PRESSURE = " << m->MemItem(gp).fPorePressure  <<std::endl;
    // - S p^n
    for (int i=0;i<nP;i++)
        ef(off_p + i,0) += weight *  fSe * pressure * dataP.phi(i,0) ;
        //ef(off_p + i,0) += weight *  fSe * m->MemItem(gp).fPorePressure * dataP.phi(i,0) ;

    // // + (1-ξ) Δt H p^n  with H = (k/μ) ∇p · ∇φ
    // const STATE coeff = (1.0 - 1) * fTimeStep * (fk/fmu);
    // if (coeff != 0.0){
    //     for (int i=0;i<nP;i++){
    //         // grad φ_i in global
    //         STATE dphix = dataP.dphix(0,i), dphiy = dataP.dphix(1,i);
    //         STATE dot = mem.fdPorePressure[0]*dphix + mem.fdPorePressure[1]*dphiy;
    //         ef(off_p + i,0) += weight * coeff * dot;
    //     }
    // }

}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::ContributeBC(
    const TPZVec<TPZMaterialDataT<STATE>>& datavec,
    REAL weight,
    TPZFMatrix<STATE>& ek,
    TPZFMatrix<STATE>& ef,
    TPZBndCondT<STATE>& bc)
{
    if (datavec.size() != 2) {
        PZError << "ContributeBC expects datavec.size()==2 (u,p)\n";
        DebugStop();
    }


    const auto& phiU = datavec[0].phi;
    const auto& phiP = datavec[1].phi;
    const int nU = phiU.Rows();
    const int nP = phiP.Rows();

    const int off_u = 0;
    const int off_p = 2*nU;

    const STATE big = 1.e12;
    // std::cout << "[BC] mid=" << bc.Id()
    // << " type=" << bc.Type()
    // << " nU=" << datavec[0].phi.Rows()
    // << " nP=" << datavec[1].phi.Rows()
    // << " off_p=" << (2*datavec[0].phi.Rows())
    // << "\n";
    const TPZVec<STATE>& V2 = bc.Val2();

    switch (bc.Type())
    {
        // ---------------------------------------------------------
        // 0 : Dirichlet em u
        case 0:
        {

            for (int i=0;i<nU;++i){
                ef(off_u+2*i+0,0) += big * (V2[0]- datavec[0].sol[0][0]) * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += big * (V2[1]- datavec[0].sol[0][0]) * phiU(i,0) * weight;
                for (int j=0;j<nU;++j){
                    const STATE pen = big * phiU(i,0) * phiU(j,0) * weight;
                    ek(off_u+2*i+0, off_u+2*j+0) += pen;
                    ek(off_u+2*i+1, off_u+2*j+1) += pen;
                }
            }

        } break;

        // 1 : Neumann em u
        case 1:
        {
            for (int i=0;i<nU;i++){
                ef(off_u+2*i+0,0) += V2[0] * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += V2[1] * phiU(i,0) * weight;
            }

        } break;

        // 2 : Dirichlet em p
        case 2:
        {
            //V2[0] valor imposto em V2[0]
            for(int in = 0 ; in < nP; in++)
            {
                //ef(in+off_p,0)	+= (V2[0]- datavec[1].sol[0][0]) *big*phiP(in,0)*weight;	// P Pressure Value
                for (int jn = 0 ; jn < nP; jn++)
                {
                    ek(in+off_p,jn+off_p)+=big*phiP(in,0)*phiP(jn,0)*weight;	// P Pressure
                }
            }

        } break;

        // 3 : Dirichlet DIRECIONAL em u.
        case 3:
        {
            for(int in = 0 ; in < nU; in++) {
                //ef(2*in+0,0) += big * (0. - datavec[0].sol[0][0]) * V2[0] * phiU(in,0) * weight;
                //ef(2*in+1,0) += big * (0. - datavec[0].sol[0][1]) * V2[1] * phiU(in,0) * weight;
                for (int jn = 0 ; jn < nU; jn++) {
                    ek(2*in+0,2*jn+0) += big * phiU(in,0) * phiU(jn,0) * weight * V2[0];
                    ek(2*in+1,2*jn+1) += big * phiU(in,0) * phiU(jn,0) * weight * V2[1];

                }//jn
            }//in

        } break;

        default:
            break;
    }
}


template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{

}
// ---------- persistência ----------
template <class TMEM>
int TPZMatPoroElastic2DMem<TMEM>::ClassId() const {
    return Hash("TPZMatPoroElastic2DMem") ^ TBase::ClassId() << 1;
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    buf.Write(&fE_fallback,1); buf.Write(&fNu_fallback,1);
    buf.Write(&falpha,1); buf.Write(&fSe,1);
    buf.Write(&fk,1); buf.Write(&fmu,1); buf.Write(&frhof,1);
    buf.Write(fG,2); /*buf.Write(fGradP0,2)*/; buf.Write(fBody,2);
    buf.Write(&fTimeStep,1);
    int pl=(int)fPlane; buf.Write(&pl,1);
    int invdt=(int)fMassInvDt; buf.Write(&invdt,1);
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Read (TPZStream &buf, void *ctx) {
    TBase::Read(buf, ctx);
    buf.Read(&fE_fallback,1); buf.Read(&fNu_fallback,1);
    buf.Read(&falpha,1); buf.Read(&fSe,1);
    buf.Read(&fk,1); buf.Read(&fmu,1); buf.Read(&frhof,1);
    buf.Read(fG,2); /*buf.Read(fGradP0,2)*/; buf.Read(fBody,2);
    buf.Read(&fTimeStep,1);
    int pl=0; buf.Read(&pl,1); fPlane=(EPlaneType)pl;
    int invdt=1; buf.Read(&invdt,1); fMassInvDt=(bool)invdt;
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Print(std::ostream & out) const {
    out << Name() << "\n";
    out << "  plane=" << (fPlane==EPlaneType::PlaneStress?"PlaneStress":"PlaneStrain") << "\n";
    out << "  E_fallback="<<fE_fallback<<" nu_fallback="<<fNu_fallback<<"\n";
    out << "  alpha="<<falpha<<" Se="<<fSe<<" k="<<fk<<" mu="<<fmu<<" rho_f="<<frhof<<"\n";
    out << "  g=("<<fG[0]<<","<<fG[1]/*<<") gradp0=("<<fGradP0[0]<<","<<fGradP0[1]*/<<")\n";
    out << "  body=("<<fBody[0]<<","<<fBody[1]<<") dt="<<fTimeStep<<" MassInvDt="<<fMassInvDt<<"\n";
}

// ---------- helpers ----------
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::ERFromMem(const TPZMaterialDataT<STATE> &data, STATE &E, STATE &nu) const {
    E=fE_fallback; nu=fNu_fallback;
    const int64_t id = data.intGlobPtIndex;
    if (id<0) return;
    const TMEM &m = this->WithMem()->MemItem(id);
    E  = m.m_ER.E();
    nu = m.m_ER.Poisson();
}

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::BuildConstitutiveMatrix(STATE E, STATE nu, TPZFMatrix<STATE> &D) const {
    D.Redim(3,3); D.Zero();
    const STATE mu = E/(2.0*(1.0+nu));
    if (fPlane==EPlaneType::PlaneStress){
        const STATE c = E/(1.0 - nu*nu);
        D(0,0)=c; D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c; D(2,2)=mu;
    } else { // PlaneStrain
        const STATE c = E/((1.0+nu)*(1.0-2.0*nu));
        D(0,0)=c*(1.0-nu); D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c*(1.0-nu);
        D(2,2)=c*(1.0-2.0*nu)/2.0; // = mu
    }
}
// ---- Bu (mecânica) : monta matriz B_u (3 x 2*nU) a partir de dphiU(2 x nU)
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu)
{
    const int nU = dphiU.Cols();
    Bu.Redim(3, 2*nU); Bu.Zero();
    for (int a=0; a<nU; ++a) {
        const STATE dNdx = dphiU(0,a), dNdy = dphiU(1,a);
        const int iu = 2*a, iv = 2*a+1;
        Bu(0,iu) = dNdx;          // exx = dudx
        Bu(1,iv) = dNdy;          // eyy = dvdy
        Bu(2,iu) = dNdy;          // gxy = dudy
        Bu(2,iv) = dNdx;          // gxy = dvdx
    }
}

// ---- Bp (fluxo) : gradiente das funções de p agrupado (2 x nP)
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::BuildBp(const TPZFMatrix<STATE>& dphiP, TPZFMatrix<STATE>& Bp)
{
    const int nP = dphiP.Cols();
    Bp.Redim(2, nP);
    for (int j=0; j<nP; ++j) {
        Bp(0,j) = dphiP(0,j);     // dNj/dx
        Bp(1,j) = dphiP(1,j);     // dNj/dy
    }
}

// ---- Np coluna (nP x 1) útil para produtos B_u^T * (Np escalar)
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::BuildNpCol(const TPZFMatrix<STATE>& phiP, TPZFMatrix<STATE>& Np)
{
    const int nP = phiP.Rows();
    Np.Redim(nP,1);
    for (int j=0; j<nP; ++j) Np(j,0) = phiP(j,0);
}

template < class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::ComputePorePressure(const TPZMaterialDataT<STATE> & data, REAL & Pp, TPZVec<REAL> & dPp)
{
    int dim = Dimension(), i;

    const int64_t id = data.intGlobPtIndex;
    if (id<0) return;
    const TMEM &m = this->WithMem()->MemItem(id);

    Pp = m.fPorePressure;
    dPp = m.fdPorePressure;

    // adding deltaP information from n+1 time
    Pp += data.sol[0][dim];
    for(i = 0; i < dim; i++)dPp[i] += data.dsol[0](i, dim);
}
template < class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::GetPrevPorePressure(const TPZMaterialDataT<STATE>& data,STATE &P0, TPZVec<STATE> &gradP0) const {
  const int64_t id = data.intGlobPtIndex;
  if (id < 0) { P0 = 0.; gradP0.Fill(0.); return; }
  const TMEM &m = this->WithMem()->MemItem(id);
  P0 = m.fPorePressure;
  gradP0 = m.fdPorePressure;
}
template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::UpdatePorePressure(const TPZMaterialDataT<STATE> & data)
{
    int dim = Dimension(), i;
    int intPt = data.intGlobPtIndex;

    const int64_t id = data.intGlobPtIndex;
    if (id<0) return;
    TMEM &m = this->WithMem()->MemItem(id);
    // updating n+1 information
    m.fPorePressure += data.sol[0][0];

    //std::cout << data.sol << std::endl;
    //std::cout << data.dsol << std::endl;

    for(i = 0; i < dim; i++)
    {
        m.fdPorePressure[i] += data.dsol[0](i, 0);
    }

}


template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::SetMem(const REAL Pp,const TPZElasticResponse ER)
{
    // 1) guarda como fallback (usado se algum IP não tiver memória escrita)
    fE_fallback  = ER.E();         // ou ER.YoungModulus()
    fNu_fallback = ER.Poisson();   // ou ER.PoissonRatio()

    // 2) grava o "default" da memória (o que novos IPs recebem ao serem criados)
    TMEM memory;
    // Se TMEM tiver outros campos, inicialize-os aqui também.
    memory.m_ER.SetEngineeringData(ER.E(), ER.Poisson());
    memory.fPorePressure=Pp;
    // Em muitas branches, SetDefaultMem é herdado de TPZMatWithMem<TMEM>,
    // exatamente como no seu exemplo plástico:
    this->SetDefaultMem(memory);
    // (se sua branch expõe via WithMem(): this->WithMem()->SetDefaultMem(memory); )
}


// --------------- UpdateMemory -----------------
template<class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::UpdateMemory(const TPZVec<TPZMaterialDataT<STATE>>& datavec)
{
    const auto& dataU = datavec[0];
    const auto& dataP = datavec[1];
    const int gp = dataU.intGlobPtIndex;
    TMEM& mem = this->MemItem(gp);

    // store p^{n+1}
    mem.fPorePressure = dataP.sol[0][0];

    // grad p in global (assuming dphix already global; adjust if needed)
    mem.fdPorePressure[0] = dataP.dsol[0](0,0);
    mem.fdPorePressure[1] = dataP.dsol[0](1,0);

    // grad u (global)
    mem.fGradSolU(0,0) = dataU.dsol[0](0,0);
    mem.fGradSolU(0,1) = dataU.dsol[0](1,0);
    mem.fGradSolU(1,0) = dataU.dsol[0](1,0);
    mem.fGradSolU(1,1) = dataU.dsol[0](1,1);

    // optional u
    mem.fSolU[0] = dataU.sol[0][0];
    mem.fSolU[1] = dataU.sol[0][1];

    //mem.m_ER.SetEngineeringData(ER.E(), ER.Poisson());
}
// ---------- instância explícita ----------
template class TPZMatPoroElastic2DMem<TPZElasticMem>;
