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
    }
    return 0;
}




template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                            int var, TPZVec<REAL> &Sol)
{
    const auto &dataU = datavec[0];
    const auto &dataP = datavec[1];

    Sol.Resize(0);
    if (var==EYoung)   { Sol.Resize(1); Sol[0]=fE_fallback;  return; }
    if (var==EPoisson) { Sol.Resize(1); Sol[0]=fNu_fallback; return; }
    if (var==EPressure){ Sol.Resize(1); Sol[0] = (dataP.sol.size()? (REAL)dataP.sol[0][0] : 0.0); return; }
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
void TPZMatPoroElastic2DMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                              REAL weight,
                                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (datavec.size() != 2) { PZError << "Needs datavec.size()==2 (u,p)\n"; DebugStop(); }

    // --------------------------- Modo "MassOnly": monta só S(p,p) ---------------------------
    if (fMassOnly) {
        const auto &phiP = datavec[1].phi;
        const int nU = datavec[0].phi.Rows();
        const int nP = phiP.Rows();
        const int off_p = 2*nU;

        const STATE mass_coeff = fMassInvDt ? (STATE(1.0)/fTimeStep) : STATE(1.0); // S/dt ou S
        const STATE coeffS     = mass_coeff * fSe * weight;                         // (|J|*w já em weight)

        for (int i=0; i<nP; ++i) {
            const STATE Ni = phiP(i,0);
            for (int j=0; j<nP; ++j) {
                ek(off_p+i, off_p+j) += coeffS * Ni * phiP(j,0);
            }
        }
        return; // somente S e sai
    }

    // --------------------------- Dados dos dois campos (u,p) ---------------------------
    const auto &dataU = datavec[0];
    const auto &dataP = datavec[1];

    const auto &phiU  = dataU.phi;   // H1 vetorial (u)
    const auto &dphiU = dataU.dphix; // grad N_u  -> dphiU(0,a)=dNa/dx, dphiU(1,a)=dNa/dy
    const auto &phiP  = dataP.phi;   // H1 escalar (p)
    const auto &dphiP = dataP.dphix; // grad N_p

    const int nU   = phiU.Rows();
    const int nP   = phiP.Rows();
    const int off_u = 0;
    const int off_p = 2*nU;

    // --------------------------- Parâmetros materiais e de tempo ---------------------------
    STATE E, nu;
    ERFromMem(dataU, E, nu);

    TPZFNMatrix<9,STATE> D(3,3,STATE(0)); // tensor constitutivo (2D: [ex,ey,gxy])
    BuildConstitutiveMatrix(E, nu, D);

    const STATE k_over_mu  = fk/fmu;
    const STATE mass_coeff = fMassInvDt ? (STATE(1.0)/fTimeStep) : STATE(1.0); // para S e Q^T
    const STATE diff_coeff = fMassInvDt ? STATE(1.0)             : fTimeStep;  // para H e RHS_p

    // --------------------------- K(u,u): ∫ B^T D B ---------------------------
    // Sem formar B: para cada par de DOFs (i,j), montamos b_i e b_j diretamente a partir de dphiU
    // b(ux_a) = [ dNa/dx, 0,       dNa/dy ]
    // b(uy_a) = [ 0,      dNa/dy,  dNa/dx ]
    for (int a=0; a<nU; ++a) {
        const STATE dNadx = dphiU(0,a);
        const STATE dNady = dphiU(1,a);

        // Coluna i = ux_a
        const STATE bi0 = dNadx;
        const STATE bi1 = STATE(0);
        const STATE bi2 = dNady;

        // Coluna i' = uy_a
        const STATE bi0p = STATE(0);
        const STATE bi1p = dNady;
        const STATE bi2p = dNadx;

        for (int b=0; b<nU; ++b) {
            const STATE dNbdx = dphiU(0,b);
            const STATE dNbdy = dphiU(1,b);

            // Coluna j = ux_b
            const STATE bj0 = dNbdx;
            const STATE bj1 = STATE(0);
            const STATE bj2 = dNbdy;

            // Coluna j' = uy_b
            const STATE bj0p = STATE(0);
            const STATE bj1p = dNbdy;
            const STATE bj2p = dNbdx;

            // t = D * b_j
            auto Kacc_ij = [&](STATE b0, STATE b1, STATE b2,
                               int row_i) {
                const STATE t0 = D(0,0)*b0 + D(0,1)*b1 + D(0,2)*b2;
                const STATE t1 = D(1,0)*b0 + D(1,1)*b1 + D(1,2)*b2;
                const STATE t2 = D(2,0)*b0 + D(2,1)*b1 + D(2,2)*b2;
                STATE dot;
                if (row_i == 0) { // i = ux_a -> bi = [bi0,bi1,bi2]
                    dot = bi0*t0 + bi1*t1 + bi2*t2;
                } else {          // i = uy_a -> bi' = [bi0p,bi1p,bi2p]
                    dot = bi0p*t0 + bi1p*t1 + bi2p*t2;
                }
                return dot;
            };

            // (i=ux_a , j=ux_b)
            ek(off_u + 2*a + 0, off_u + 2*b + 0) += weight * Kacc_ij(bj0,  bj1,  bj2,  0);
            // (i=ux_a , j=uy_b)
            ek(off_u + 2*a + 0, off_u + 2*b + 1) += weight * Kacc_ij(bj0p, bj1p, bj2p, 0);
            // (i=uy_a , j=ux_b)
            ek(off_u + 2*a + 1, off_u + 2*b + 0) += weight * Kacc_ij(bj0,  bj1,  bj2,  1);
            // (i=uy_a , j=uy_b)
            ek(off_u + 2*a + 1, off_u + 2*b + 1) += weight * Kacc_ij(bj0p, bj1p, bj2p, 1);
        }
    }

    // --------------------------- f_u: força de corpo (u) ---------------------------
    for (int a=0; a<nU; ++a) {
        const STATE Na = phiU(a,0);
        ef(off_u + 2*a + 0, 0) += weight * Na * fBody[0];
        ef(off_u + 2*a + 1, 0) += weight * Na * fBody[1];
    }

    // --------------------------- -Q(u,p) e +Q^T(p,u) ---------------------------
    // -Q(u,p)   = -α ∫ Np * [dNa/dx ; dNa/dy]
    // +Q^T(p,u) = +α/Δt (incremental) ou +α (direto)
    const STATE cQ  = (-falpha) * weight;                                   // -Q
    const STATE cQT = (fMassInvDt ? (falpha/fTimeStep) : falpha) * weight;  // +Q^T

    for (int a=0; a<nU; ++a) {
        const STATE dNadx = dphiU(0,a);
        const STATE dNady = dphiU(1,a);

        for (int j=0; j<nP; ++j) {
            const STATE Nj = phiP(j,0);

            const STATE qx = Nj * dNadx;  // parte x
            const STATE qy = Nj * dNady;  // parte y

            // -Q  (linha u, coluna p)
            ek(off_u + 2*a + 0, off_p + j) += cQ  * qx;
            ek(off_u + 2*a + 1, off_p + j) += cQ  * qy;

            // +Q^T (linha p, coluna u)
            ek(off_p + j,       off_u + 2*a + 0) += cQT * qx;
            ek(off_p + j,       off_u + 2*a + 1) += cQT * qy;
        }
    }

    // --------------------------- H(p,p): difusão (k/μ) ∫ ∇N·∇N ---------------------------
    {
        const STATE factorH = diff_coeff * k_over_mu * weight; // (Δt)*(k/μ) no direto; (1)*(k/μ) no incremental
        for (int i=0; i<nP; ++i) {
            const STATE dNix = dphiP(0,i);
            const STATE dNiy = dphiP(1,i);
            for (int j=0; j<nP; ++j) {
                const STATE dNjx = dphiP(0,j);
                const STATE dNjy = dphiP(1,j);
                ek(off_p + i, off_p + j) += factorH * (dNix*dNjx + dNiy*dNjy);
            }
        }
    }

    // --------------------------- S(p,p): armazenamento S_e ∫ N N ---------------------------
    {
        const STATE coeffS = mass_coeff * fSe * weight; // S/dt (incremental) ou S (direto)
        for (int i=0; i<nP; ++i) {
            const STATE Ni = phiP(i,0);
            for (int j=0; j<nP; ++j) {
                ek(off_p + i, off_p + j) += coeffS * Ni * phiP(j,0);
            }
        }
    }

    // --------------------------- RHS de p: fontes volumétricas s(x) (opcional) ------------
    if (fForcingP) {
        TPZManVector<STATE> res(1,0.0);
        for (int i=0; i<nP; ++i) {
            fForcingP->Execute(dataP.x, res); // s(x) no ponto
            ef(off_p + i, 0) += diff_coeff * weight * phiP(i,0) * res[0];
        }
    }

    // Observação: gravidade do fluxo (ρ_f g) entra via BC de fluxo em ContributeBC; não há termo
    // volumétrico quando k, μ e g são constantes (∇·(k/μ ρ_f g)=0).
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

#ifndef gBigNumber
    const STATE big = 1.e12;
#else
    const STATE big = gBigNumber;
#endif

    // fator de tempo p/ termos de FLUXO (Neumann/Robin) de pressão
    //const STATE flux_scale = (this->fMassInvDt ? STATE(1.0) : this->fTimeStep);
    const STATE diff_coeff = fMassInvDt ? 1.0             : fTimeStep;
    // ---- leitura dos valores de contorno ----
    auto read_vals = [&](STATE& tx, STATE& ty, STATE& pD)
    {
        tx = ty = pD = 0.;
        if (bc.HasForcingFunctionBC()) {
            TPZManVector<STATE> res(3,0.0);
            //bc.ForcingFunctionBC()->Execute(datavec[1].x, res); // res[0]=tx, res[1]=ty, res[2]=pD
            if (res.size() > 0) tx = res[0];
            if (res.size() > 1) ty = res[1];
            if (res.size() > 2) pD = res[2];
        } else {
            const TPZVec<STATE>& V2 = bc.Val2(); // vetor
            const int nv = V2.NElements();
            if (nv > 0) tx = V2[0];
            if (nv > 1) ty = V2[1];
            if (nv > 2) pD = V2[2];
        }
    };

    switch (bc.Type())
    {
        // ---------------------------------------------------------
        // 0 : Dirichlet completo em u e p  (penalização)
        case 0:
        {
            STATE ux, uy, pD; read_vals(ux,uy,pD);

            // u = (ux,uy) por penalização
            for (int i=0;i<nU;++i){
                ef(off_u+2*i+0,0) += big * ux * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += big * uy * phiU(i,0) * weight;
                for (int j=0;j<nU;++j){
                    const STATE pen = big * phiU(i,0) * phiU(j,0) * weight;
                    ek(off_u+2*i+0, off_u+2*j+0) += pen;
                    ek(off_u+2*i+1, off_u+2*j+1) += pen;
                }
            }
            // p = pD por penalização (SEM Δt)
            for (int i=0;i<nP;++i){
                ef(off_p+i,0) += big * pD * phiP(i,0) * weight;
                for (int j=0;j<nP;++j){
                    ek(off_p+i, off_p+j) += big * phiP(i,0) * phiP(j,0) * weight;
                }
            }
        } break;

        // ---------------------------------------------------------
        // 1 : Neumann em u (tração), nada em p
        case 1:
        {
            STATE tx, ty, dummy; read_vals(tx,ty,dummy);
            for (int i=0;i<nU;++i){
                ef(off_u+2*i+0,0) += tx * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += ty * phiU(i,0) * weight;
            }
        } break;

        // ---------------------------------------------------------
        // 2 : Dirichlet apenas em p  (penalização)
        case 2:
        {
            STATE dummyx, dummyy, pD; read_vals(dummyx,dummyy,pD);
            for (int i=0;i<nP;++i){
                ef(off_p+i,0) += big * pD * phiP(i,0) * weight;
                for (int j=0;j<nP;++j){
                    ek(off_p+i, off_p+j) += big * phiP(i,0) * phiP(j,0) * weight;
                }
            }
        } break;

        // ---------------------------------------------------------
        // 3 : Neumann em p (fluxo normal q_n) -> só RHS, com escala temporal
        case 3:
        {
            STATE tx, ty, qn; read_vals(tx,ty,qn);
            for (int i=0;i<nP;++i){
                //ef(off_p+i,0) += flux_scale * qn * phiP(i,0) * weight;
                ef(off_p+i,0) += diff_coeff * weight * qn * phiP(i,0);
            }
        } break;

        // ---------------------------------------------------------
        // 4 : Dirichlet DIRECIONAL em u (máscara [mask_x,mask_y] em Val2);
        //     nada em p.
        case 4:
        {
            for (int i=0;i<nU;++i){
                //ef(off_u+2*i+0,0) += big * mask_x * ux * phiU(i,0) * weight;
                //ef(off_u+2*i+1,0) += big * mask_y * uy * phiU(i,0) * weight;
                for (int j=0;j<nU;++j){
                    ek(off_u+2*i+0, off_u+2*j+0) +=big * bc.Val2()[0] * phiU(i,0) * phiU(j,0) * weight;
                    ek(off_u+2*i+1, off_u+2*j+1) +=big * bc.Val2()[1] * phiU(i,0) * phiU(j,0) * weight;
                }
            }
            // for(in = 0 ; in < phr; in++) {
            //     //                ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
            //     //                ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
            //     const auto &v2 = bc.Val2();
            //     for (jn = 0 ; jn < phr; jn++) {
            //         ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
            //         ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
            //     }//jn
            // }//in
        } break;

        // ---------------------------------------------------------
        // 5 : Robin em p  -> β(p - pD) ~ q_n  (ek += β*NiNj , ef += β*pD*Ni)
        case 5:
        {
            const STATE beta = (bc.Val1().Rows()>2 ? bc.Val1()(2,0) : 0.);
            const TPZVec<STATE>& V2 = bc.Val2();
            const STATE pD = (V2.NElements()>2 ? V2[2] : 0.);
            const STATE coeff = diff_coeff * beta;
            for (int i=0;i<nP;++i){
                ef(off_p+i,0) += coeff * pD * phiP(i,0) * weight;
                for (int j=0;j<nP;++j){
                    ek(off_p+i, off_p+j) += coeff * phiP(i,0) * phiP(j,0) * weight;
                }
            }
        } break;

        // ---------------------------------------------------------
        // 10 : Neumann em u  +  Dirichlet em p (atalho)
        case 10:
        {
            STATE tx, ty, pD; read_vals(tx,ty,pD);
            // Neumann-u
            for (int i=0;i<nU;++i){
                ef(off_u+2*i+0,0) += tx * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += ty * phiU(i,0) * weight;
            }
            // Dirichlet-p
            for (int i=0;i<nP;++i){
                ef(off_p+i,0) += big * pD * phiP(i,0) * weight;
                for (int j=0;j<nP;++j){
                    ek(off_p+i, off_p+j) += big * phiP(i,0) * phiP(j,0) * weight;
                }
            }
        } break;

        default:
            break;
    }
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

template <class TMEM>
void TPZMatPoroElastic2DMem<TMEM>::BuildBMatrix(const TPZFMatrix<STATE> &dphix, TPZFMatrix<STATE> &B) const {
    const int n = dphix.Cols();
    B.Redim(3, 2*n); B.Zero();
    for (int a=0;a<n;a++){
        const STATE dNdx = dphix(0,a);
        const STATE dNdy = dphix(1,a);
        const int iu = 2*a, iv = 2*a+1;
        B(0,iu) = dNdx;   // exx = dudx
        B(1,iv) = dNdy;   // eyy = dvdy
        B(2,iu) = dNdy;   // gxy = dudy
        B(2,iv) = dNdx;   // gxy = dvdx
    }
}

// ---------- instância explícita ----------
template class TPZMatPoroElastic2DMem<TPZElasticMem>;
