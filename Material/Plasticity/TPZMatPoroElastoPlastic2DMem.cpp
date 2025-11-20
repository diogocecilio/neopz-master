#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "Elasticity/TPZElasticMem.h"
#include "TPZMatPoroElastoPlastic2DMem.h"


template<class T, class TMEM>
TPZMatPoroElastoPlastic2DMem<T,TMEM>::TPZMatPoroElastoPlastic2DMem() : TBase()
{
    fPlaneStrain = 1;      // ou 0, conforme o seu padrão
    fFactor = 1.0;
    //fRhoBulkDensity = 0.0;
    fBody[0]=0.;
    fBody[1]=0.;
}


template <class T, class TMEM>
TPZMatPoroElastoPlastic2DMem<T,TMEM>::TPZMatPoroElastoPlastic2DMem(int id): TBase(id)
{
    fPlaneStrain = 1;      // ou 0, conforme o seu padrão
    fFactor = 1.0;
    //fRhoBulkDensity = 0.0;
    fBody[0]=0.;
    fBody[1]=0.;
}


template <class T, class TMEM>
TPZMatPoroElastoPlastic2DMem<T,TMEM>::TPZMatPoroElastoPlastic2DMem(const TPZMatPoroElastoPlastic2DMem &cp)
: fPlaneStrain(cp.fPlaneStrain), fFactor(cp.fFactor),
/*fRhoBulkDensity(cp.fRhoBulkDensity),*/fPlasticityModel(cp.fPlasticityModel), m_PER(cp.m_PER)
{
    fBody[0]=cp.fBody[0];
    fBody[1]=cp.fBody[1];
}

template <class T, class TMEM>
TPZMatPoroElastoPlastic2DMem<T,TMEM>::~TPZMatPoroElastoPlastic2DMem()
{

}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Print(std::ostream & out) const {
    out << Name() << "\n";
    out << "\n with template argurment T = " << fPlasticityModel.Name();
    out << "\n Base material Data:\n";
    WithMem()->PrintMem(out);
    out << "\n Localy defined members:";
    out << "\n Body Forces: " << fBody[0] << " ,  "<< fBody[1];

    out << "\n Internal plasticity <T> member:\n";
    fPlasticityModel.Print(out);
}

template <class T, class TMEM>
int TPZMatPoroElastoPlastic2DMem<T,TMEM>::VariableIndex(const std::string &name) const
{
    if (name=="Displacement") return EDisplacement;
    if (name=="Pressure")     return EPressure;
    if (name=="PlasticStrain")     return EPlasticStrain;
    if (name=="ElasticStrain")     return ElasticStrain;
    if (name=="StrainPlasticJ2")     return EStrainPlasticJ2;
    if (name=="StrainPlasticI1")     return EStrainPlasticI1;
    if (name=="Coesion")     return ECoesion;
    if (name=="Atrito")     return EAtrito;

    return 0;
}

template <class T, class TMEM>
int TPZMatPoroElastoPlastic2DMem<T,TMEM>::NSolutionVariables(int var) const
{
    switch (var) {
        case EDisplacement: return 3;
        case EPlasticStrain: return 3;
        case ElasticStrain:     return 3;
        case EPressure:     return 1;
        case EStrainPlasticJ2: return 1;
        case EStrainPlasticI1:     return 1;
        case ECoesion: return 1;
        case EAtrito:     return 1;

    }
    return 0;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,int var, TPZVec<REAL> &Sol)
{
    const auto &dataU = datavec[0];
    const auto &dataP = datavec[1];

    auto m = this->WithMem();
    const int gp = dataU.intGlobPtIndex;


    TMEM &Memory = m->MemItem(gp);


    TPZTensor<REAL> & totalStrain = Memory.m_elastoplastic_state.m_eps_t;
    TPZTensor<REAL> & plasticStrain = Memory.m_elastoplastic_state.m_eps_p;
    TPZTensor<REAL> & elasticStrain = totalStrain;
    elasticStrain-=plasticStrain;
    int sz= Memory.m_elastoplastic_state.fmatprop.size();
    if ( sz==0 ) {
        Memory.m_elastoplastic_state.fmatprop.Resize ( 3 );
        Memory.m_elastoplastic_state.fmatprop[0]=0.;
        Memory.m_elastoplastic_state.fmatprop[1]=0.;
        Memory.m_elastoplastic_state.fmatprop[2]=0.;
    }
    int sz2= Memory.m_elastoplastic_state.fflux.size();
    if ( sz2==0 ) {
        Memory.m_elastoplastic_state.fflux.Resize ( 3 );
        Memory.m_elastoplastic_state.fflux[0]=0.;
        Memory.m_elastoplastic_state.fflux[1]=0.;
        Memory.m_elastoplastic_state.fflux[2]=0.;
    }
    if (var==EDisplacement){
        const REAL ux = dataU.sol[0][0];
        const REAL uy = dataU.sol[0][1];
        Sol.Resize(3); Sol[0]=ux; Sol[1]=uy; Sol[2]=0.0;
        return;
    }
    if (var==EPressure)
    {
        Sol.Resize(1);
        Sol[0] = dataP.sol[0][0];
        return;

    }
    if (var==ElasticStrain)
    {
        Sol.Resize(3);
        Sol[0] = elasticStrain.XX();
        Sol[1] = elasticStrain.YY();
        Sol[2] = elasticStrain.XY();
        return;

    }
    if (var==EPlasticStrain)
    {
        Sol.Resize(3);
        Sol[0] = plasticStrain.XX();
        Sol[1] = plasticStrain.YY();
        Sol[2] = plasticStrain.XY();
        return;

    }
    if (var==EStrainPlasticJ2)
    {
        Sol.Resize(1);
         Sol[0] = plasticStrain.J2();
        return;

    }
    if (var==EStrainPlasticI1)
    {
        Sol.Resize(1);
         Sol[0] = plasticStrain.I1();
        return;

    }
    if (var==ECoesion)
    {
        Sol.Resize(1);
         Sol[0] = Memory.m_elastoplastic_state.fmatprop[0];
        return;

    }
    if (var==EAtrito)
    {
        Sol.Resize(1);
         Sol[0] = Memory.m_elastoplastic_state.fmatprop[1];
        return;

    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    const int n = datavec.size();
    for(int i=0;i<n;i++){
        datavec[i].SetAllRequirements(true);
        // datavec[i].fNeedsSol = true;
        // datavec[i].fNeedsNormal = false;
        // datavec[i].fNeedsHSize = false;
        // datavec[i].fNeedsNeighborCenter = false;
    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::FillBoundaryConditionDataRequirements(int, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    //DebugStop();
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{

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

    TPZFNMatrix<9>  Deriv(3,3);
    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);

    //dataU.Print(std::cout);
    this->ComputeDeltaStrainVector(dataU, DeltaStrain);
    this->ApplyDeltaStrainComputeDep(dataU, DeltaStrain, Stress, Dep);

    TPZFMatrix<STATE> D(3,3,0.),Sigma2D(3,1,0.);
    if(fConsistentMatrix)
    {
        D(0, 0)=Dep(_XX_, _XX_);D(0, 1)=Dep(_XX_, _YY_);D(0,2)=Dep(_XX_, _XY_);
        D(1, 0)=Dep(_YY_, _XX_);D(1, 1)=Dep(_YY_, _YY_);D(1,2)=Dep(_YY_, _XY_);
        D(2, 0)=Dep(_XY_, _XX_);D(2, 1)=Dep(_XY_, _YY_);D(2,2)=Dep(_XY_, _XY_);
    }else{
        const int64_t id = dataU.intGlobPtIndex;
        const TMEM &m = this->WithMem()->MemItem(id);
        STATE E  = m.m_ER.E();
        STATE nu = m.m_ER.Poisson();
        const STATE c = E/((1.0+nu)*(1.0-2.0*nu));
        D(0,0)=c*(1.0-nu); D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c*(1.0-nu);
        D(2,2)=c*(1.0-2.0*nu)/2.0; // = mu
    }



    TPZFMatrix<STATE> Bu,But, Bp, Ke, Qe, He, Se,temp,Nu,Nut,fint1,fint2;

    //Matriz de rigidez
    BuildBu(dphiU, Bu);            // Bu: (3 x 2*nU) em Voigt
    BuildBp(dphiP, Bp);            // Bp: (dim x nP)

    Bu.Transpose(&But);
    But.Multiply(D,temp);
    temp.Multiply(Bu,Ke);

    //Stress.Print("Stress");
    //Forças internas
    Sigma2D(0,0)=Stress(_XX_, 0);
    Sigma2D(1,0)=Stress(_YY_, 0);
    Sigma2D(2,0)=Stress(_XY_, 0);
    //Sigma2D*=0.;
    //Sigma2D.Print("Sigma2D");

    BuildNu(phiU,Nu);
    Nu.Transpose(&Nut);

    But.Multiply(Sigma2D,fint1);

    //Forças de volume
    TPZFMatrix<STATE> gvec(dim,1,0.);
    gvec(0,0) = fBody[0]; gvec(1,0) = fBody[1];

    Nut.Multiply(gvec,fint2);

    TPZFMatrix<STATE> Bpt; Bp.Transpose(&Bpt);
    Bpt.Multiply(Bp, He);
    He *= (fk/fmu);
    //std::cout << "He = "<<std::endl;
    //He.Print("He");
    // --- Se = Se * (phiP^T phiP)
    TPZFMatrix<STATE> phiPt; phiP.Transpose(&phiPt);
    phiP.Multiply(phiPt, Se);      // ateno: phiP * phiP^T
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



    // =================== Montagem em ek ===================
    const int nueqs = dim * nU;   // 2*nU
    const int npeqs = nP;

    // TPZFMatrix<REAL> axesT, dphiXYZ;
    // // rotating the shape functions to the XYZ coordinates
    // const TPZFMatrix<REAL> &axes = dataU.axes;
    // axes.Transpose(&axesT);
    // axesT.Multiply(dphiU,dphiXYZ);
    // for(int i=0;i<nU;i++)
    // {
    //
    //     STATE val=0.;
    //     val  =  fBody[0] * phiU(i,0); // fb
    //     val -= Stress(_XX_,0) * dphiXYZ(0,i); // |
    //     val -= Stress(_XY_,0) * dphiXYZ(1,i);
    //     ef(2*i+0,0) += weight * val;
    //
    //     //Second equation: fb and fk
    //     val  = fBody[1]  * phiU(i,0); // fb
    //     val -= Stress(_XY_,0) * dphiXYZ(0,i); // |
    //     val -= Stress(_YY_,0) * dphiXYZ(1,i); // fk
    //
    //     ef(2*i+1,0) += weight * val;
    // }


    //euler implicito
    // |K  -Q      | |un1|= |fu           |
    // |QT  S+dt H | |pn1|  |fq+QT un+S pn|

    switch (fWhichassemble)
    {
        case EK:
        {
            for (int i = 0; i < nueqs; ++i)
            {

                ef(off_u + i,0) += weight * (fint2(i,0)-fint1(i,0));
                for (int j = 0; j < nueqs; ++j)
                {
                    ek(off_u + i, off_u + j) += Ke(i, j)* weight;
                }

        }
        break;
        case EQ:
        {

            for (int i = 0; i < nueqs; ++i)
            {
                for (int j = 0; j < npeqs; ++j)
                {
                    ek(off_u + i, off_p + j) += -Qe(i, j)* weight;
                }

            }
        }
        break;
        case EQT:
        {
            for (int i = 0; i < npeqs; ++i)
            {
                for (int j = 0; j < nueqs; ++j)
                {
                    ek(off_p + i, off_u + j) +=Qet(i, j)* weight;
                }

            }
        }
        break;
        case ES:
        {
            for (int i = 0; i < npeqs; ++i)
            {
                for (int j = 0; j < npeqs; ++j)
                {
                    ek(off_p + i, off_p + j) +=Se(i, j)* weight;
                }
            }
        }
        break;
        case EH:
        {
            for (int i = 0; i < npeqs; ++i)
            {
                for (int j = 0; j < npeqs; ++j)
                {
                    ek(off_p + i, off_p + j) +=fTimeStep* He(i, j)* weight;

                }
            }
        }
        break;

        default:
            break;
    }
    }


    TPZTensor<STATE> gradu;
    TPZFMatrix<STATE> gradp(dim,1,0.);

    auto m = this->WithMem();

    const int gp = dataU.intGlobPtIndex;


    //gradp=dataP.dsol[0];
    for (int k=0;k<dim;k++) gradp(k,0)=m->MemItem(gp).m_elastoplastic_state.fdPorePressure[k];

    //gradu=dataU.dsol[0];
    gradu=m->MemItem(gp).m_elastoplastic_state.fGradSolU;

    //REAL pressure=dataP.sol[0][0];
    REAL pressure=m->MemItem(gp).m_elastoplastic_state.fpressure;


    TPZFMatrix<STATE> qH;
    Bpt.Multiply(gradp, qH); // (nP x 1) (k/mu) Bp^T gradp0

    TPZFMatrix<STATE> qh;
    Bpt.Multiply(gvec, qh); //rhof_f (k/mu) Bp^T g


    //std::cout << "pressure = "<<pressure <<std::endl;
    for (int i = 0; i < npeqs; ++i)
    {
        ef(off_p + i,0) += weight *(fk/fmu)* (  - qH(i,0) )*fTimeStep ;
        //ef(off_p + i,0) += weight *(fk/fmu)* (qh(i,0)*frhof  - qH(i,0) )*fTimeStep ;
        // ef(off_p + i,0) += weight *   falpha * (gradu.I1()) * dataP.phi(i,0) ;
        // ef(off_p + i,0) += weight *  fSe * pressure * dataP.phi(i,0) ;//- S p^n

    }

    #ifdef PZ_LOG
    {
        std::ostringstream oss;
            const int dim   = 2;
            const int nU    = phiU.Rows();
            const int nP    = phiP.Rows();
            const int off_u = 0;
            const int off_p = dim*nU;

            oss.setf(std::ios::scientific);
            oss.precision(6);
/*
            oss << "\n[KQHS] gp=" << dataU.intGlobPtIndex
            << "  K("  << Ke.Rows()  << "x" << Ke.Cols()  << ")"
            << "  Q("  << Qe.Rows()  << "x" << Qe.Cols()  << ")"
            << "  QT(" << Qet.Rows() << "x" << Qet.Cols() << ")"
            << "  H("  << He.Rows()  << "x" << He.Cols()  << ")"
            << "  S("  << Se.Rows()  << "x" << Se.Cols()  << ")\n";*/
/*
            // Amostras completas (PZ): cuidado, podem ser grandes!
            Ke.Print("Ke", oss);
            Qe.Print("Qe", oss);
            Qet.Print("Qet", oss);
            He.Print("He", oss);
            Se.Print("Se", oss);*/

            // // Verificar Qet == (Qe)^T
            // TPZFMatrix<STATE> QeT; Qe.Transpose(&QeT);
            // QeT -= Qet;
            // QeT.Print("Qe^T - Qet (deve ser ~0)", oss);

            // Resíduos elementares (L1 simples) para u e p
            const int npeqs = nP;
            TPZFMatrix<STATE> fu1(nueqs,1,0.), fp1(npeqs,1,0.);
            for (int i=0;i<nueqs;++i) fu1(i,0) = ef(off_u+i,0);
            for (int i=0;i<nP;++i)    fp1(i,0) = ef(off_p+i,0);

            oss << "Conribute(ek,ef) "
            << "|fu|_1=" << Norm(fu1) << "  |fp|_1=" << Norm(fp1)
            << "|Bt Sigma|=" << Norm(fint1)*weight << "  |Nt b|=" << Norm(fint2)*weight
            << "  dt=" << fTimeStep
            << "  alpha=" << falpha
            << "  k/mu=" << (fk/fmu)
            << "  Se=" << fSe
            << "\n";

            LOGPZ_DEBUG(logger_poro, oss.str());





    }
    #endif // PZ_PORO_AUDIT

}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                                      REAL weight, TPZFMatrix<STATE> &ef)
{

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

    TPZFNMatrix<9>  Deriv(3,3);
    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);

    //dataU.Print(std::cout);
    this->ComputeDeltaStrainVector(dataU, DeltaStrain);
    this->ApplyDeltaStrainComputeDep(dataU, DeltaStrain, Stress, Dep);

    TPZFMatrix<STATE> D(3,3,0.),Sigma2D(3,1,0.);
    if(fConsistentMatrix)
    {
        D(0, 0)=Dep(_XX_, _XX_);D(0, 1)=Dep(_XX_, _YY_);D(0,2)=Dep(_XX_, _XY_);
        D(1, 0)=Dep(_YY_, _XX_);D(1, 1)=Dep(_YY_, _YY_);D(1,2)=Dep(_YY_, _XY_);
        D(2, 0)=Dep(_XY_, _XX_);D(2, 1)=Dep(_XY_, _YY_);D(2,2)=Dep(_XY_, _XY_);
    }else{
        const int64_t id = dataU.intGlobPtIndex;
        const TMEM &m = this->WithMem()->MemItem(id);
        STATE E  = m.m_ER.E();
        STATE nu = m.m_ER.Poisson();
        const STATE c = E/((1.0+nu)*(1.0-2.0*nu));
        D(0,0)=c*(1.0-nu); D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c*(1.0-nu);
        D(2,2)=c*(1.0-2.0*nu)/2.0; // = mu
    }


    TPZFMatrix<STATE> Bu,But, Bp, Ke, Qe, He, Se,temp,Nu,Nut,fint1,fint2;

    //Matriz de rigidez
    BuildBu(dphiU, Bu);            // Bu: (3 x 2*nU) em Voigt
    BuildBp(dphiP, Bp);            // Bp: (dim x nP)

    Bu.Transpose(&But);
    But.Multiply(D,temp);
    temp.Multiply(Bu,Ke);

    //Stress.Print("Stress");
    //Forças internas
    Sigma2D(0,0)=Stress(_XX_, 0);
    Sigma2D(1,0)=Stress(_YY_, 0);
    Sigma2D(2,0)=Stress(_XY_, 0);
    //Sigma2D*=0.;
    //Sigma2D.Print("Sigma2D");

    BuildNu(phiU,Nu);
    Nu.Transpose(&Nut);

    But.Multiply(Sigma2D,fint1);

    //Forças de volume
    TPZFMatrix<STATE> gvec(dim,1,0.);
    gvec(0,0) = fBody[0]; gvec(1,0) = fBody[1];

    Nut.Multiply(gvec,fint2);

    TPZFMatrix<STATE> Bpt;
    Bp.Transpose(&Bpt);

    const int nueqs = dim * nU;   // 2*nU
    const int npeqs = nP;


    for (int i = 0; i < nueqs; ++i)
    {
        ef(off_u + i,0) += weight * (fint2(i,0)-fint1(i,0));
    }

    TPZTensor<STATE> gradu;
    TPZFMatrix<STATE> gradp(dim,1,0.);

    auto m = this->WithMem();

    const int gp = dataU.intGlobPtIndex;


    //gradp=dataP.dsol[0];
    for (int k=0;k<dim;k++) gradp(k,0)=m->MemItem(gp).m_elastoplastic_state.fdPorePressure[k];

    //gradu=dataU.dsol[0];
    gradu=m->MemItem(gp).m_elastoplastic_state.fGradSolU;

    //REAL pressure=dataP.sol[0][0];
    REAL pressure=m->MemItem(gp).m_elastoplastic_state.fpressure;


    TPZFMatrix<STATE> qH;
    Bpt.Multiply(gradp, qH); // (nP x 1) (k/mu) Bp^T gradp0

    TPZFMatrix<STATE> qh;
    Bpt.Multiply(gvec, qh); //rhof_f (k/mu) Bp^T g


    //std::cout << "pressure = "<<pressure <<std::endl;
    for (int i = 0; i < npeqs; ++i)
    {
        ef(off_p + i,0) += weight *(fk/fmu)* (  - qH(i,0) )*fTimeStep ;
        //ef(off_p + i,0) += weight *(fk/fmu)* (qh(i,0)*frhof  - qH(i,0) )*fTimeStep ;
        ef(off_p + i,0) += weight *   falpha * (gradu.I1()) * dataP.phi(i,0) ;
        ef(off_p + i,0) += weight *  fSe * pressure * dataP.phi(i,0) ;//- S p^n

    }

    #ifdef PZ_LOG
    {
        std::ostringstream oss;
        const int dim   = 2;
        const int nU    = phiU.Rows();
        const int nP    = phiP.Rows();
        const int off_u = 0;
        const int off_p = dim*nU;

        oss.setf(std::ios::scientific);
        oss.precision(6);

        // Resíduos elementares (L1 simples) para u e p
        const int npeqs = nP;
        TPZFMatrix<STATE> fu1(nueqs,1,0.), fp1(npeqs,1,0.);
        for (int i=0;i<nueqs;++i) fu1(i,0) = ef(off_u+i,0);
        for (int i=0;i<nP;++i)    fp1(i,0) = ef(off_p+i,0);

        oss << "Conribute(ek,ef) "
        << "|fu|_1=" << Norm(fu1) << "  |fp|_1=" << Norm(fp1)
        << "|Bt Sigma|=" << Norm(fint1)*weight << "  |Nt b|=" << Norm(fint2)*weight
        << "  dt=" << fTimeStep
        << "  alpha=" << falpha
        << "  k/mu=" << (fk/fmu)
        << "  Se=" << fSe
        << "\n";

        LOGPZ_DEBUG(logger_poro, oss.str());


    }
    #endif // PZ_PORO_AUDIT
}
// template <class T, class TMEM>
// void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
//                                                       REAL weight, TPZFMatrix<STATE> &ef)
// {
//     auto &dataU  = datavec[0];
//     auto &dataP  = datavec[1];
//     auto &phiU   = dataU.phi;      // H1 vetorial (u)
//     auto &dphiU  = dataU.dphix;    // grad N_u  (linhas = {dx,dy})
//     auto &phiP   = dataP.phi;      // H1 escalar (p)
//     auto &dphiP  = dataP.dphix;    // grad N_p
//
//     const int nU   = phiU.Rows();
//     const int nP   = phiP.Rows();
//     const int dim  = 2;
//     const int off_u = 0;
//     const int off_p = dim*nU;
//
//     TPZFNMatrix<9>  Deriv(3,3);
//     TPZFNMatrix<36> Dep(6,6);
//     TPZFNMatrix<6>  DeltaStrain(6,1);
//     TPZFNMatrix<6>  Stress(6,1);
//
//     // --- parte mecânica (mesma do ek,ef) ---
//     this->ComputeDeltaStrainVector(dataU, DeltaStrain);
//     this->ApplyDeltaStrainComputeDep(dataU, DeltaStrain, Stress, Dep);
//
//     TPZFMatrix<STATE> D(3,3,0.), Sigma2D(3,1,0.);
//     if(fConsistentMatrix)
//     {
//         D(0,0)=Dep(_XX_,_XX_); D(0,1)=Dep(_XX_,_YY_); D(0,2)=Dep(_XX_,_XY_);
//         D(1,0)=Dep(_YY_,_XX_); D(1,1)=Dep(_YY_,_YY_); D(1,2)=Dep(_YY_,_XY_);
//         D(2,0)=Dep(_XY_,_XX_); D(2,1)=Dep(_XY_,_YY_); D(2,2)=Dep(_XY_,_XY_);
//     }
//     else
//     {
//         const int64_t id = dataU.intGlobPtIndex;
//         const TMEM &m = this->WithMem()->MemItem(id);
//         STATE E  = m.m_ER.E();
//         STATE nu = m.m_ER.Poisson();
//         const STATE c = E/((1.0+nu)*(1.0-2.0*nu));
//         D(0,0)=c*(1.0-nu); D(0,1)=c*nu; D(1,0)=c*nu; D(1,1)=c*(1.0-nu);
//         D(2,2)=c*(1.0-2.0*nu)/2.0;
//     }
//
//     TPZFMatrix<STATE> Bu, But, Nu, Nut, fint1, fint2,temp,Ke;
//     BuildBu(dphiU, Bu);
//
//
//     Bu.Transpose(&But);
//     But.Multiply(D,temp);
//     temp.Multiply(Bu,Ke);
//
//     Bu.Transpose(&But);
//
//
//     // forças internas (mesmo sinal do ek,ef)
//     Sigma2D(0,0) = Stress(_XX_,0);
//     Sigma2D(1,0) = Stress(_YY_,0);
//     Sigma2D(2,0) = Stress(_XY_,0);
//
//     BuildNu(phiU, Nu);
//     Nu.Transpose(&Nut);
//
//     But.Multiply(Sigma2D, fint1);
//
//     // forças de volume mecânicas
//     TPZFMatrix<STATE> gvec(dim,1,0.);
//     gvec(0,0) = fBody[0];
//     gvec(1,0) = fBody[1];
//     Nut.Multiply(gvec, fint2);
//
//     // monta Ru
//     const int nueqs = dim*nU;
//     for (int i=0; i<nueqs; ++i)
//     {
//         ef(off_u + i,0) += weight * (fint2(i,0)-fint1(i,0));
//     }
//
//     // // ----------------- parte de pressão (seguindo exatamente o ek,ef) -----------------
//     // TPZFMatrix<STATE> Bp, Bpt;
//     // BuildBp(dphiP, Bp);
//     // Bp.Transpose(&Bpt);
//     //
//     // // (trial) vindos do material data
//     // TPZFMatrix<STATE> gradp_trial = dataP.dsol[0];                 // ∇p^{trial}
//     // const STATE       p_trial     = dataP.sol[0][0];               // p^{trial}
//     // const STATE       divu_trial  = dataU.dsol[0](0,0) + dataU.dsol[0](1,1); // tr(grad u^{trial})
//     //
//     // // histórico (n) vindos da memória
//     // auto m  = this->WithMem();
//     // const int gp = dataU.intGlobPtIndex;
//     // TMEM &mem = m->MemItem(gp);
//     //
//     // TPZFMatrix<STATE> gradp_n(dim,1,0.);
//     // for (int k=0; k<dim; ++k) gradp_n(k,0) = mem.m_elastoplastic_state.fdPorePressure[k];
//     // const STATE p_n    = mem.m_elastoplastic_state.fpressure;
//     // const STATE divu_n = mem.m_elastoplastic_state.fGradSolU.I1();
//     //
//     // // monta vetores nP×1 (sem gravidade, para bater com seu ek,ef)
//     // TPZFMatrix<STATE> Sp_trial(nP,1,0.), Hgradp_trial(nP,1,0.), Qt_u_trial(nP,1,0.);
//     // TPZFMatrix<STATE> Sp_n(nP,1,0.),    Hgradp_n(nP,1,0.),       Qt_u_n(nP,1,0.), tmp;
//     //
//     // // S p^{trial}  (S = fSe * φ φ^T)  -> como vetor: fSe * p_trial * φ
//     // for (int i=0; i<nP; ++i) Sp_trial(i,0) = fSe * p_trial * phiP(i,0);
//     //
//     // // dt H ∇p^{trial}  (H = (k/μ) Bp^T Bp)
//     // Bpt.Multiply(gradp_trial, Hgradp_trial);
//     // Hgradp_trial *= fTimeStep * (fk/fmu);
//     //
//     // // Q^T u^{trial} = α * div(u^{trial}) * φ
//     // for (int i=0; i<nP; ++i) Qt_u_trial(i,0) = falpha * divu_trial * phiP(i,0);
//     //
//     // // históricos
//     // for (int i=0; i<nP; ++i) Sp_n(i,0) = fSe * p_n * phiP(i,0);
//     //
//     // Bpt.Multiply(gradp_n, Hgradp_n);
//     // Hgradp_n *= fTimeStep * (fk/fmu);
//     //
//     // for (int i=0; i<nP; ++i) Qt_u_n(i,0) = falpha * divu_n * phiP(i,0);
//     //
//     // // F_p = (S+dtH) p^{trial} + Q^T u^{trial}  − [ S p^n + dtH ∇p^n + Q^T u^n ]
//     // TPZFMatrix<STATE> Fp(nP,1,0.);
//     // Fp  = Sp_trial;
//     // Fp += Hgradp_trial;
//     // Fp += Qt_u_trial;
//     // Fp -= Sp_n;
//     // Fp -= Hgradp_n;
//     // Fp -= Qt_u_n;
//     //
//     // // acumula no residual
//     // for (int i=0; i<nP; ++i)
//     // {
//     //     ef(off_p + i, 0) += weight * Fp(i,0);
//     // }
//     //
//     // #ifdef PZ_LOG
//     // // impressão leve de auditoria (um PG basta)
//     // //if (dataU.intGlobPtIndex == 0) {
//     //     std::ostringstream oss;
//     //     oss.setf(std::ios::scientific);
//     //     oss.precision(6);
//     //     const int npeqs = nP;
//     //     TPZFMatrix<STATE> fu1(nueqs,1,0.), fp1(npeqs,1,0.);
//     //     for (int i=0;i<nueqs;++i) fu1(i,0) = ef(off_u+i,0);
//     //     for (int i=0;i<nP;++i)    fp1(i,0) = ef(off_p+i,0);
//     //
//     //     oss << "Conribute(ef) "
//     //     << "|fu|_1=" << Norm(fu1) << "  |fp|_1=" << Norm(fp1)
//     //     << "|Bt Sigma|=" << Norm(fint1) << "  |Nt b|=" << Norm(fint2)
//     //     << "  dt=" << fTimeStep
//     //     << "  alpha=" << falpha
//     //     << "  k/mu=" << (fk/fmu)
//     //     << "  Se=" << fSe
//     //     << "  p_tr=" << p_trial << "  p_n=" << p_n
//     //     << "  divu_tr=" << divu_trial << "  divu_n=" << divu_n
//     //     << "\n";
//     //     LOGPZ_DEBUG(logger_poro, oss.str());
//     // //}
//     // #endif
// }


template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    const auto& phiU = datavec[0].phi;
    const auto& phiP = datavec[1].phi;
    const int nU = phiU.Rows();
    const int nP = phiP.Rows();

    const int off_u = 0;
    const int off_p = 2*nU;

    const STATE big = 1.e12;

    const TPZVec<STATE>& V2 = bc.Val2();

    switch (bc.Type())
    {
        // ---------------------------------------------------------
        // 0 : Dirichlet em u
        case 0:
        {

            for (int i=0;i<nU;++i){
                ef(off_u+2*i+0,0) += big * (V2[0]) * phiU(i,0) * weight;
                ef(off_u+2*i+1,0) += big * (V2[1]) * phiU(i,0) * weight;
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
                ef(in+off_p,0)	+= (V2[2]) *big*phiP(in,0)*weight;	// P Pressure Value
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
               // ef(2*in+0,0) += big * (0. - datavec[0].sol[0][0]) * V2[0] * phiU(in,0) * weight;
               // ef(2*in+1,0) += big * (0. - datavec[0].sol[0][1]) * V2[1] * phiU(in,0) * weight;
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

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeBC(datavec, weight, fakeek, ef, bc);
}
template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Errors(const TPZMaterialDataT<STATE>&data,
                                         TPZVec<REAL> &values)
{
    //TBase::Errors(data,values);
}
template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain)
{
    TPZFNMatrix<9> DSolXYZ(3, 3, 0.);
    data.axes.Multiply(data.dsol[0], DSolXYZ, 1/*transpose*/);
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }

    DeltaStrain(_XX_, 0) = DSolXYZ(0, 0);
    DeltaStrain(_YY_, 0) = DSolXYZ(1, 1);
    DeltaStrain(_XY_, 0) =  1/2.* (DSolXYZ(1, 0) + DSolXYZ(0, 1));
    DeltaStrain(_XZ_, 0) = 0.;
    DeltaStrain(_YZ_, 0) = 0.;
    DeltaStrain(_ZZ_, 0) = 0.;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
                                        TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticityModel);

    // #ifdef PZ_LOG
    // if (logger_poro.isDebugEnabled()) {
    //     std::ostringstream oss;
    //     oss << "intPt=" <<intPt  << "\n";
    //     LOGPZ_INFO(logger_poro, oss.str());
    // }
    // #endif
    /// Access to memory data
    plasticloc.SetState(this->MemItem(intPt).m_elastoplastic_state);
    TPZTensor<REAL> eps_t, sigma(this->MemItem(intPt).m_sigma);
    eps_t.CopyFrom(DeltaStrain);
    eps_t.Add(plasticloc.GetState().m_eps_t, 1.);

    plasticloc.SetElasticResponse(this->MemItem(intPt).m_ER);

    plasticloc.ApplyStrainComputeSigma(eps_t, sigma, &Dep);

    sigma.CopyTo(Stress);
    TPZVec<STATE> phi;
    plasticloc.Phi(eps_t,phi);


    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        TMEM &mem=this->MemItem(intPt);
        mem.m_sigma        = sigma;
        mem.m_elastoplastic_state = plasticloc.GetState(); // unidades consistentes
        mem.m_plastic_steps= plasticloc.IntegrationSteps();      // rad    // rad
        mem.m_elastoplastic_state.fmatprop =plasticloc.GetState().fmatprop;
        mem.m_phi =phi[0];
        mem.m_ER = plasticloc.GetElasticResponse();

    }
}
template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
                              TPZFMatrix<REAL> & Stress)
{
    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticityModel);
    // #ifdef PZ_LOG
    // if (logger_poro.isDebugEnabled()) {
    //     std::ostringstream oss;
    //     oss << "intPt=" <<intPt  << "\n";
    //     LOGPZ_INFO(logger_poro, oss.str());
    // }
    // #endif
    /// Access to memory data
    plasticloc.SetState(this->MemItem(intPt).m_elastoplastic_state);
    TPZTensor<REAL> eps_t, sigma(this->MemItem(intPt).m_sigma);
    eps_t.CopyFrom(DeltaStrain);
    eps_t.Add(plasticloc.GetState().m_eps_t, 1.);


    plasticloc.SetElasticResponse(this->MemItem(intPt).m_ER);

    plasticloc.ApplyStrainComputeSigma(eps_t, sigma);


    TPZVec<STATE> phi;
    plasticloc.Phi(eps_t,phi);
    sigma.CopyTo(Stress);

    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        TMEM &mem=this->MemItem(intPt);
        mem.m_sigma        = sigma;
        mem.m_elastoplastic_state = plasticloc.GetState(); // unidades consistentes
        mem.m_plastic_steps= plasticloc.IntegrationSteps();      // rad    // rad
        mem.m_elastoplastic_state.fmatprop =plasticloc.GetState().fmatprop;
        mem.m_phi =phi[0];
        mem.m_ER = plasticloc.GetElasticResponse();

    }
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::SetElasticResponse(const TPZElasticResponse &ER)
{
    TMEM m;
    m.m_ER.SetEngineeringData(ER.E(), ER.Poisson());
    this->SetDefaultMem(m);
}
template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::BuildBu(const TPZFMatrix<STATE>& dphiU, TPZFMatrix<STATE>& Bu)
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

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::BuildNu(const TPZFMatrix<STATE>& phiU, TPZFMatrix<STATE>& Nu)
{
    const int nU = phiU.Rows();
    Nu.Redim(2, 2*nU);
    Nu.Zero();
    for(int a=0;a<nU;a++)
    {
        const REAL N = phiU(a,0);
        Nu(0, 2*a ) = N; // ux shape
        Nu(1, 2*a + 1) = N; // uy shape
    }
}


template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::BuildBp(const TPZFMatrix<STATE>& dphiP, TPZFMatrix<STATE>& Bp)
{
    const int nP = dphiP.Cols();
    Bp.Redim(2, nP);
    for (int j=0; j<nP; ++j) {
        Bp(0,j) = dphiP(0,j);     // dNj/dx
        Bp(1,j) = dphiP(1,j);     // dNj/dy
    }
}


template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::SetPlasticModel(T & plasticmodel)
{
    TMEM memory;
    fPlasticityModel = plasticmodel;
    T plastloc(fPlasticityModel);

    memory.m_elastoplastic_state = plastloc.GetState();
    //memory.m_ER=plastloc.GetElasticResponse();
    plastloc.ApplyStrainComputeSigma(memory.m_elastoplastic_state.m_eps_t, memory.m_sigma);

    this->SetDefaultMem(memory);
    //memory.Print(std::cout);

}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::SetBulkDensity(REAL & RhoB)
{
    //fRhoBulkDensity = RhoB;
}

template <class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::SetPorousElasticity(TPZElasticResponse & PER){
    m_PER = PER;

}

template <class T, class TMEM>
TPZElasticResponse TPZMatPoroElastoPlastic2DMem<T,TMEM>::GetPorousElasticity(TPZElasticResponse & PER){
    return PER;
}

template <class T, class TMEM>
T & TPZMatPoroElastoPlastic2DMem<T,TMEM>::GetPlasticModel(){
    return fPlasticityModel;
}


// --- Name() ---
template<class T, class TMEM>
std::string TPZMatPoroElastoPlastic2DMem<T,TMEM>::Name() const {
    return "TPZMatPoroElastoPlastic2DMem";
}

// --- NewMaterial() ---
template<class T, class TMEM>
TPZMaterial * TPZMatPoroElastoPlastic2DMem<T,TMEM>::NewMaterial() const {
    return new TPZMatPoroElastoPlastic2DMem<T,TMEM>(*this);
}

// --- ClassId() --- (qualquer inteiro determinístico está ok)
template<class T, class TMEM>
int TPZMatPoroElastoPlastic2DMem<T,TMEM>::ClassId() const {
    return Hash("TPZMatPoroElastoPlastic2DMem");
}

// --- Write/Read mínimos (se quiser evitar DebugStop) ---
template<class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Write(TPZStream &buf, int withclassid) const {
    TBase::Write(buf, withclassid);
    buf.Write(&fPlaneStrain,1);
    buf.Write(&fFactor,1);
    buf.Write(fBody);
    // se necessário, serialize fPlasticityModel e m_PER
}

template<class T, class TMEM>
void TPZMatPoroElastoPlastic2DMem<T,TMEM>::Read(TPZStream &buf, void *context) {
    TBase::Read(buf, context);
    buf.Read(&fPlaneStrain,1);
    buf.Read(&fFactor,1);
    buf.Read(fBody);

    // se necessário, deserializar fPlasticityModel e m_PER
}


template class TPZMatPoroElastoPlastic2DMem<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>,TPZPorousElastoPlasticMem>;
template class TPZMatPoroElastoPlastic2DMem<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>,TPZElastoPlasticMem>;
