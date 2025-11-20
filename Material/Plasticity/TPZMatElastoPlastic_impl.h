#include "TPZMatElastoPlastic.h"
#include "TPZBndCondT.h"
#include "TPZMaterialDataT.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger elastoplasticLogger("pz.material.pzElastoPlastic");
static TPZLogger updatelogger("pz.material.pzElastoPlastic.update");
static TPZLogger ceckconvlogger("checkconvmaterial");
#endif


template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic() : TBase(), m_force(), m_rho_bulk(0), m_PostProcessDirection(), m_tol(1.e-6)
{
    m_force.Resize(3,0);
    m_force[1] = -9.8; // proper gravity acceleration in m/s^2
    m_PostProcessDirection.Resize(3,0);
    m_PostProcessDirection[0] = 1.;
    m_use_non_linear_elasticity_Q = false;

    #ifdef PZ_LOG
    if(elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif

}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic(int id) : TBase(id), m_force(), m_rho_bulk(0), m_PostProcessDirection(), m_tol(1.e-6)
{
    m_force.Resize(3,0);
    //m_force[1] = -9.8; // proper gravity acceleration in m/s^2 -> 1=y 0=x 2=z
    m_PostProcessDirection.Resize(3,0);
    m_PostProcessDirection[0] = 1.;
    m_use_non_linear_elasticity_Q = false;
    TPZPlasticState<STATE> def;


    #ifdef PZ_LOG
    if (elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>(int id) constructor called with id = " << id << " ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif

}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic(const TPZMatElastoPlastic &other) : TBase(other),
m_force(other.m_force), m_rho_bulk(other.m_rho_bulk), m_PostProcessDirection(other.m_PostProcessDirection),
m_plasticity_model(other.m_plasticity_model), m_tol(other.m_tol), m_PER(other.m_PER)
{
    #ifdef PZ_LOG
    if(elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() copy constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif
    m_use_non_linear_elasticity_Q = other.m_use_non_linear_elasticity_Q;
}
template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetPlasticityModel(T & plasticity)
{
    TMEM memory;
    m_plasticity_model = plasticity;
    T plastloc(m_plasticity_model);
    memory.m_elastoplastic_state = plastloc.GetState();
    memory.m_ER=plastloc.GetElasticResponse();
    plastloc.ApplyStrainComputeSigma(memory.m_elastoplastic_state.m_eps_t, memory.m_sigma);
    this->SetDefaultMem(memory);

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetBulkDensity(REAL & RhoB)
{
    m_rho_bulk = RhoB;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetPorousElasticity(TPZPorousElasticResponse & PER){
    m_PER = PER;
    m_use_non_linear_elasticity_Q = true;
}

template <class T, class TMEM>
TPZPorousElasticResponse & TPZMatElastoPlastic<T,TMEM>::GetPorousElasticity(TPZPorousElasticResponse & PER){
    return PER;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetPlasticModel(T & plasticity_model){
    m_plasticity_model = plasticity_model;
}

template <class T, class TMEM>
T & TPZMatElastoPlastic<T,TMEM>::GetPlasticModel(){
    return m_plasticity_model;
}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::~TPZMatElastoPlastic()
{

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Print(std::ostream &out, const int memory) const
{
    out << this->Name();
    out << "\n with template argurment T = " << m_plasticity_model.Name();
    out << "\n Base material Data:\n";
    TPZMatWithMem<TMEM>::PrintMem(out, memory);
    out << "\n Localy defined members:";
    out << "\n Body Forces: " << m_force;
    out << "\n Bulk density = " << m_rho_bulk;
    out << "\n Post process direction: " << m_PostProcessDirection;
    out << "\n Tolerance for internal post processing iterations: " << m_tol;
    out << "\n Directive for the use of nonlinear elasticity: " << m_use_non_linear_elasticity_Q;
    out << "\n Internal plasticity <T> member:\n";
    m_plasticity_model.Print(out);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << this->Name();
    TPZMatWithMem<TMEM>::Print(out);
    out << "\n Body Forces: " << m_force;
    out << "\n Bulk density = " << m_rho_bulk;
    out << "\n Post process direction: " << m_PostProcessDirection;
    out << "\n Tolerance for internal post processing iterations: " << m_tol;
    out << "\n Directive for the use of nonlinear elasticity: " << m_use_non_linear_elasticity_Q;
    out << "\n Internal plasticity <T> member:\n";
    m_plasticity_model.Print(out);
}

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::VariableIndex(const std::string &name) const
{
    if(!strcmp("DisplacementDoF",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EDisplacementDoFx;
    if(!strcmp("Displacement",name.c_str()))            return TPZMatElastoPlastic<T,TMEM>::EDisplacement;
    if(!strcmp("Strain",name.c_str()))                  return TPZMatElastoPlastic<T,TMEM>::EStrain;
    if(!strcmp("StressXX",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::ESX;
    if(!strcmp("StressYY",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::ESY;
    if(!strcmp("StressZZ",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::ESZ;
    if(!strcmp("StrainPlasticZZ",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EEPZ;
    if(!strcmp("StrainTotalZZ",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EEPZT;
    if(!strcmp("StrainElasticZZ",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EEEZ;
    if(!strcmp("StrainPlasticXX",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EEPX;
    if(!strcmp("StrainPlasticYY",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EEPY;
    if(!strcmp("StrainPlastic",name.c_str()))           return TPZMatElastoPlastic<T,TMEM>::EStrainPlastic;
    if(!strcmp("Yield",name.c_str()))                   return TPZMatElastoPlastic<T,TMEM>::EYield;
    if(!strcmp("VolHardening",name.c_str()))            return TPZMatElastoPlastic<T,TMEM>::EVolHardening;
    if(!strcmp("StrainPValues",name.c_str()))           return TPZMatElastoPlastic<T,TMEM>::EStrainPValues;
    if(!strcmp("StressPValues",name.c_str()))           return TPZMatElastoPlastic<T,TMEM>::EStressPValues;
    if(!strcmp("StrainElasticPValues",name.c_str()))    return TPZMatElastoPlastic<T,TMEM>::EStrainElasticPValues;
    if(!strcmp("StrainPlasticPValues",name.c_str()))    return TPZMatElastoPlastic<T,TMEM>::EStrainPlasticPValues;
    if(!strcmp("StrainI1",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EStrainI1;
    if(!strcmp("StressI1",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EStressI1;
    if(!strcmp("StrainElasticI1",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EStrainElasticI1;
    if(!strcmp("StrainPlasticI1",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EStrainPlasticI1;
    if(!strcmp("StrainJ2",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EStrainJ2;
    if(!strcmp("StressJ2",name.c_str()))                return TPZMatElastoPlastic<T,TMEM>::EStressJ2;
    if(!strcmp("StrainElasticJ2",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EStrainElasticJ2;
    if(!strcmp("StrainPlasticJ2",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EStrainPlasticJ2;
    if(!strcmp("FailureType",name.c_str()))             return TPZMatElastoPlastic<T,TMEM>::EFailureType;
    if(!strcmp("DamageVariable",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>:: EDamageVar;
    if(!strcmp("Exact",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EEXACT;
    if(!strcmp("EBodyForce",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EBodyForce;
    if(!strcmp("EOrder",name.c_str()))         return TPZMatElastoPlastic<T,TMEM>::EOrder;
    PZError << "TPZMatElastoPlastic<T,TMEM>:: VariableIndex Error\n";
    return TPZMatElastoPlastic<T,TMEM>::ENone;
}

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::NSolutionVariables(int var) const
{
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementDoFx) return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacement) return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrain) return 6;
    if(var == TPZMatElastoPlastic<T,TMEM>::ESX) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::ESY) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::ESZ) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEPZ) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEPZT) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEPX) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEPY) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEEZ) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainPlastic) return 6;
    if(var == TPZMatElastoPlastic<T,TMEM>::EYield) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EVolHardening) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainPValues) return 6;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStressPValues) return 6;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainElasticPValues) return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainPlasticPValues) return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainI1) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainElasticI1) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainPlasticI1) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainJ2) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStressJ2) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainElasticJ2) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EStrainPlasticJ2) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EFailureType) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EEXACT) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDamageVar) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EBodyForce) return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EOrder) return 1;
    if(var == 100) return 1;
    return TBase::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDirection(TPZFMatrix<REAL> &vectorTensor, TPZVec<REAL> &Out)
{
    Out.Resize(3);
    TPZVec<REAL> &Dir = this->m_PostProcessDirection;
    Out[0] = Dir[0] * vectorTensor(_XX_,0) + Dir[1] * vectorTensor(_XY_,0) + Dir[2] * vectorTensor(_XZ_,0);
    Out[1] = Dir[0] * vectorTensor(_XY_,0) + Dir[1] * vectorTensor(_YY_,0) + Dir[2] * vectorTensor(_YZ_,0);
    Out[2] = Dir[0] * vectorTensor(_XZ_,0) + Dir[1] * vectorTensor(_YZ_,0) + Dir[2] * vectorTensor(_ZZ_,0);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) {

        int intPt = data.intGlobPtIndex;
        if (intPt == -1 || TPZMatElastoPlastic<T, TMEM>::GetMemory()->NElements() == 0) {
                return;
        }
        TMEM &Memory = this->MemItem(intPt);
    /// Displacements from Degree of Freedom



    T plasticloc(m_plasticity_model);
    plasticloc.SetState(Memory.m_elastoplastic_state);
    switch (var) {
        case EDisplacement:
        {
            Solout.Resize(3);
            //for (int i = 0; i < 3; i++) {
                //Solout[i] = Memory.m_u[i];
                const REAL ux = data.sol[0][0];
                const REAL uy = data.sol[0][1];
                const REAL uz = data.sol[0][2];
                Solout[0]=ux; Solout[1]=uy; Solout[2]=uz;
            //}
        }
        break;
        case EDisplacementDoFx:
        {
                Solout.Resize(3);
                Solout[0] = Memory.m_u[0];
                Solout[1] = Memory.m_u[1];
                Solout[2] = Memory.m_u[2];

        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::ESX:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & sigma = Memory.m_sigma;
            Solout[0] = sigma.XX();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::ESY:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & sigma = Memory.m_sigma;
            Solout[0] = sigma.YY();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::ESZ:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & sigma = Memory.m_sigma;
            Solout[0] = sigma.ZZ();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEPZ:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.ZZ();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEPZT:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_t;
            Solout[0] = eps_p.ZZ();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEPX:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.XX();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEPY:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.YY();


        }
        break;
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEEZ:
        {
            Solout.Resize(1);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            TPZTensor<REAL> & eps_t = Memory.m_elastoplastic_state.m_eps_t;
            TPZTensor<REAL> eps_e=eps_t-eps_p;
            Solout[0] = eps_e.ZZ();


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EDamageVar:
        {
            Solout.Resize(1);
            STATE hardening= Memory.m_elastoplastic_state.m_hardening;
            Solout[0] =hardening;


        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainPlastic:
        {
            Solout.Resize(6);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.XX();
            Solout[1] = eps_p.YY();
            Solout[2] = eps_p.ZZ();
            Solout[3] = eps_p.XY();
            Solout[4] = eps_p.XZ();
            Solout[5] = eps_p.YZ();

        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainPValues:
        {
            Solout.Resize(3);
            TPZTensor<REAL> & eps = Memory.m_elastoplastic_state.m_eps_t;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            eps.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++)Solout[i] = eigensystem.fEigenvalues[i];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStressPValues:
        {
            Solout.Resize(3);
            TPZTensor<REAL> & Sigma = Memory.m_sigma;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            Sigma.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++)Solout[i] = eigensystem.fEigenvalues[i];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainElasticPValues:
        {
            Solout.Resize(3);
            TPZTensor<REAL> eps_e(Memory.m_elastoplastic_state.m_eps_t);
            eps_e -= Memory.m_elastoplastic_state.m_eps_p;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            eps_e.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++) Solout[i] = eigensystem.fEigenvalues[i];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainPlasticPValues:
        {
            Solout.Resize(3);
            TPZTensor<REAL> & eps_p = Memory.m_elastoplastic_state.m_eps_p;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            eps_p.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++) Solout[i] = eigensystem.fEigenvalues[i];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainI1:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_t = Memory.m_elastoplastic_state.m_eps_t;
            Solout[0] = eps_t.I1();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStressI1:
        {
            Solout.Resize(1);
            TPZTensor<REAL> sigma = Memory.m_sigma;
            Solout[0] = sigma.I1();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainElasticI1:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_e(Memory.m_elastoplastic_state.m_eps_t);
            eps_e -= Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_e.I1();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainPlasticI1:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.I1();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainJ2:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_t = Memory.m_elastoplastic_state.m_eps_t;
            Solout[0] = eps_t.J2();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStressJ2:
        {
            Solout.Resize(1);
            TPZTensor<REAL> sigma = Memory.m_sigma;
            Solout[0] = sigma.J2();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainElasticJ2:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_e(Memory.m_elastoplastic_state.m_eps_t);
            eps_e -= Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_e.J2();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EStrainPlasticJ2:
        {
            Solout.Resize(1);
            TPZTensor<REAL> eps_p = Memory.m_elastoplastic_state.m_eps_p;
            Solout[0] = eps_p.J2();
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EYield:
        {

            TPZTensor<REAL> & EpsT = Memory.m_elastoplastic_state.m_eps_t;
            TPZTensor<STATE> epsElastic(EpsT);
            epsElastic -= Memory.m_elastoplastic_state.m_eps_p;
            plasticloc.Phi(epsElastic, Solout);
        }//EVolPlasticSteps - makes sense only if the evaluated point refers to an identified integration point
        break;
        case TPZMatElastoPlastic<T, TMEM>::EFailureType:
        {
            Solout.Resize(1);
            int m_type = Memory.m_elastoplastic_state.m_m_type;
            Solout[0] = m_type;
        }
        break;
        break;
        case TPZMatElastoPlastic<T, TMEM>::EEXACT:
        {
            TPZVec<STATE> uex(2, 0.0);
            TPZFMatrix<STATE> duex(2, 2, 0.0);
            fExactSolution(data.x, uex, duex);
            Solout.Resize(2);
            // Comparar deslocamentos FEM e exatos
            TPZVec<STATE> uh(2);
            Solout[0] = data.sol[0][0];
            Solout[1] = data.sol[1][0];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EBodyForce:
        {
            Solout.Resize(1);
            Solout[0] = m_force[1];
        }
        break;
        case TPZMatElastoPlastic<T, TMEM>::EOrder:
        {
            Solout.Resize(1);
            Solout[0] = data.p;
        }
        break;
        default:
        {
            TBase::Solution(data, var, Solout);
        }
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Contribute(const TPZMaterialDataT<STATE> &data,
                                             REAL weight,
                                             TPZFMatrix<REAL> &ek,
                                             TPZFMatrix<REAL> &ef)
{

    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &phi  = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    const TPZManVector<REAL,3> &x = data.x;
    TPZFMatrix<REAL> axesT, dphiXYZ;
    // rotating the shape functions to the XYZ coordinates
    axes.Transpose(&axesT);
    axesT.Multiply(dphi,dphiXYZ);

    TPZFNMatrix<9>  Deriv(3,3);
    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);

    this->ComputeDeltaStrainVector(data, DeltaStrain);
    this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);

    TPZManVector<STATE, 3> ForceLoc(3,0.0);
    if(this->fForcingFunction)
    {
        this->fForcingFunction(data.x,ForceLoc);
    }
    // força de corpo (gravidade), já em unidades de volume
    TPZFMatrix<STATE> flocal(3,1,0.);
    flocal(0,0) = ForceLoc[0];
    flocal(1,0) = ForceLoc[1];
    flocal(2,0) = ForceLoc[2];
    flocal *= m_rho_bulk; // se o ForceLoc é aceleração (g), multiplica pela densidade

    // constrói B e N
    TPZFMatrix<STATE> B, N, Bt, Nt, temp, ektemp, fint1, fint2;
    BuildBN(dphiXYZ, phi, B, N);

    // Bt * Dep * B
    B.Transpose(&Bt);                 // (3*phr x 6)
    Bt.Multiply(Dep, temp);           // (3*phr x 6) * (6 x 6) = (3*phr x 6)
    temp.Multiply(B, ektemp);         // (3*phr x 6) * (6 x 3*phr) = (3*phr x 3*phr)
    ektemp *= weight;
    ek += ektemp;

    // forças: ef += (N * fbody  -  Bt * sigma) * weight
    N.Multiply(flocal, fint2);        // (3*phr x 3) * (3 x 1) = (3*phr x 1)
    Bt.Multiply(Stress, fint1);       // (3*phr x 6) * (6 x 1) = (3*phr x 1)
    fint2 -= fint1;
    ef += fint2 * weight;

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                               REAL weight,
                                               TPZFMatrix<REAL> &ek,
                                               TPZFMatrix<REAL> &ef,
                                               TPZBndCondT<STATE> &bc)
{
    #ifdef PZ_LOG
    if(elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::ContributeBC *** with bc.Type()=" << bc.Type();
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif
    const TPZFMatrix<REAL> &phi = data.phi;

    const REAL BIGNUMBER  = 1.e16;

    int dim = Dimension();
    int nstate = NStateVariables();

    const int phr = phi.Rows();
    int in,jn,idf,jdf;
    REAL v2[3];
    v2[0] = bc.Val2()[0];
    v2[1] = bc.Val2()[1];
    v2[2] = bc.Val2()[2];


    const TPZFMatrix<REAL> &v1 = bc.Val1();
    //bc.Print(cout);
    //std::cout << "val2:  " << v2[2] << std::endl;
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (v2[2] - data.sol[0][2]) * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                }//jn
            }//in
            break;

        case 1: // Neumann condition
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
            }
            break;

        case 2: // Mixed condition
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
                for(jn=0; jn<phi.Rows(); jn++)
                {
                    for(idf=0; idf<3; idf++) for(jdf=0; jdf<3; jdf++)
                    {
                        ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf);
                    }
                }
            }//in
            break;

        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (0. - data.sol[0][2]) * v2[2] * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[2];
                }//jn
            }//in
            break;

        case 4: // stressField Neumann condition
            for(in = 0; in < dim; in ++)
                v2[in] = - ( v1(in,0) * data.normal[0] +
                v1(in,1) * data.normal[1] +
                v1(in,2) * data.normal[2] );
            // The normal vector points towards the neighbour. The negative sign is there to
            // reflect the outward normal vector.
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
                //    cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                //    cout << "val2:  " << v2[0]  << endl;
            }
            break;

        case 5://PRESSAO
            for(in = 0 ; in < phi.Rows(); in++)
            {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight * (data.normal[0]);
                ef(nstate*in+1,0) += v2[0] * phi(in,0) * weight * (data.normal[1]);
                ef(nstate*in+2,0) += v2[0] * phi(in,0) * weight * (data.normal[2]);
            }
            break;

        case 6: // Directional Dirichlet - displacement is set to u_D in the non-null vector component direction

            REAL v_null[3];
            v_null[0] = bc.Val1()(0,0);
            v_null[1] = bc.Val1()(1,1);
            v_null[2] = bc.Val1()(2,2);

            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * v_null[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * v_null[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (v2[2] - data.sol[0][2]) * v_null[2] * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[0];
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[1];
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[2];
                }//jn
            }//in
            break;

        default:
            #ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
        #endif
        PZError << "TPZMatElastoPlastic::ContributeBC error - Wrong boundary condition type" << std::endl;
    }//switch

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Contribute(const TPZMaterialDataT<STATE> &data,
                                             REAL weight, TPZFMatrix<REAL> &ef)
{
    #ifdef PZ_LOG
    if(elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Point index = " << data.intGlobPtIndex;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif

    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZFMatrix<REAL> &phi  = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    const TPZManVector<REAL,3> &x = data.x;
    TPZFMatrix<REAL> axesT, dphiXYZ;
    // rotating the shape functions to the XYZ coordinates
    axes.Transpose(&axesT);
    axesT.Multiply(dphi,dphiXYZ);
    const int phr = phi.Rows();

    TPZFNMatrix<9>  Deriv(3,3);
    //    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);

    this->ComputeDeltaStrainVector(data, DeltaStrain);
    this->ApplyDeltaStrain(data,DeltaStrain,Stress);

    int nstate = NStateVariables();
    REAL val;

    TPZManVector<STATE, 3> ForceLoc(nstate,0.0);
    if(this->fForcingFunction)
    {
        this->fForcingFunction(data.x,ForceLoc);
    }

    int in;
    for(in = 0; in < phr; in++) { //in: test function index

        // m_force represents the gravity acceleration
        //First equation: fb and fk
        val  = m_rho_bulk * ForceLoc[0] * phi(in,0); // fb
        val -= Stress(_XX_,0) * dphiXYZ(0,in); // |
        val -= Stress(_XY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_XZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+0,0) += weight * val;

        //Second equation: fb and fk
        val  = m_rho_bulk * ForceLoc[1] * phi(in,0); // fb
        val -= Stress(_XY_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_YZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+1,0) += weight * val;

        //third equation: fb and fk
        val  = m_rho_bulk * ForceLoc[2] * phi(in,0); // fb
        val -= Stress(_XZ_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YZ_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_ZZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+2,0) += weight * val;

    }//in

    #ifdef PZ_LOG
    if(elastoplasticLogger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        //sout << " Resultant rhs vector:\n" << ef;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    #endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                               REAL weight,
                                               TPZFMatrix<REAL> &ef,
                                               TPZBndCondT<STATE> &bc)
{
    TPZFMatrix<REAL> fakeek(ef.Rows(),ef.Rows(),0.);
    ContributeBC(data,weight,fakeek,ef,bc);
    //TBase::ContributeBC(data, weight, ef, bc);//not efficient but here to remember reimplementing it when ContributeBC becomes robust
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Errors(const TPZMaterialDataT<STATE>&data,
                                         TPZVec<REAL> &values)
{
    const auto &x = data.x;
    const auto &u = data.sol[0];
    const auto &dudx = data.dsol[0];

    #ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
    #endif
    TPZManVector<STATE,3> u_exact(3);
    TPZFNMatrix<9,STATE> du_exact(3,3,0.);
    this->fExactSol(x,u_exact,du_exact);
    int i, j;

    /** L2 norm */
    REAL L2 = 0.;
    for(i = 0; i < 3; i++) L2 += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);

    /** H1 semi-norm */
    REAL SemiH1 = 0.;
    for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) SemiH1 += (dudx(i,j) - du_exact(i,j)) * (dudx(i,j) - du_exact(i,j));

    /** H1 norm */
    REAL H1 = L2 + SemiH1;

    //values[1] : eror em norma L2
    values[1]  = L2;

    //values[2] : erro em semi norma H1
    values[2] = SemiH1;

    //values[0] : erro em norma H1 <=> norma Energia
    values[0]  = H1;

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &Strain)
{
    ComputeDeltaStrainVector(data, Strain);

    TPZTensor<REAL> & EpsT = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_t;

    int i;
    for( i = 0; i < 6; i++ )Strain(i,0) = Strain(i,0) + EpsT.fData[i];
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain)
{
    TPZFNMatrix<9> DSolXYZ(3,3,0.);
    data.axes.Multiply(data.dsol[0],DSolXYZ,1/*transpose*/);
    //std::cout << "\n data dsol \n";
    //data.dsol[0].Print("data.sol");
    DeltaStrain.Redim(6,1);
    DeltaStrain(_XX_,0) = DSolXYZ(0,0);
    DeltaStrain(_YY_,0) = DSolXYZ(1,1);
    DeltaStrain(_ZZ_,0) = DSolXYZ(2,2);
    DeltaStrain(_XY_,0) =( DSolXYZ(1,0) + DSolXYZ(0,1) );
    DeltaStrain(_XZ_,0) = ( DSolXYZ(2,0) + DSolXYZ(0,2) );
    DeltaStrain(_YZ_,0) = ( DSolXYZ(2,1) + DSolXYZ(1,2) );
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeStressVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &Stress)
{

    TPZFNMatrix<6> DeltaStrain;
    ComputeDeltaStrainVector(data, DeltaStrain);
    Stress.Redim(6,1);
    ApplyDeltaStrain(data, DeltaStrain, Stress);

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::CheckConvergence(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain)
{
    // f[x_] := ComputeN[x]; df[x_] := GradNdSig[x];
    // x0 = RandomReal[{0, 6}, 6]; dx = RandomReal[{0, 0.01}, 6];
    // err[\[Alpha]_] :=
    // Abs[f[x0 + \[Alpha] dx] - (f[x0] + \[Alpha] (df[x0] . dx))];
    // \[Alpha]1 = RandomReal[{0, 0.01}]; \[Alpha]2 = RandomReal[{0, 0.01}];
    // p = Log[err[\[Alpha]2]/err[\[Alpha]1]]/Log[\[Alpha]2/\[Alpha]1]; {x0,
    //     dx, err[\[Alpha]1], err[\[Alpha]2], p} // TableForm
    int intPt = data.intGlobPtIndex;//, plasticSteps;
    T plasticloc(m_plasticity_model);
    plasticloc.SetState(this->MemItem(intPt).m_elastoplastic_state);
    TPZTensor<REAL> deps;
    deps.CopyFrom(DeltaStrain);

    REAL alfa =1.e-6;
    REAL alfa2 = 2.e-6;
    TPZTensor<REAL> part1,part2,part3,temp;
    TPZTensor<REAL> Alfa1DeltaEps, Alfa2DeltaEps,Eps(plasticloc.GetState().m_eps_t);
    Alfa1DeltaEps.CopyFrom(DeltaStrain);
    Alfa2DeltaEps.CopyFrom(DeltaStrain);
    TPZFNMatrix<36,REAL> DEP(6,6);

    Alfa1DeltaEps*=alfa;
    Alfa2DeltaEps*=alfa2;
    temp=Eps;
    temp+=Alfa1DeltaEps;
    plasticloc.ApplyStrainComputeSigma(temp,part1);
    plasticloc.ApplyStrainComputeDep(Eps,part2,DEP);
    TPZFNMatrix<6,REAL> part3temp(6,1),tempAlfa1DeltaEps(6,1);
    for(int i=0;i<6;i++)
    {
        tempAlfa1DeltaEps(i,0)=Alfa1DeltaEps.fData[i];
    }
    DEP.Multiply(tempAlfa1DeltaEps, part3temp);
    part3.CopyFrom(part3temp);
    TPZTensor<REAL> e1(part1);
    e1-=part2;
    e1-=part3;

    part1*=0.;
    part3*=0.;
    part3temp*=0.;
    temp*=0.;

    temp=Eps;
    temp+=Alfa2DeltaEps;
    plasticloc.ApplyStrainComputeSigma(temp,part1);
    for(int i=0;i<6;i++)
    {
        tempAlfa1DeltaEps(i,0)=Alfa2DeltaEps.fData[i];
    }
    DEP.Multiply(tempAlfa1DeltaEps, part3temp);
    part3.CopyFrom(part3temp);
    TPZTensor<REAL> e2(part1);
    e2-=part2;
    e2-=part3;
    REAL n = (log10(Norm(e1))-log10(Norm(e2)))/(log10(alfa)-log10(alfa2));
    TPZManVector<REAL,3> phi(3,1);
    plasticloc.Phi(Eps,phi);
    std::cout << "DEP "<< DEP << std::endl;
    std::cout << "tempAlfa1DeltaEps "<< tempAlfa1DeltaEps << std::endl;
    std::cout << "Phi "<< phi << std::endl;
    std::cout << "Integration Point "<< intPt << std::endl;
    std::cout << "n = " << n << std::endl;
    #ifdef PZ_LOG
    if(ceckconvlogger.isDebugEnabled())
    {
        std::stringstream sout;
        TPZManVector<REAL,3> phi(3,1);
        plasticloc.Phi(Eps,phi);
        sout << "DEP "<< DEP << std::endl;
        sout << "tempAlfa1DeltaEps "<< tempAlfa1DeltaEps << std::endl;
        sout << "Phi "<< phi << std::endl;
        sout << "Integration Point "<< intPt << std::endl;
        sout << "n = " << n << std::endl;
        LOGPZ_DEBUG(ceckconvlogger, sout.str())
    }
    #endif



}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
                                                             TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{

    int intPt = data.intGlobPtIndex;
    T plasticloc(m_plasticity_model);

    /// Access to memory data
    plasticloc.SetState(this->MemItem(intPt).m_elastoplastic_state);
    TPZTensor<REAL> eps_t, sigma(this->MemItem(intPt).m_sigma);
    eps_t.CopyFrom(DeltaStrain);
    eps_t.Add(plasticloc.GetState().m_eps_t, 1.);

    if (m_use_non_linear_elasticity_Q) {
        //        TPZTensor<REAL> & last_eps_t = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_t;
        TPZTensor<REAL> & last_eps_p = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_p;
        TPZTensor<REAL> eps_e = eps_t - last_eps_p;
        this->MemItem(intPt).m_ER = m_PER.EvaluateElasticResponse(eps_e);
    }
    plasticloc.SetElasticResponse(this->MemItem(intPt).m_ER);

    UpdateMaterialCoeficients(data.x,plasticloc);
    //std::cout<< "eps_t = " << eps_t <<std::endl;
    plasticloc.ApplyStrainComputeSigma(eps_t, sigma, &Dep);

    sigma.CopyTo(Stress);

    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).m_sigma        = sigma;
        //this->MemItem(intPt).m_elastoplastic_state = plasticloc.GetState();
        this->MemItem(intPt).m_elastoplastic_state.m_eps_t = plasticloc.GetState().m_eps_t;
        this->MemItem(intPt).m_elastoplastic_state.m_eps_p = plasticloc.GetState().m_eps_p;
        this->MemItem(intPt).m_elastoplastic_state.m_hardening = plasticloc.GetState().m_hardening;
        this->MemItem(intPt).m_plastic_steps = plasticloc.IntegrationSteps();
        this->MemItem(intPt).m_ER = plasticloc.GetElasticResponse();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).m_u[i] = data.sol[0][i];
        }
    }

}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::UpdateMaterialCoeficients(const TPZVec<REAL> &x,T & plasticity)
{
    /// Nothing to do
    return;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & Strain,
                                                   TPZFMatrix<REAL> & Stress)
{

    int intPt = data.intGlobPtIndex;
    T plasticloc(m_plasticity_model);

    /// Access to memory data
    plasticloc.SetState(this->MemItem(intPt).m_elastoplastic_state);
    TPZTensor<REAL> eps_t, sigma(this->MemItem(intPt).m_sigma);
    eps_t.CopyFrom(Strain);
    eps_t.Add(plasticloc.GetState().m_eps_t, 1.);

    if (m_use_non_linear_elasticity_Q) {
        //        TPZTensor<REAL> & last_eps_t = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_t;
        TPZTensor<REAL> & last_eps_p = this->MemItem(data.intGlobPtIndex).m_elastoplastic_state.m_eps_p;
        TPZTensor<REAL> eps_e = eps_t - last_eps_p;
        this->MemItem(intPt).m_ER = m_PER.EvaluateElasticResponse(eps_e);
    }
    plasticloc.SetElasticResponse(this->MemItem(intPt).m_ER);

    UpdateMaterialCoeficients(data.x,plasticloc);
    plasticloc.ApplyStrainComputeSigma(eps_t, sigma);
    sigma.CopyTo(Stress);

    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).m_sigma        = sigma;
        //this->MemItem(intPt).m_elastoplastic_state = plasticloc.GetState();
        this->MemItem(intPt).m_elastoplastic_state.m_eps_t = plasticloc.GetState().m_eps_t;
        this->MemItem(intPt).m_elastoplastic_state.m_eps_p = plasticloc.GetState().m_eps_p;
        this->MemItem(intPt).m_elastoplastic_state.m_hardening = plasticloc.GetState().m_hardening;
        this->MemItem(intPt).m_plastic_steps = plasticloc.IntegrationSteps();
        this->MemItem(intPt).m_ER = plasticloc.GetElasticResponse();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).m_u[i] = data.sol[0][i];
        }
    }


}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::EigenValues(TPZFMatrix<REAL> & vectorTensor, TPZVec<REAL> & ev)
{
    TPZFNMatrix<9> Tensor(3,3);
    ev.Resize(3);
    this->vectorToTensor(vectorTensor, Tensor);
    int64_t numiterations = 1000;

    #ifdef PZDEBUG
    bool result = Tensor.SolveEigenvaluesJacobi(numiterations, m_tol, &ev);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << m_tol << std::endl;
        #ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::EigenValues *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
        #endif
    }
    #else
    Tensor.SolveEigenvaluesJacobi(numiterations, m_tol, &ev);
    #endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::EigenVectors(TPZFMatrix<REAL> &vectorTensor, TPZVec< REAL > &Solout, int direction)
{
    TPZFNMatrix<9> Tensor(3,3);
    this->vectorToTensor(vectorTensor, Tensor);

    TPZManVector<REAL,3> Eigenvalues(3);
    TPZFNMatrix<9> Eigenvectors(3,3);

    int64_t numiterations = 1000;
    #ifdef PZDEBUG
    bool result = Tensor.SolveEigensystemJacobi(numiterations, m_tol, Eigenvalues, Eigenvectors);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << m_tol << std::endl;
        #ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::EigenVectors *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
        #endif
    }
    #else
    Tensor.SolveEigensystemJacobi(numiterations, m_tol, Eigenvalues, Eigenvectors);
    #endif
    Solout.Resize(3);
    for(int i = 0; i < 3; i++) Solout[i] = Eigenvectors(direction,i);
}

template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlastic<T,TMEM>::NewMaterial() const
{
    return new TPZMatElastoPlastic<T,TMEM>(*this);
}
/*
 * void TPZMatElastoPlastic::SetData(std::istream &data)
 * {
 * TPZMaterial::SetData(data);
 * data >> fDeltaT; // to be removed in the elastoplastic material and readded to the poroelastoplastic material
 * }*/

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template <class T, class TMEM>
std::string TPZMatElastoPlastic<T,TMEM>::Name() const {
    return "TPZMatElastoPlastic<T,TMEM>";
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TPZMatWithMem<TMEM>::Write(buf, withclassid);

    buf.Write(&m_force[0], 3);
    buf.Write(&m_PostProcessDirection[0], 3);
    m_plasticity_model.Write(buf, withclassid);
    buf.Write(&m_tol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Read(TPZStream &buf, void *context) {
    //    TPZSavable::Read(buf, context);

    TPZMatWithMem<TMEM>::Read(buf, context);

    buf.Read(&m_force[0], 3);
    buf.Read(&m_PostProcessDirection[0], 3);
    m_plasticity_model.Read(buf, context);
    buf.Read(&m_tol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetTol(const REAL & tol)
{
    m_tol = tol;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetBulkDensity(const REAL & bulk)
{
    m_rho_bulk = bulk;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::vectorToTensor(const TPZFMatrix<REAL> & vectorTensor, TPZFMatrix<REAL> & Tensor)
{
    TPZTensor<REAL> vecT;
    vecT.CopyFrom(vectorTensor);
    vecT.CopyToTensor(Tensor);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::FillDataRequirements(TPZMaterialData &data) const
{

    TBase::FillDataRequirements(data);
    data.fNeedsSol = true;
    data.fNeedsNormal = false;
    data.fNeedsHSize = false;
    data.fNeedsNeighborCenter = false;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::FillBoundaryConditionDataRequirements(int type,TPZMaterialData &data) const
{
    TBase::FillDataRequirements(data);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}




