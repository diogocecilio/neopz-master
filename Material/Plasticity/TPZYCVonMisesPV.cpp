#include "TPZYCVonMisesPV.h"
#include "pzerror.h"
#include "pzlog.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZStream.h"
#include <cmath>
#include <iomanip>



// --------- util local: obtém λ e μ a partir de E e ν ----------
static inline void LameFromER(const TPZElasticResponse &ER, REAL &lambda, REAL &mu)
{
    // Ajuste os nomes abaixo conforme sua versão de TPZElasticResponse
    const REAL E  = ER.E();          // ou ER.E()
    const REAL nu = ER.Poisson();   // ou ER.Nu()

    mu     = E/(2.0*(1.0 + nu));
    lambda = (E*nu)/((1.0 + nu)*(1.0 - 2.0*nu));
}

// ======================================================================
//                              CTORs / DTORs
// ======================================================================

TPZYCVonMisesPV::TPZYCVonMisesPV() :
fYieldStress(0.0)//, fER()
{
    // vazio
}

TPZYCVonMisesPV::TPZYCVonMisesPV(REAL yieldstress, TPZElasticResponse &ER) :
fYieldStress(yieldstress)//, fER(ER)
{
    // vazio
}

TPZYCVonMisesPV::TPZYCVonMisesPV(const TPZYCVonMisesPV &cp) :
TPZPlasticCriterion(cp),
fYieldStress(cp.fYieldStress)//,
//fER(cp.fER)
{
    // vazio
}

TPZYCVonMisesPV & TPZYCVonMisesPV::operator=(const TPZYCVonMisesPV &cp)
{
    if(this != &cp)
    {
        TPZPlasticCriterion::operator=(cp);
        fYieldStress = cp.fYieldStress;
        //fER = cp.fER;
    }
    return *this;
}

// ======================================================================
//                         Estado/local e parâmetros
// ======================================================================

void TPZYCVonMisesPV::SetUp(REAL yieldstress, TPZElasticResponse &ER)
{
    fYieldStress = yieldstress;
    //fER = ER;
}

void TPZYCVonMisesPV::SetLocalMatState(TPZPlasticState<REAL> & /*state*/)
{
    // Caso sua formulação precise guardar/alterar estado interno local,
    // implemente aqui. Por ora, sem estado interno específico.
     DebugStop(); // remova o DebugStop do header para compilar.
}

TPZPlasticState<REAL> TPZYCVonMisesPV::GetLocalMatState()
{
    // Retorne o estado interno local, se houver. Aqui devolvemos um default.
    TPZPlasticState<REAL> st;
    return st;
}


void TPZYCVonMisesPV::ChangeLocalMatParameters(TPZPlasticState<REAL> & /*state*/, REAL /*factor*/)
{
    // Atualize parâmetros locais (ex.: hardening), se necessário.
    // Mantido como stub intencional.
}



int TPZYCVonMisesPV::ClassId() const
{
    // Ajuste se sua infraestrutura exigir um ID fixo/Hash específico.
    // Retornar um valor estável é suficiente para muitos casos.
        return Hash("TPZYCVonMisesPV") ;
}

void TPZYCVonMisesPV::Read(TPZStream& buf, void* /*context*/)
{
    buf.Read(&fYieldStress,1);
    //fER.Read(buf,nullptr);
}

void TPZYCVonMisesPV::Write(TPZStream& buf, int /*withclassid*/) const
{
    buf.Write(&fYieldStress,1);
    //fER.Write(buf,0);
}

void TPZYCVonMisesPV::Print(std::ostream &out) const
{
    out << "----- TPZYCVonMisesPV -----\n";
    out << "Yield stress (sigma_y): " << std::setprecision(12) << fYieldStress << "\n";
    //out << "Elastic response (E, nu): E=" << fER.E() << "  nu=" << fER.Poisson() << "\n";
    out << "NYield = " << NYield << "\n";
    out << "---------------------------\n";
}

// ======================================================================
//                 Phi (função de escoamento) e helpers elásticos
// ======================================================================

void TPZYCVonMisesPV::Phi(TPZVec<STATE> sig_vec, STATE alpha, TPZVec<STATE> &phi) const
{
    // Von Mises: f = q - sqrt(2/3)*(sigma_y + alpha)
    // q = sqrt(3/2) ||s|| ; s = sigma - (tr(sigma)/3) * I, em principais
    if(phi.size() != as_integer(NYield)) phi.Resize(as_integer(NYield), 0.0);

    const STATE I1   = sig_vec[0] + sig_vec[1] + sig_vec[2];
    const STATE p    = I1 / 3.0;
    const STATE s0   = sig_vec[0] - p;
    const STATE s1   = sig_vec[1] - p;
    const STATE s2   = sig_vec[2] - p;
    const STATE sJ2  = 0.5*((s0 - s1)*(s0 - s1) + (s1 - s2)*(s1 - s2) + (s2 - s0)*(s2 - s0))/3.0 * 3.0;
    // A expressão acima é equivalente a J2 = 1/2 s:s; em principais pode-se usar:
    // J2 = ( (s0^2 + s1^2 + s2^2) )/2   (com s0+s1+s2=0). Para evitar confusão, reescrevemos explicitamente:
    const STATE J2   = 0.5*(s0*s0 + s1*s1 + s2*s2);
    const STATE q    = std::sqrt(3.0*J2); // q = sqrt(3*J2) = sqrt(3/2) ||s||

    const STATE sigY = fYieldStress + alpha; // alpha: variável escalar de encruamento (se aplicável)
    const STATE rhs  = std::sqrt(2.0/3.0) * sigY;

    phi[0] = q - rhs;
}

template<class T>
TPZVec<T> TPZYCVonMisesPV::SigmaElastPV(const TPZVec<T> &deform) const
{
    // σ_i = λ tr(ε) + 2 μ ε_i  (em principais)
    TPZVec<T> sigma(3, T(0));
    REAL lambda = 0.0, mu = 0.0;
    //LameFromER(fER, lambda, mu);

    const T tr = deform[0] + deform[1] + deform[2];
    const T common = T(lambda)*tr;
    sigma[0] = common + T(2.0*mu)*deform[0];
    sigma[1] = common + T(2.0*mu)*deform[1];
    sigma[2] = common + T(2.0*mu)*deform[2];
    return sigma;
}


void TPZYCVonMisesPV::ProjectSigma(const TPZTensor<STATE> & epst,const TPZTensor<STATE> & epsp,STATE k_prev)
{
    TPZFMatrix<REAL>  Ce,invCe;
    //fER.CMatrix(Ce) ;
    Ce.Print("sd");
    //fER.InvCMatrix(invCe) ;
    invCe.Print("invCe");
    // Usa o overload que calcula Dep (6x6) e reaproveita o resultado
    #ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> TPZYCVonMisesPV::ProjectSigma (trial)\n";
        // sout << "E=" << E << " nu=" << nu << " sigy=" << sigy << "\n";
        // sout << "sigma_trial_principal = " << sigma_trial_principal << "\n";
        // sout << "p=" << p << "  q=" << q << "  J2=" << J2t << "  f=" << fval << "\n";
        LOGPZ_DEBUG(loggerVonMIsesPV, sout.str().c_str());
    }
    #endif

    #ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << ">>> TPZYCVonMisesPV::ProjectSigma (proj)\n";
        // sout << "dgamma=" << dgamma << " denom=" << denom << "\n";
        // sout << "N=" << N << "\n";
        // sout << "sigma_proj_principal = " << sigma_principal_out << "\n";
        // sout << "Dep (6x6):\n";
        // Dep_out.Print(sout);
        LOGPZ_DEBUG(loggerVonMIsesPV, sout.str().c_str());
    }
    #endif
}


