#include "TPZElasticMem.h"

// IMPORTANTE: TPZStream::Write/Read tem overloads específicos.
// Para PODs use buf.Write(double) e buf.Read(double&).
// Evite tentar buf.Write(&bool,1) etc.

void TPZElasticMem::Write(TPZStream &buf, int /*withclassid*/) const {
    // serializa apenas o que interessa: E e nu
    // Ajuste os getters conforme a API do TPZElasticResponse da sua branch.
    const REAL E  = m_ER.E();   // ou m_ER.YoungModulus()
    const REAL nu = m_ER.Poisson();  // ou m_ER.PoissonRatio()
    buf.Write(E);
    buf.Write(nu);
}

void TPZElasticMem::Read(TPZStream &buf, void* /*context*/) {
    REAL E = 0., nu = 0.;
    buf.Read(&E);
    buf.Read(&nu);
    // Ajuste o setter conforme sua API (SetEngineeringData / SetUp / Set(E,nu) etc.)
    m_ER.SetEngineeringData(E, nu);
}

void TPZElasticMem::Print(std::ostream &out) const {
    out << Name() << " { E=" << m_ER.E() << ", nu=" << m_ER.Poisson() << " }";
}

