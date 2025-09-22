//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "TPZElasticMem.h"
#ifdef FIX_PLASTIC_TRANSLATORS
#include "TPZElastoPlasticMemTranslator.h"
#endif
#include "pzadmchunk.h"

TPZElasticMem::TPZElasticMem()
{
}

TPZElasticMem::TPZElasticMem(const TPZElasticMem & other): m_ER(other.m_ER),fPorePressure(other.fPorePressure), fdPorePressure(other.fdPorePressure),fSolU(other.fSolU),fGradSolU(other.fGradSolU) {

}


TPZElasticMem::~TPZElasticMem(){

}

void TPZElasticMem::Write(TPZStream &buf, int withclassid) const
{
    m_ER.Write(buf, withclassid);
    buf.Write(&fPorePressure,1);
    buf.Write(&fdPorePressure[0],3);
    buf.Write(&fSolU[0],3);
    fGradSolU.Write(buf,withclassid);
}

void TPZElasticMem::Read(TPZStream &buf, void *context)
{
    m_ER.Read(buf, context);
    buf.Read(&fPorePressure,1);
    fdPorePressure.Resize(3);
    buf.Read(&fdPorePressure[0],3);
    buf.Read(&fSolU[0],3);
    fGradSolU.Read(buf,context);
}

void TPZElasticMem::Print(std::ostream &out)const
{
    out << Name();
    out << "\n fPorePressure = " << fPorePressure;
    out << "\n fdPorePressure = " << fdPorePressure;
    out << "\n fSolU = " << fSolU;
    out << "\n fGradSolU = " << fGradSolU;
    m_ER.Print(out);

}

const std::string TPZElasticMem::Name()const
{
    return "TPZElasticMem";
}

int TPZElasticMem::ClassId() const{
    return Hash("TPZElasticMem");
}

const TPZElasticMem & TPZElasticMem::operator=(const TPZElasticMem & other)
{

    /// check for self-assignment
    if(&other == this){
        return *this;
    }

    m_ER = other.m_ER;
    fPorePressure = other.fPorePressure;
    fdPorePressure = other.fdPorePressure;
    fSolU = other.fSolU;
    fGradSolU = other.fGradSolU;

    return *this;
}

#ifdef FIX_PLASTIC_TRANSLATORS
template class TPZRestoreClassWithTranslator<TPZElasticMem, TPZElastoPlasticMemTranslator>;
#endif
template class TPZRestoreClass<TPZAdmChunkVector<TPZElasticMem>>;

