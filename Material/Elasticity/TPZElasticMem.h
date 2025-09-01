//$Id: pzelasticmem.h, 2025-25-98 11:16:00 diogo Exp $

#ifndef PZELASTICMEM_H
#define PZELASTICMEM_H

#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZElastoPlasticMem.h"
#ifdef FIX_PLASTIC_TRANSLATORS
#include "TPZElastoPlasticMemTranslator.h"
#endif
#include "pzadmchunk.h"
#include "Plasticity/TPZElasticResponse.h"

/**
 * This class defines the material memory for a standar elastic calculation.
 */
class TPZElasticMem : public TPZSavable
{

private:

public:

    int ClassId() const override;

    TPZElasticMem();

    TPZElasticMem(const TPZElasticMem & other);

    const TPZElasticMem & operator = (const TPZElasticMem & other);

    virtual ~TPZElasticMem();

    const std::string Name() const;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;

    virtual void Print(std::ostream &out = std::cout) const;

    friend std::ostream& operator<<( std::ostream& Out, const TPZElasticMem & s )
    {
        s.Print(Out);
        return Out;
    }

    /// Elastoplastic response (It is required when elasti response depends on spatial variables)
    TPZElasticResponse m_ER;

};

#endif

