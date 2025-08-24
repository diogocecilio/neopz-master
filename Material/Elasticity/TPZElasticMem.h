#pragma once
#include "TPZSavable.h"
#include "TPZElasticResponse.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"
#include "TPZHash.h" // Hash("...")

/**
 * Memória por ponto de integração para material elástico
 * guarda E e nu via TPZElasticResponse
 */
class TPZElasticMem : public TPZSavable {
public:
    TPZElasticResponse m_ER;

    // --- boilerplate de persistência ---
    TPZElasticMem() = default;
    TPZElasticMem(const TPZElasticMem&) = default;
    TPZElasticMem& operator=(const TPZElasticMem&) = default;
    ~TPZElasticMem() override = default;

    int ClassId() const override {
        return Hash("TPZElasticMem");
    }

    const std::string Name() const {
        return "TPZElasticMem";
    }

    void Write(TPZStream &buf, int withclassid) const override;
    void Read (TPZStream &buf, void *context) override;

    void Print(std::ostream &out = std::cout) const;

    friend std::ostream& operator<<(std::ostream& Out, const TPZElasticMem &s) {
        s.Print(Out);
        return Out;
    }
};

