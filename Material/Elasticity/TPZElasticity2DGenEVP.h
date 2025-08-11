#ifndef TPZELASTICITY2DGENEVP_H
#define TPZELASTICITY2DGENEVP_H

#include "Elasticity/TPZElasticity2D.h"
#include "TPZMatGeneralisedEigenVal.h"
//#include "pzsave.h" // para Hash e TPZSavable

/// Material para EVP generalizado K u = λ M u
/// Matriz A = rigidez (elasticidade), Matriz B = massa consistente
class TPZElasticity2DGenEVP :
    public TPZElasticity2D,
    public TPZMatGeneralisedEigenVal
{
public:
    enum EWhichMatrix { EMatrixA_K, EMatrixB_M };

    TPZElasticity2DGenEVP(int id, REAL E, REAL nu,
                          REAL fx, REAL fy,
                          REAL rho, REAL thickness = 1.);

    // --- Interface TPZMatGeneralisedEigenVal ---
    void SetMatrixA() override { fWhich = EMatrixA_K; }
    void SetMatrixB() override { fWhich = EMatrixB_M; }

    // --- Overrides para evitar conflito de múltipla herança ---
    std::string Name() const override { return "TPZElasticity2DGenEVP"; }
    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;

    // --- Montagem ---
    void Contribute(const TPZMaterialDataT<STATE> &data,
                    REAL weight,
                    TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZMaterialDataT<STATE> &data,
                      REAL weight,
                      TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    // --- variáveis para pós-processo (VTK) ---
    enum EPostVarEVP
    {
        EVP_UX = 0,
        EVP_UY,
        EVP_UMAG,
        EVP_UVEC,   // (opcional) deslocamento vetorial
        EVP_PORDER
    };

    int  VariableIndex(const std::string &name) const override;

    int  NSolutionVariables(int var) const override;

    // single-physics (este material é single)
    void Solution(const TPZMaterialDataT<STATE> &data,int var,TPZVec<STATE> &sol) override;

    // (opcional) compat forward para multiphysics, caso alguém chame com datavec
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,int var,TPZVec<STATE> &sol) ;

    void FillDataRequirements(TPZMaterialData &data);


private:
    /// Monta matriz de massa consistente
    void ContributeMass(const TPZMaterialDataT<STATE> &data,
                        REAL weight,
                        TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef);

    void ContributeBCMass(const TPZMaterialDataT<STATE> &data,
                                                 REAL weight,
                                                 TPZFMatrix<STATE> &ek,
                                                 TPZFMatrix<STATE> & /*ef*/,
                                                 TPZBndCondT<STATE> &bc);

    REAL fRho;            ///< densidade
    REAL fThick;          ///< espessura
    EWhichMatrix fWhich;  ///< qual matriz montar (K ou M)
};

#endif // TPZELASTICITY2DGENEVP_H

