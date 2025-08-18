#pragma once

#include "TPZStructMatrixT.h"
#include "TPZStrMatParInterface.h"
#include "pzcmesh.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
//#include "pzsave.h"
#include <ostream>
#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
template<class T>
class TPZElementMatrixT;
class TPZBaseMatrix;
class TPZStructMatrix;
/**
 * pzdoublestrmatriz<TVar>
 * - Deriva de TPZStructMatrixT<TVar> e implementa TPZStrMatParInterface.
 * - Gerencia montagem global de A=C (Nyström) e B=M (massa), consultando o estado do material TPZMatKLKernel.
 * - Não implementa física; só coordena a montagem (material fornece blocos locais).
 */
template<class TVar>
class pzdoublestrmatriz : public TPZStructMatrixT<TVar>, public virtual TPZStrMatParInterface
{
public:
    using Super = TPZStructMatrixT<TVar>;

    pzdoublestrmatriz(TPZCompMesh *mesh);
    pzdoublestrmatriz(TPZAutoPointer<TPZCompMesh> cmesh);

    // ---- TPZStructMatrixT ----
    TPZMatrix<TVar>* Create() override;

    // ---- TPZStrMatParInterface ----
    void Assemble(TPZBaseMatrix &stiffness, TPZBaseMatrix &rhs,
                  TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    void Assemble(TPZBaseMatrix &rhs,
                  TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    TPZBaseMatrix* CreateAssemble(TPZBaseMatrix &rhs,
                                  TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    // Clone
    TPZStructMatrix* Clone() override;

    // Config
    void SetMassOrder(int qmass) { fQMass = qmass; }

    // IO / RTTI
    int ClassId() const override;
    void Write(TPZStream &buf, int withclassid) const override;
    void Read(TPZStream &buf, void *context) override;

private:
    // Montagens globais (preenchem TPZFMatrix<TVar>)
    void AssembleC(TPZFMatrix<TVar> &C);
    void AssembleB(TPZFMatrix<TVar> &B);

    // Mapeia índices locais (shapes) -> equações globais
    void GetDestIndex(long iel, int nshape,
                      TPZManVector<long> &source, TPZManVector<long> &dest);

    // Lê o alvo (A/C ou B/M) do material TPZMatKLKernel
    bool CurrentTargetIsA() const;

private:
    int fQMass = 3; // ordem para massa consistente
};

// ==== Instanciações explicitas usuais ====
extern template class pzdoublestrmatriz<STATE>;
extern template class pzdoublestrmatriz<CSTATE>;
