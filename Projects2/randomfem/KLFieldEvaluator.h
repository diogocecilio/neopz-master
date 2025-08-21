// ================================
// File: KLFieldEvaluator.h  (versão em classe)
// ================================
#pragma once
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "TPZMatWithMem.h"
#include "pzmanvector.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include <random>
#include <cstdint>
#include <string>
#include <fstream>
#include <vector>

// Memória por QP: exemplo mínimo
struct TKLPointMem {
    STATE fE = 0.;
    int ClassId() const { return Hash("TKLPointMem"); }
    void Write(TPZStream& buf, int withclassid) const { buf.Write(&fE,1); }
    void Read (TPZStream& buf, void* ctx)             { buf.Read(&fE,1); }
};

// Manifesto simples para persistência de Alpha
struct KLBaseManifest {
    uint32_t version = 1;
    uint32_t M = 0;
    uint64_t ndof_base = 0;
    uint64_t gmesh_nodes = 0;
    uint64_t gmesh_elements = 0;
    uint64_t seed = 0;
};

// Mapeamento QP-alvo -> (gel_base, qsi_base)
struct KLQPMapEntry {
    long gel_id = -1;
    TPZManVector<REAL,3> qsi = TPZManVector<REAL,3>(3,0.);
};

class KLFieldEvaluator {
public:
    KLFieldEvaluator(const TPZCompMesh* base_cmesh,
                     const TPZFMatrix<STATE>& PHI_nodal,
                     const TPZVec<STATE>& lambdas);

    // Gera Alpha (M x S), com α_k = sqrt(λ_k) · ξ_k
    TPZFMatrix<STATE> DrawAlpha(int samples, uint64_t seed = 123456) const;

    // Persistência de Alpha (binário simples)
    static void SaveAlphaBIN(const std::string& path,
                             const TPZFMatrix<STATE>& Alpha,
                             const KLBaseManifest& man);
    static bool LoadAlphaBIN(const std::string& path,
                             TPZFMatrix<STATE>& Alpha,
                             KLBaseManifest& man);

    // Pré-avalia φ_k(x_qp) nos QPs da malha-alvo; também computa VarG por QP
    void BuildPhiTableAtTargetQPs(TPZCompMesh* target_cmesh,
                                  TPZMatWithMem<TKLPointMem>* target_mat,
                                  TPZFMatrix<float>& PhiQP_out,
                                  TPZVec<STATE>& varG_out) const;

    // Pré-mapeia QPs do alvo para coordenadas na base (para acelerar rebuilds)
    void BuildQPMap(TPZCompMesh* target_cmesh,
                    TPZMatWithMem<TKLPointMem>* target_mat,
                    std::vector<KLQPMapEntry>& out_map) const;

    // Usa um QPMap pronto para montar PhiQP e VarG
    void BuildPhiTableFromQPMap(const std::vector<KLQPMapEntry>& qpmap,
                                TPZFMatrix<float>& PhiQP_out,
                                TPZVec<STATE>& varG_out) const;

    // Aplica uma amostra s no MEM do material (gaussiano ou lognormal)
    void ApplySampleToMemory(const TPZFMatrix<float>& PhiQP,
                             const TPZVec<STATE>& varG,
                             const TPZFMatrix<STATE>& Alpha,
                             int sidx,
                             bool use_lognormal,
                             STATE mean_mu,
                             STATE cov_c,
                             TPZMatWithMem<TKLPointMem>* target_mat) const;

    // Acessos
    int M() const { return fM; }
    long NDofBase() const { return fPHI.Rows(); }

private:
    // Helpers internos
    static void ComputeShapeGlobalDOFIndices(TPZInterpolationSpace* intel,
                                             TPZCompMesh const* cmesh,
                                             TPZVec<long>& dofIdx);
    static void EvalPhiK_at_qsi_on_base(TPZInterpolationSpace* intel_base,
                                        TPZCompMesh const* base_cmesh,
                                        const TPZFMatrix<STATE>& PHI_nodal,
                                        const TPZManVector<REAL,3>& qsi,
                                        TPZVec<long> const& dofIdx,
                                        TPZVec<STATE>& phiK_out);

private:
    const TPZCompMesh* fBaseCmesh = nullptr;
    TPZFMatrix<STATE>  fPHI;      // [ndof_base x M]
    TPZVec<STATE>      fL;        // lambdas[0..M-1]
    int                fM = 0;
};

