// SeepageProblem.h
//
// Problema hidráulico desacoplado do rebaixamento rápido (seção 3 do artigo): excesso de poropressão u = p - γw y'
// (y' = profundidade abaixo da crista) solução de div(-K grad u) = 0 em Ω, com
//   u = 0                 na crista,
//   u = -γw min(y', hw)   na face do talude,
//   u = -γw hw            na superfície do pé,
// e fluxo nulo na base e nos lados (eqs. 20-21). K = kv e_y⊗e_y + kh (1 - e_y⊗e_y), kh = α kv (eq. 17).
// As forças de percolação são -grad u (força de corpo, kN/m³).
//
// Implementação nativa: TPZDarcyFlow (H1) com a permeabilidade escalar substituída pelo tensor anisotrópico por
// elemento (TPZDarcyFlowAnisotropic, que só sobrescreve Contribute), condições de Dirichlet com
// TPZBndCondT::SetForcingFunctionBC, TPZLinearAnalysis com TPZSkylineStructMatrix e saída por
// DefineGraphMesh/PostProcess.
//
#ifndef SEEPAGEPROBLEM_H
#define SEEPAGEPROBLEM_H

#include <memory>
#include <string>
#include <vector>

#include "DarcyFlow/TPZDarcyFlow.h"
#include "SlopeGeometry.h"
#include "TPZLinearAnalysis.h"
#include "pzcmesh.h"

/// TPZDarcyFlow com K = kv(e) diag(α, 1) constante por elemento (kv relativo: o problema só depende das razões)
class TPZDarcyFlowAnisotropic : public TPZDarcyFlow {
public:
    TPZDarcyFlowAnisotropic(int id, int dim) : TPZDarcyFlow(id, dim) {}
    TPZMaterial *NewMaterial() const override { return new TPZDarcyFlowAnisotropic(*this); }
    std::string Name() const override { return "TPZDarcyFlowAnisotropic"; }

    void SetAnisotropy(REAL alpha) { fAlpha = alpha; }
    /// kv por identificador de elemento geométrico (vazio: kv = 1)
    void SetElementPermeability(const std::vector<REAL> &kvById) { fKv = kvById; }
    const std::vector<REAL> &ElementPermeability() const { return fKv; }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

private:
    REAL fAlpha = 1.;
    std::vector<REAL> fKv;
};

class TSeepageProblem {
public:
    struct TParams {
        REAL gammaW = 10.;  ///< peso específico da água (kN/m³)
        REAL hw = 5.;       ///< rebaixamento (m)
        REAL alpha = 1.;    ///< kh/kv
        int porder = 2;
    };

    TSeepageProblem(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TParams &par);
    ~TSeepageProblem();

    /// kv por elemento geométrico (índice do elemento; tamanho NElements do gmesh) — campo aleatório
    void SetElementPermeability(const std::vector<REAL> &kvByIndex);
    void Solve();

    /// u e grad u (componentes globais x, y) no ponto paramétrico qsi do elemento geométrico gel
    void Evaluate(int64_t gelIndex, const TPZVec<REAL> &qsi, REAL &u, TPZManVector<REAL, 2> &gradu);
    /// Funcional hidráulico J(u) = 1/2 ∫ grad u · K · grad u dΩ (eq. 22; v^d = 0) com K = kv diag(α, 1), kv
    /// relativo (1 sem campo aleatório): J/(k_h H² γw²) da Fig. 5 é Functional()/(α H² γw²)
    REAL Functional();

    void DefineVTK(const std::string &file);
    void WriteVTK(int step);
    TPZCompMesh *Mesh() { return fCMesh; }

private:
    TPZGeoMesh *fGMesh;
    TSlopeGeometry fGeo;
    TParams fPar;
    TPZCompMesh *fCMesh = nullptr;
    TPZDarcyFlowAnisotropic *fMat = nullptr;
    std::unique_ptr<TPZLinearAnalysis> fAn;
    std::vector<TPZCompEl *> fCelOfGel;  ///< elemento computacional (Darcy) de cada elemento geométrico
    bool fVTKDefined = false;
};

#endif
