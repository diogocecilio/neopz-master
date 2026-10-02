// PoroCamClayNativo.h
//
// Problema u-p 3D com Cam-Clay modificado montado com a estrutura nativa do NeoPZ:
//
//   cmeshU  H1 vetorial com memória (SetAllCreateFunctionsContinuousWithMem), material
//           TPZMatCamClayPostProc (TPZMatElastoPlastic<TPZModifiedCamClay>) que só é usado no
//           pós-processamento e compartilha a memória dos pontos de integração com o material u-p;
//   cmeshP  H1 escalar de ordem 1 (TPZNullMaterial);
//   mphys   TPZMultiphysicsCompMesh com TPZMatPoroElastoPlastic3DMem<TPZModifiedCamClay>;
//   an      TPZLinearAnalysis (renumeração de banda, TPZSkylineNSymStructMatrix, LU) com o laço de Newton;
//   post    TPZPostProcAnalysis sobre cmeshU (tensões e variáveis de estado guardadas nos pontos de
//           integração), como no SlopeAnalysis/GeoMecDeterm;
//   VTK     an.DefineGraphMesh / PostProcess (u, p, excesso de p, fluxo) e post.DefineGraphMesh /
//           TransferSolution / PostProcess (σ', σ, p', q, p_c, α, ...), e TPZVTKGeoMesh para a malha.
//
#ifndef POROCAMCLAYNATIVO_H
#define POROCAMCLAYNATIVO_H

#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLinearAnalysis.h"
#include "pzpostprocanalysis.h"
#include "TPZModifiedCamClay.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElastoPlasticMem.h"
#include "TPZMatPoroElastoPlastic3DMem.h"

using TMatUP = TPZMatPoroElastoPlastic3DMem<TPZModifiedCamClay, TPZElastoPlasticMem>;
using TModelFn = std::function<void(const TPZVec<REAL> &x, TPZModifiedCamClay &model)>;
using TScalarFn = std::function<REAL(const TPZVec<REAL> &x)>;

/// Material da malha de u usado no pós-processamento (TPZPostProcAnalysis): lê a memória dos pontos de
/// integração (compartilhada com o material u-p) e acrescenta variáveis do adensamento com Cam-Clay às do
/// TPZMatElastoPlastic (Displacement, VolHardening, StressXX, ...).
class TPZMatCamClayPostProc : public TPZMatElastoPlastic<TPZModifiedCamClay, TPZElastoPlasticMem> {
    using TBase = TPZMatElastoPlastic<TPZModifiedCamClay, TPZElastoPlasticMem>;

public:
    enum {
        EPorePressure = 80, EExcessPorePressure, EMeanEffectiveStress, EDeviatoricStress, EPreconsolidation,
        EPlasticPoint, EEffectiveStress, ETotalStress
    };
    explicit TPZMatCamClayPostProc(int id) : TBase(id) {}
    TPZMaterial *NewMaterial() const override { return new TPZMatCamClayPostProc(*this); }
    std::string Name() const override { return "TPZMatCamClayPostProc"; }

    void SetModelUpdate(const TModelFn &f) { fModelUpdate = f; }
    void SetHydrostatic(const TScalarFn &f) { fHydro = f; }
    void SetAlpha(REAL a) { fAlpha = a; }
    void SetIntegrationOrder(int order) { fIntegrationOrder = order; }

    int IntegrationRuleOrder(const int elPMaxOrder) const override;
    int VariableIndex(const std::string &name) const override;
    int NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;

private:
    TModelFn fModelUpdate;
    TScalarFn fHydro;
    REAL fAlpha = 1.;
    int fIntegrationOrder = -1;
};

class TPoroCamClayNativo {
public:
    /// Condição de contorno do material u-p (ver TPZMatPoroElastoPlastic3DMem.h): Val1 3x3, Val2 {v_x, v_y, v_z, p}
    struct TBC {
        int id;
        int type;
        std::array<REAL, 3> v1diag;
        std::array<REAL, 4> v2;
    };
    struct TParams {
        int orderU = 2;              ///< 1: Q1-Q1 (hexaedro de 8 nós); 2: u quadrático (Q2) e p linear
        int integrationOrder = -1;   ///< <= 0: 2 p (completa); > 0: ordem da regra (3 = 2x2x2: instável com Q2)
        TPZModifiedCamClay model;    ///< modelo de referência
        TModelFn modelUpdate;        ///< ajuste por ponto (σ'0, v0, ...); opcional
        TScalarFn p0;                ///< poropressão inicial; opcional
        TScalarFn hydro;             ///< poropressão hidrostática (excesso); opcional
        REAL alpha = 1., Se = 0., k = 0., mu = 1., rhof = 0.;
        std::array<REAL, 3> g = {0., 0., 0.}, body = {0., 0., 0.};
        std::vector<TBC> bcs;
        int vtkResolution = 0;
    };

    TPoroCamClayNativo(TPZGeoMesh *gmesh, const TParams &par);
    ~TPoroCamClayNativo();

    TPZMultiphysicsCompMesh *MPhys() { return fMPhys; }
    TPZCompMesh *MeshU() { return fCMeshU; }
    TPZCompMesh *MeshP() { return fCMeshP; }
    TMatUP *Material() { return fMat; }
    TPZLinearAnalysis &Analysis() { return *fAn; }
    int64_t NEquations() const { return fMPhys->NEquations(); }

    /// Altera Val2 de uma condição de contorno do material u-p
    void SetBCVal2(int id, const std::array<REAL, 4> &v2);

    /// Um passo de Euler implícito com Newton (ek Δx = ef). Converge quando ||R_livre|| < tol max(||F_ext||, 1)
    /// (como no fe3d_up.py); dirichletChanged força ao menos uma solução. guess: chute inicial (vetor do
    /// multifísico). Em caso de falha restaura o estado e lança std::runtime_error (ou ReturnMappingError).
    /// Se convergir, atualiza a memória (SetUpdateMem + AssembleResidual) e retorna o número de iterações.
    int Step(REAL dt, bool flow, bool dirichletChanged, const TPZFMatrix<STATE> *guess = nullptr,
             REAL tol = 1.e-8, int maxit = 30);

    /// Solução atual do multifísico (u e p totais) e máscara das equações de u
    const TPZFMatrix<STATE> &Solution() { return fAn->Solution(); }
    const std::vector<bool> &UEquationMask() const { return fIsU; }
    /// Escreve valores nas equações dos vértices (vetor do multifísico): u_c(nó) = f(x)
    void SetVertexU(TPZFMatrix<STATE> &vec, const std::function<std::array<REAL, 3>(const TPZVec<REAL> &)> &f) const;

    // ---------------------------------------------------------------- monitores
    REAL NodalU(int64_t node, int comp) const;
    REAL NodalP(int64_t node) const;
    /// Média da poropressão nos 8 vértices do elemento (valor de "zona" do FLAC3D)
    REAL ZonePressure(TPZGeoEl *gel) const;
    REAL MaxAbsPressure() const;
    /// Elemento 3D que contém x e as coordenadas paramétricas
    TPZGeoEl *FindVolumeElement(const TPZVec<REAL> &x, TPZVec<REAL> &qsi) const;
    /// σ' interpolada dos pontos de integração (polinômios de Lagrange nas abscissas de Gauss)
    TPZTensor<REAL> StressAt(TPZGeoEl *gel, const TPZVec<REAL> &qsi) const;
    /// Soma das reações (F_int - Q p - F_ext) na componente comp dos vértices dados
    REAL SumReactions(const std::vector<int64_t> &nodes, int comp);
    /// ||F_ext|| (forças de corpo e cargas de contorno)
    REAL ExternalForceNorm();
    /// max |R_u| nos graus de liberdade livres para a solução atual (verificação do equilíbrio inicial)
    REAL FreeResidualU();

    // ---------------------------------------------------------------- VTK
    /// Define os arquivos <base>_up.scal_vec.N.vtk (multifísico) e <base>_tensoes.scal_vec.N.vtk (memória)
    void DefineVTK(const std::string &base);
    void WriteVTK(int step);

private:
    void BuildVertexMaps();
    void InitializeMemory();
    int64_t MPhysEquation(int64_t uconnect, int comp) const;

    TParams fPar;
    TPZGeoMesh *fGMesh = nullptr;
    TPZCompMesh *fCMeshU = nullptr, *fCMeshP = nullptr;
    TPZMultiphysicsCompMesh *fMPhys = nullptr;
    TMatUP *fMat = nullptr;
    TPZMatCamClayPostProc *fMatU = nullptr;
    TPZLinearAnalysis *fAn = nullptr;
    TPZPostProcAnalysis *fPost = nullptr;
    std::vector<bool> fFree, fIsU;
    std::map<int64_t, std::pair<TPZCompEl *, int>> fVertexU, fVertexP;  // nó -> (elemento, vértice local)
    std::map<TPZGeoEl *, TPZCompEl *> fGelToU;
};

#endif
