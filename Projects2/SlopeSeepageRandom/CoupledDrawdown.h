// CoupledDrawdown.h
//
// Rebaixamento do nível d'água com o acoplamento hidromecânico (Biot) que o artigo não considera: lá o fluxo
// estacionário de Darcy é resolvido à parte e só entra como força de percolação -grad u no esqueleto. Aqui o
// problema u-p em deformação plana é resolvido no tempo com o material nativo TPZMatPoroElastoPlastic3DMem
// (dim = 2), com a mesma lei elastoplástica do esqueleto (Mohr-Coulomb ou Cam-Clay modificado):
//
//   cmeshU  H1 vetorial (ordem p) com memória nos pontos de integração; material TPZMatElastoPlastic2DSeepage<T>
//           (pós-processamento), que compartilha a memória com o material u-p;
//   cmeshP  H1 escalar de ordem 1 (TPZNullMaterial);
//   mphys   TPZMultiphysicsCompMesh com TPZMatPoroElastoPlastic3DMem<T>(ESoil, 2);
//   análise TPZLinearAnalysis (TPZSkylineNSymStructMatrix, LU) com o laço de Newton e Euler implícito.
//
// Tensões totais com p > 0 em compressão (poropressão total), força de corpo γ_sat, fluxo
// q = -(K/γw)(∇p + γw e_y) com K = diag(α k_v, k_v). Contorno: base fixa, laterais com u_x = 0 e sem fluxo;
// crista, face e pé com p = p_w = γw max(z_w(t) - y, 0) e a carga do reservatório t = -p_w n, onde
// z_w(t) = D - h_w min(t / t_d, 1) é o nível d'água (na crista em t = 0; rebaixado de h_w em t_d).
// Para t -> ∞ a poropressão tende à solução estacionária do artigo (mesmas condições de contorno).
//
// Estado inicial (nível na crista, poropressão hidrostática): Mohr-Coulomb — peso próprio aplicado em passos
// drenados; Cam-Clay — σ'0 da análise elástica geostática (γ'), em equilíbrio com p hidrostática.
//
// Tempo adimensional T = c_v t / H², c_v = k_v M / γw, M = E (1 - ν) / ((1 + ν)(1 - 2ν)) (módulo edométrico).
//
#ifndef COUPLEDDRAWDOWN_H
#define COUPLEDDRAWDOWN_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "Plasticity/TPZMatPoroElastoPlastic3DMem.h"
#include "SlopeStability.h"
#include "TPZMultiphysicsCompMesh.h"

template <class T>
class TCoupledDrawdown {
public:
    using TMatUP = TPZMatPoroElastoPlastic3DMem<T, TPZElastoPlasticMem>;
    using TPoint = typename TSlopeFEM<T>::TPoint;

    struct TParams {
        REAL hw = 5.;        ///< rebaixamento do nível d'água a partir da crista (m)
        REAL kv = 1.e-5;     ///< permeabilidade vertical (m/s); artigo: k_v/γw = 1e-6 m⁴/(kN s)
        REAL alpha = 1.;     ///< anisotropia k_h / k_v
        REAL Se = 0.;        ///< armazenamento 1/M_b (1/kPa); 0: fluido e grãos incompressíveis
        REAL Td = 0.1;       ///< duração adimensional do rebaixamento (T = c_v t / H²)
        int porderU = 2;     ///< ordem de u (p é linear)
        int nGravity = 5;    ///< passos drenados do peso próprio (Mohr-Coulomb)
        int nDrawdown = 10;  ///< passos durante o rebaixamento
        REAL growth = 1.5;   ///< fator de crescimento do passo após o rebaixamento
        REAL tol = 1.e-8;    ///< ||R_livre|| <= tol max(||F_ext||, 1)
        int maxIter = 25;
        int maxCuts = 6;     ///< cortes do passo de tempo antes de declarar colapso
        int verbose = 0;
    };

    TCoupledDrawdown(TPZGeoMesh *gmesh, const TSlopeGeometry &geo, const TSoil &soil, const TParams &par);
    ~TCoupledDrawdown();

    int64_t NEquations() const { return fMPhys->NEquations(); }
    /// Pontos de integração da malha de u (ordem dos índices de memória)
    const std::vector<TPoint> &Points() const { return fPoints; }
    /// c (kPa) e φ (rad) por ponto (antes de Initialize)
    void SetStrength(const std::vector<REAL> &c, const std::vector<REAL> &phi);
    /// Cam-Clay: σ'0 por ponto (antes de Initialize), em equilíbrio com γ' e p hidrostática
    void SetInitialStress(const std::vector<TPZTensor<REAL>> &sigma0);

    REAL Cv() const;                                    ///< coeficiente de adensamento (m²/s)
    REAL TimeScale() const;                             ///< H² / c_v (s)
    REAL Time() const { return fTime / TimeScale(); }   ///< tempo adimensional atual
    REAL WaterLevel() const { return fZw; }             ///< z_w atual (m)
    REAL WaterLevelAt(REAL Tad) const;                  ///< z_w no tempo adimensional Tad

    /// Equilíbrio inicial com o nível d'água na crista. Retorna false se não convergir.
    bool Initialize();
    /// Avança até o tempo adimensional T (passos de Euler implícito). false: colapso (Newton sem convergência
    /// mesmo com maxCuts cortes); o estado fica no último passo convergido.
    bool AdvanceTo(REAL Ttarget, REAL dtMax);

    /// Excesso de poropressão u = p - γw (D - y) do artigo (em relação à hidrostática com o nível na crista) e
    /// o seu gradiente (componentes globais) no ponto qsi do elemento geométrico gel
    void ExcessPorePressure(int64_t gel, const TPZVec<REAL> &qsi, REAL &u, TPZManVector<REAL, 2> &gradu) const;
    /// Transfere u e grad u (congelados) aos pontos de integração da análise de estabilidade
    template <class T2>
    void TransferSeepage(TSlopeFEM<T2> &fem) const {
        const auto &pts = fem.Points();
        std::vector<REAL> u(pts.size(), 0.);
        std::vector<TPZManVector<REAL, 2>> g(pts.size(), TPZManVector<REAL, 2>(2, 0.));
        for (size_t i = 0; i < pts.size(); i++)
            if (pts[i].gel >= 0) ExcessPorePressure(pts[i].gel, pts[i].qsi, u[i], g[i]);
        fem.SetSeepage(u, g);
    }

    /// Monitores: poropressão no ponto x, maior |Δu| desde o fim de Initialize, número de pontos plásticos
    REAL PorePressureAt(const TPZVec<REAL> &x) const;
    REAL MaxDisplacementIncrement() const;
    int64_t NPlasticPoints() const;

    void DefineVTK(const std::string &base);
    void WriteVTK(int step);

private:
    /// Um passo de Euler implícito (ou drenado com dt grande); false se não convergir (estado restaurado)
    bool Step(REAL dt, bool flow, int &iterations);
    void SetGravityFactor(REAL f);
    void InitializeMemory();

    TPZGeoMesh *fGMesh;
    TSlopeGeometry fGeo;
    TSoil fSoil;
    TParams fPar;
    TPZCompMesh *fCMeshU = nullptr, *fCMeshP = nullptr;
    TPZMultiphysicsCompMesh *fMPhys = nullptr;
    TMatUP *fMat = nullptr;
    TPZMatElastoPlastic2DSeepage<T> *fMatU = nullptr;
    std::unique_ptr<TPZLinearAnalysis> fAn;
    class TPZPostProcAnalysis *fPost = nullptr;
    std::vector<TPoint> fPoints;
    std::vector<TPZCompEl *> fCelP;  ///< elemento de p por índice de elemento geométrico
    std::vector<bool> fFree;
    std::vector<TPZTensor<REAL>> fSigma0;
    std::vector<REAL> fC, fPhi;
    std::vector<std::array<REAL, 2>> fU0;  ///< deslocamentos nos pontos no fim de Initialize
    REAL fTime = 0.;                 ///< tempo físico (s)
    REAL fDtAfter = 0.;              ///< passo atual após o rebaixamento (s)
    REAL fZw = 0.;                   ///< nível d'água atual
    REAL fGravity = 1.;              ///< fator do peso próprio e das cargas hidrostáticas
};

#endif
