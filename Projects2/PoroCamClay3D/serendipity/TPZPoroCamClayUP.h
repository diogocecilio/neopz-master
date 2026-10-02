// TPZPoroCamClayUP.h
//
// Elementos finitos acoplados u-p (Biot) 3D com o Cam-Clay modificado (TPZModifiedCamClay) como lei da
// tensão efetiva. Porte para o NeoPZ de fe3d_up.py:
//
//   EHex8   hexaedro de 8 nós   (u e p trilineares, Q1-Q1), integração 2x2x2
//   EHex20  hexaedro de 20 nós  (u serendipity quadrático, p trilinear nos 8 vértices: Q2-Q1), 3x3x3
//   EHex20R hexaedro de 20 nós com integração reduzida 2x2x2 (análogo do C3D20RP do Abaqus)
//
// A malha é um TPZGeoMesh do NeoPZ: elementos 3D TPZGeoCube (8 nós) ou TPZQuadraticCube (20 nós) e faces
// de contorno 2D (TPZGeoQuad / TPZQuadraticQuad) cujo material id é o marcador da face. As funções de
// forma de u são as do mapeamento geométrico (TPZCube / TPZQuadraticCube: o elemento é isoparamétrico),
// as de p são as trilineares dos vértices; as regras de integração são TPZIntCube3D / TPZIntQuad; o
// sistema linear (não simétrico) é resolvido com TPZSkylNSymMatrix, com os nós renumerados por
// TPZCutHillMcKee.
//
// Formulação:  K u - Q p = f_u ;  Qᵀ u' + S p' + H p = f_p, Euler implícito multiplicado por Δt:
//   [[K_T, -Q], [Qᵀ, S + Δt H]] {Δu, Δp} = -{R_u, R_p}  (Newton), com
//   R_u = F_int(σ') - Q P - F_ext,   R_p = Qᵀ(U - U_n) + S(P - P_n) + Δt (H P - f_g),
//   Q = α ∫ Bᵀ m N_p dΩ,  S = ∫ (1/M) N_pᵀ N_p dΩ,  H = ∫ k ∇N_pᵀ ∇N_p dΩ,  f_g = ∫ k ∇N_pᵀ (γ_w g⃗) dΩ.
// Convenções: tração positiva, p > 0 em compressão, σ = σ' - α p I; Voigt do Cam-Clay
// {xx, xy, xz, yy, yz, zz} com distorções de engenharia.
//
#ifndef TPZPOROCAMCLAYUP_H
#define TPZPOROCAMCLAYUP_H

#include <array>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "pzgmesh.h"
#include "pzfmatrix.h"
#include "TPZTensor.h"
#include "TPZModifiedCamClay.h"

class TPZPoroCamClayUP {
public:
    enum EElement { EHex8 = 0, EHex20 = 1, EHex20R = 2 };

    static std::string ElementName(EElement t) {
        return t == EHex8 ? "hex8" : (t == EHex20 ? "hex20" : "hex20r");
    }

    /// Parâmetros do acoplamento hidromecânico
    struct TCoupling {
        REAL alpha = 1.;                     ///< coeficiente de Biot
        REAL invBiotModulus = 0.;            ///< 1/M (0: fluido e grãos incompressíveis)
        REAL perm = 0.;                      ///< mobilidade k/μ (m²/(kPa s))
        REAL gammaW = 0.;                    ///< peso específico do fluido (termo gravitacional do fluxo)
        std::array<REAL, 3> body = {0., 0., 0.};  ///< força de corpo por volume (p.ex. (0, 0, -γ_sat))
    };

    /// Estado de um ponto de integração no fim do último passo convergido
    struct TGPState {
        TPZTensor<REAL> epsp, sig, eps;  ///< ε^p, σ' e ε total
        REAL alpha = 0.;
        bool plastic = false;            ///< plástico no último passo
    };

    /// Face de contorno (nós na ordem de TPZGeoQuad / TPZQuadraticQuad) e marcador (material id)
    struct TFace {
        std::vector<int64_t> nodes;
        int marker = 0;
    };

    using TParFunc = std::function<TPZModifiedCamClay(const TPZVec<REAL> &x)>;
    using TScalarFunc = std::function<REAL(const TPZVec<REAL> &x)>;
    using TFaceSelect = std::function<bool(const TFace &f)>;
    using TVecFunc = std::function<std::array<REAL, 3>(const TPZVec<REAL> &x)>;

    /// parOfX: parâmetros (e σ'0) do Cam-Clay no ponto x; p0OfX: poropressão inicial nos nós.
    /// initAtCentroid = false: estado inicial em cada ponto de Gauss (equilíbrio inicial exato);
    /// true: um estado por elemento, avaliado no centróide (como as zonas do FLAC3D).
    TPZPoroCamClayUP(TPZGeoMesh *gmesh, EElement type, const TParFunc &parOfX, const TCoupling &coupling,
                     const TScalarFunc &p0OfX, bool initAtCentroid = false);

    // ----------------------------------------------------------------------- acesso
    int64_t NNodes() const { return fNN; }
    int64_t NU() const { return 3 * fNN; }
    int64_t NP() const { return fNP; }
    int64_t NElements() const { return int64_t(fEl.size()); }
    EElement Type() const { return fType; }
    const std::array<REAL, 3> &Coord(int64_t node) const { return fX[node]; }
    int64_t PDof(int64_t node) const { return fPMap[node]; }
    const std::vector<int64_t> &ElementNodes(int64_t el) const { return fEl[el].nodes; }
    const std::vector<TFace> &Faces() const { return fFaces; }
    const std::vector<REAL> &U() const { return fU; }
    const std::vector<REAL> &P() const { return fP; }
    const std::vector<REAL> &Fb() const { return fFb; }
    const std::vector<std::vector<TGPState>> &State() const { return fState; }
    const TPZModifiedCamClay &GPModel(int64_t el, int ip) const { return fPar[fEl[el].par[ip]]; }
    int NGP(int64_t el) const { return int(fEl[el].gps.size()); }
    /// posição do ponto de integração
    const std::array<REAL, 3> &GPX(int64_t el, int ip) const { return fEl[el].gps[ip].x; }
    /// Forças internas F_int(σ') do último passo convergido
    const std::vector<REAL> &Fint() const { return fFint; }

    /// Nós que satisfazem um critério (p.ex. coordenadas)
    std::vector<int64_t> NodesWhere(const std::function<bool(const std::array<REAL, 3> &)> &crit) const;
    /// Nós das faces com um dado marcador
    std::vector<int64_t> NodesOnMarker(int marker) const;

    // ----------------------------------------------------------------------- cargas
    /// Forças nodais consistentes de uma tração uniforme t nas faces selecionadas
    std::vector<REAL> FaceLoad(const TFaceSelect &select, const std::array<REAL, 3> &t) const;
    /// Forças nodais de uma pressão (positiva comprimindo): F_a = -∫ N_a p n dA, n orientada por outward(x)
    std::vector<REAL> FacePressure(const TFaceSelect &select, REAL pressure, const TVecFunc &outward) const;

    // ----------------------------------------------------------------------- solução
    /// Forças internas e verificação do estado atual (sem alterar o estado)
    std::vector<REAL> InternalForces() const;

    /// Um passo de Euler implícito com Newton. fixedU: {gdl de u: valor total}, fixedP: {gdl de p: valor}.
    /// Uguess: chute inicial dos deslocamentos (o estado convergido continua sendo U_n, P_n do passo).
    /// Lança std::runtime_error (ou TPZModifiedCamClay::ReturnMappingError) se não convergir; nesse caso
    /// o modelo não é alterado. Retorna o número de iterações (1 = resíduo inicial já abaixo da tolerância).
    int Step(const std::vector<REAL> &Fext, const std::map<int64_t, REAL> &fixedU,
             const std::map<int64_t, REAL> &fixedP, REAL dt, bool flow, REAL tol = 1.e-8, int maxit = 30,
             const std::vector<REAL> *Uguess = nullptr, bool verbose = false);

    /// R = F_int(σ') - Q P - F_ext do último passo: zero nos gdl livres, reações nos restritos
    std::vector<REAL> Reactions(const std::vector<REAL> &Fext) const;

    // ----------------------------------------------------------------------- pós-processamento
    /// Elemento que contém x e as coordenadas paramétricas (Newton no mapeamento isoparamétrico)
    bool Locate(const std::array<REAL, 3> &x, int64_t &el, std::array<REAL, 3> &xi) const;
    /// σ' interpolada dos pontos de Gauss (polinômios de Lagrange nas abscissas de Gauss)
    TPZTensor<REAL> StressAt(int64_t el, const std::array<REAL, 3> &xi) const;
    /// Poropressão média dos 8 vértices do elemento (valor de "zona", como no FLAC3D)
    REAL ZonePressure(int64_t el) const;
    /// Poropressão interpolada no nó (vértices: gdl; nós de meio de aresta: média das extremidades)
    std::vector<REAL> NodalPressure() const;

    /// Arquivo VTK legado (UNSTRUCTURED_GRID, hexaedros lineares ou quadráticos de 20 nós) com
    /// deslocamentos, poropressão (e excesso em relação a hydro, se dado), tensões efetivas e totais
    /// (extrapoladas dos pontos de Gauss e médias nos nós) e, por elemento, médias de σ', p', q, α,
    /// p_c, fração de pontos plásticos/plastificados e poropressão de zona.
    void WriteVTK(const std::string &file, const std::string &title,
                  const std::vector<std::pair<std::string, REAL>> &fieldData = {},
                  const TScalarFunc &hydro = nullptr) const;

private:
    struct TGP {
        std::vector<REAL> Nu;            // funções de forma de u
        TPZFNMatrix<60, REAL> dNu;       // ∇N_u (3 x nu)
        REAL Np[8];                      // funções de forma de p
        REAL dNp[3][8];                  // ∇N_p
        REAL wdJ = 0.;
        std::array<REAL, 3> x, xi;
    };
    struct TElem {
        std::vector<int64_t> nodes;      // nu nós (ordem do NeoPZ)
        std::vector<TGP> gps;
        std::vector<int> par;            // índice em fPar por ponto de Gauss
        TPZFNMatrix<480, REAL> Qe;       // 3nu x 8
        TPZFNMatrix<64, REAL> Se, He;    // 8 x 8
    };

    void ShapeU(TPZVec<REAL> &xi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) const;
    void BuildEquations(const std::map<int64_t, REAL> &fixedU, const std::map<int64_t, REAL> &fixedP);
    /// Forças internas (e, se K != nullptr, a tangente consistente de cada elemento) para U
    void Assemble(const std::vector<REAL> &U, std::vector<REAL> &F, std::vector<std::vector<TGPState>> &trial,
                  std::vector<TPZFMatrix<REAL>> *Ke) const;
    /// Pesos de interpolação dos valores dos pontos de Gauss para o ponto paramétrico xi
    std::vector<REAL> GaussWeights(int64_t el, const std::array<REAL, 3> &xi) const;

    EElement fType;
    int fNu = 8;                          // nós de u por elemento
    int64_t fNN = 0, fNP = 0;
    std::vector<std::array<REAL, 3>> fX;
    std::vector<TElem> fEl;
    std::vector<TFace> fFaces;
    std::vector<int64_t> fPMap;           // nó -> gdl de p (-1: sem pressão)
    std::vector<int64_t> fPNode;          // gdl de p -> nó
    std::vector<TPZModifiedCamClay> fPar;
    TCoupling fCoupling;
    std::vector<REAL> fU, fP, fFg, fFb, fFint;
    std::vector<std::vector<TGPState>> fState;
    std::vector<REAL> fGauss1D;           // abscissas de Gauss (1D) da regra de volume

    // numeração das equações (gdl livres) e perfil do sistema
    std::vector<int64_t> fNodePerm;       // nó -> nova posição (Cuthill-McKee)
    std::vector<int64_t> fEqU, fEqP;      // gdl -> equação (-1: prescrito)
    std::vector<int64_t> fSkyline;
    int64_t fNEq = 0;
    std::vector<int64_t> fFixedKeyU, fFixedKeyP;
};

#endif
