// Estabilidade do talude com refinamento adaptativo hp: em cada nivel monta as malhas, resolve a agua (percolacao do MEF
// ou forca de percolacao dada), calcula o fator (GIM ou SRM, por arc-length ou bisseccao) e refina em h e p a zona
// plastica. Ganchos para o futuro: chuva (transiente do material nao saturado, WaterModel::Rain) e campos aleatorios
// (propriedades por ponto de integracao, SetPointProperties)
#ifndef SLOPEANALYSIS_H
#define SLOPEANALYSIS_H

#include "SlopeModel.h"
#include "SlopeSolver.h"

#include <utility>

namespace slope {

enum class MaterialKind { UP, Unsat }; // TPZMatPoroElastoPlasticUP ou TPZMatUnsaturatedUPwPa
const char *Name(MaterialKind k);

// O que a agua faz no talude
enum class WaterModel {
    Dry,           // sem agua: peso total, pressoes nulas
    Drawdown,      // rebaixamento: percolacao permanente do MEF (esqueleto rigido), pressoes congeladas x lambda
    DrawdownForce, // a mesma percolacao como forca de corpo -grad u, com o peso submerso (deve dar o mesmo Drawdown)
    ForceField,    // forca de percolacao dada por SetSeepageForce (o campo semi-analitico), com o peso submerso
    Rain           // chuva (transiente do material nao saturado): a implementar
};
const char *Name(WaterModel w);

struct Refinement {
    int nunif = 1;   // refinamentos uniformes da malha inicial
    int nref = 3;    // niveis adaptativos
    int pmax = 3;    // ordem maxima de u
    REAL mark = 0.1; // marca os triangulos com sqrt(J2(eps_p)) >= mark x maximo
};

// Uma linha da tabela final
struct LevelSummary {
    int level = 0;
    int64_t ntri = 0, neq = 0, nfree = 0;
    std::map<int, int> orders; // ordem de u -> triangulos
    REAL value = 0.;           // Gamma (GIM) ou FS (SRM)
    int steps = 0;
    double seconds = 0.;
};

class SlopeAnalysis {
public:
    using ForceField = std::function<std::array<REAL, 2>(const TPZVec<REAL> &x)>; // forca de percolacao (kN/m3)
    // propriedades de um ponto pela posicao (campos aleatorios): {c, phi, psi} em mem.m_elastoplastic_state.fmatprop
    using PointProperties = std::function<void(const TPZVec<REAL> &x, TPZElastoPlasticMem &mem)>;

    SlopeAnalysis(const Data &d, MaterialKind kind, WaterModel water, FactorMethod factor);
    void SetRefinement(const Refinement &r) { fRef = r; }
    void SetSeepageForce(const ForceField &f) { fForce = f; }
    void SetPointProperties(const PointProperties &f) { fProps = f; }
    void SetOutputPrefix(const std::string &prefix) { fPrefix = prefix; } // arquivos por nivel; vazio: nenhum
    void SetTarget(REAL target) { fTarget = target; }                     // para o erro impresso
    // nivel 0: m0 com nunif refinamentos uniformes; depois nref niveis com h e p nos marcados
    std::vector<LevelSummary> Run(FlatMesh m0);

private:
    struct LevelResult {
        REAL value = 0., fmc = 0.;
        int64_t neq = 0, nfree = 0;
        int steps = 0;
        std::string end;
        std::vector<REAL> ind, pw; // indicador plastico por elemento; p_w por no (percolacao)
    };
    LevelResult SolveLevelUP(TPZGeoMesh *gmesh, const std::vector<int> &ordens, const std::string &header,
                             const std::string &base);
    LevelResult SolveLevelUnsat(TPZGeoMesh *gmesh, const std::vector<int> &ordens, const std::string &header,
                                const std::string &base);
    bool UsesForce() const { return fWater == WaterModel::DrawdownForce || fWater == WaterModel::ForceField; }

    Data fData;
    MaterialKind fKind;
    WaterModel fWater;
    FactorMethod fFactor;
    Refinement fRef;
    ForceField fForce;
    PointProperties fProps;
    std::string fPrefix;
    REAL fTarget = 0.;
};

// Uma rodada da tabela final: nome, alvo (para o erro) e os niveis
struct SummaryRun {
    std::string name;
    REAL target = 0.;
    std::vector<LevelSummary> levels;
};
void PrintSummary(const std::string &title, const std::vector<SummaryRun> &runs, double total);

} // namespace slope

#endif
