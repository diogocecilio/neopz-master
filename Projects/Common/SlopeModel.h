// Modelo do talude: dados, condicoes da agua, materiais (u-p e nao saturado) e malhas computacionais
#ifndef SLOPEMODEL_H
#define SLOPEMODEL_H

#include "SlopeMesh.h"

#include "TPZMultiphysicsCompMesh.h"
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>

namespace slope {

// Mohr-Coulomb (superficie PV2, associado) com resposta elastica linear
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMatUP = TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>; // Biot saturado u-p
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;      // nao saturado u-p_w-p_a

// Dados do talude (kN, m, s, kPa); padrao: o talude do relatorio (H = 5 m, 45 graus) com as medias da Tabela 2
struct Data {
    REAL c = 10.;                 // coesao efetiva (kPa)
    REAL phi = 30. * M_PI / 180.; // atrito efetivo (rad); fluxo associado
    REAL gamma = 20.;             // peso especifico do solo saturado (kN/m3)
    REAL E = 2.e4, nu = 0.3;      // elasticidade (nao muda a carga de colapso)
    REAL kv = 1.e-6;              // permeabilidade k_v (m4/kN/s)
    REAL gammaw = 9.81;           // peso especifico da agua (kN/m3)
    REAL H = 5.;                  // altura do talude (m)
    REAL beta = M_PI / 4.;        // inclinacao da face (rad)
    REAL xCrista = 10.;           // borda da crista (o ponto O do artigo)
    REAL yCrista = 10.;
    REAL hw = 5.;                 // rebaixamento h_w (m)
    REAL Nivel() const { return yCrista - hw; } // cota do nivel externo depois do rebaixamento
};

// Poropressao imposta no contorno id (Eq. 1 do artigo; laterais: hidrostatica da crista e do reservatorio)
REAL ContourPressure(int id, REAL y, const Data &d);
// Excesso de poropressao u (Eq. 2)
REAL ExcessPorePressure(REAL p, REAL y, const Data &d);
TMC MohrCoulomb(const Data &d);

// Material do talude: regra de integracao de ordem 2 (p_u - 1) (3 pontos em P2, 6 em P3, 12 em P4) e uma forca de corpo
// por ponto de integracao (forca de percolacao), vezes fForceFactor (o lambda do GIM). Nao altera os materiais base
template <class TMat>
class SlopeMat : public TMat {
public:
    using TMat::TMat;
    std::shared_ptr<std::vector<std::array<REAL, 2>>> fPointForce; // por item de memoria; nulo: sem forca
    REAL fForceFactor = 1.;

    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override
    {
        return std::max(2, 2 * (elPMaxOrder[0] - 1));
    }
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override
    {
        TMat::Contribute(datavec, weight, ek, ef);
        AddForce(datavec, weight, ef);
    }
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override
    {
        TMat::Contribute(datavec, weight, ef);
        AddForce(datavec, weight, ef);
    }

private:
    // ef = -R: a forca externa entra com sinal +, nas linhas de u (as primeiras, 2 por funcao)
    void AddForce(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) const
    {
        if (!fPointForce || fForceFactor == 0.) return;
        const TPZMaterialDataT<STATE> &du = datavec[0];
        const auto &f = (*fPointForce)[du.intGlobPtIndex];
        for (int a = 0; a < du.phi.Rows(); a++)
            for (int k = 0; k < 2; k++) ef(2 * a + k, 0) += weight * fForceFactor * f[k] * du.phi.GetVal(a, 0);
    }
};
using MatUP = SlopeMat<TMatUP>;
using MatUnsat = SlopeMat<TUnsat>;

// Malha atomica H1 com 'nstate' componentes; ordens: ordem de cada elemento do solo (p-refinamento), pelo indice geometrico
TPZCompMesh *CreateAtomicCompMesh(TPZGeoMesh *gmesh, int nstate, int ordem, const std::set<int> &bcids,
                                  const std::vector<int> &ordens = {});
// Malha u-p: base engastada, laterais com u_x = 0, poropressao do contorno (EDirichletP) na superficie e nos gemeos
TPZMultiphysicsCompMesh *CreateUPCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d, const std::vector<int> &ordens = {});
// Malha nao saturada (solo saturado em qualquer pressao); unitWeight: peso especifico (gamma ou o submerso)
TPZMultiphysicsCompMesh *CreateUnsatCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d, REAL unitWeight,
                                             const std::vector<int> &ordens = {});
void EraseMesh(TPZMultiphysicsCompMesh *mphys);

// Poropressao imposta em cada no da superficie e das laterais (a base fica sem fluxo)
std::map<int64_t, REAL> BoundaryPressures(TPZGeoMesh *gmesh, const Data &d);

// Pontos de integracao do solo: f(elemento, x, indice do item de memoria), na ordem da regra de cada elemento
int RuleOrder(TPZCompEl *cel);
void ForEachPoint(TPZMultiphysicsCompMesh *cmesh,
                  const std::function<void(TPZCompEl *, const TPZVec<REAL> &, int64_t)> &f);

REAL SqrtJ2Plastic(const TPZElastoPlasticMem &mem); // sqrt(J2) da deformacao plastica

// Maior sqrt(J2(eps_p)) dos pontos de cada elemento (o indicador do refinamento), pelo indice geometrico
template <class TMat>
std::vector<REAL> PlasticIndicator(TPZMultiphysicsCompMesh *cmesh, TMat *mat)
{
    std::vector<REAL> ind(cmesh->Reference()->NElements(), 0.);
    ForEachPoint(cmesh, [&](TPZCompEl *cel, const TPZVec<REAL> &, int64_t im) {
        REAL &v = ind[cel->Reference()->Index()];
        v = std::max(v, SqrtJ2Plastic(mat->MemItem(im)));
    });
    return ind;
}

// Maior f/(2c cos phi) de Mohr-Coulomb nas tensoes da memoria: 0 = sobre a superficie, > 0 = fora
template <class TMat>
REAL MaxYield(TMat *mat, const Data &d)
{
    auto &mem = *mat->GetMemory();
    REAL fmax = -1.e30;
    for (int64_t i = 0; i < mem.NElements(); i++) {
        const TPZTensor<REAL> &s = mem[i].m_sigma;
        const REAL c = 0.5 * (s.XX() + s.YY()), r = std::hypot(0.5 * (s.XX() - s.YY()), s.XY());
        const REAL smax = std::max({c + r, c - r, s.ZZ()}), smin = std::min({c + r, c - r, s.ZZ()});
        const REAL f = (smax - smin) + (smax + smin) * std::sin(d.phi) - 2. * d.c * std::cos(d.phi);
        fmax = std::max(fmax, f / (2. * d.c * std::cos(d.phi)));
    }
    return fmax;
}

// CSV com x, y e sqrt(J2(eps_p)) de cada ponto de integracao
template <class TMat>
void WritePlasticPoints(TPZMultiphysicsCompMesh *cmesh, TMat *mat, const std::string &file)
{
    std::ofstream out(file);
    out << "x,y,sqrtJ2\n";
    ForEachPoint(cmesh, [&](TPZCompEl *, const TPZVec<REAL> &x, int64_t im) {
        out << x[0] << "," << x[1] << "," << SqrtJ2Plastic(mat->MemItem(im)) << "\n";
    });
}

} // namespace slope

#endif
