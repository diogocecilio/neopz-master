// Material para o refinamento hp: a regra de integracao de cada elemento vem da ordem de u dele, ordem 2 (p_u - 1)
// (3 pontos em P2, 6 em P3, 12 em P4), como o SlopeMat do Common. Com a ordem fixa (SetIntegrationOrder), um elemento
// P3 seria integrado com a regra de um P2.
// Generico: TBase e qualquer material do UnsatModel (u-p ou nao saturado, com qualquer lei plastica), e o
// HPMaterial<TBase> entra no lugar dele como o TMat do UnsatModel e da UnsatAnalysis. Tudo o mais (BCs, memoria,
// PlasticModel(), SetGravityFactor...) e o do TBase
#pragma once

#include "pzvec.h"

#include <algorithm>

class TPZMaterial;

template <class TBase>
class HPMaterial : public TBase {
public:
    using TBase::TBase; // os mesmos construtores do material base

    // elPMaxOrder[0]: a maior ordem de u no elemento (o espaco 0). O numero de pontos de integracao (e de itens de
    // memoria) de cada elemento sai daqui quando a malha multifisica e montada
    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override
    {
        return std::max(2, 2 * (elPMaxOrder[0] - 1));
    }

    // copia continua sendo um HPMaterial (senao uma copia feita pelo neopz perderia a regra)
    TPZMaterial *NewMaterial() const override { return new HPMaterial(*this); }
};
