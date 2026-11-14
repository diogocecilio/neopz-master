// Modelo de um problema poro-elastoplastico do neopz: malha geometrica, material do dominio, condicoes de contorno e
// malha computacional multifisica.
//
// Generico no material TMat, que ja traz a lei plastica TPlastic (Mohr-Coulomb, Modified Cam Clay, ...):
//   saturado u-p              TPZMatPoroElastoPlasticUP<TPlastic>   espacos u, p_w
//   nao saturado u-p_w-p_a    TPZMatUnsaturatedUPwPa<TPlastic>      espacos u, p_w, p_a
// O que e de um problema (ids da malha, dados do solo, pressoes do contorno) fica fora da classe, no main.
// Os materiais aceitos estao no fim do UnsatModel.cpp (uma linha por material).
#pragma once // o Common tem um SlopeModel.h com o guarda SLOPEMODEL_H; aqui nao ha guarda com nome

#include "TPZMultiphysicsCompMesh.h"
#include "TPZMatTypes.h" // ForcingFunctionBCType
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzquad.h" // TPZIntPoints
#include "pzerror.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Um espaco da formulacao: o deslocamento (nstate = 2) ou uma pressao (nstate = 1), com a sua ordem polinomial
struct UnsatField {
    int nstate = 1;
    int order = 1;
};

template <class TMat>
class UnsatModel {
public:
    // Condicao de Dirichlet, imposta pela analise (equacoes eliminadas): o espaco (0 = u, 1 = p_w, 2 = p_a), as
    // componentes prescritas e o valor f(x) de cada uma (sem f: zero)
    struct DirichletBC {
        int id, type, field;
        std::vector<int> comps;
        ForcingFunctionBCType<STATE> f;
    };
    // Condicao natural (tracao, fluxo): so o material contribui
    struct NaturalBC {
        int id, type;
        TPZManVector<STATE, 3> val2;
        ForcingFunctionBCType<STATE> f;
    };

    UnsatModel() = default;
    ~UnsatModel();
    // a classe e dona das malhas: uma copia apagaria as mesmas malhas duas vezes
    UnsatModel(const UnsatModel &) = delete;
    UnsatModel &operator=(const UnsatModel &) = delete;

    // ---- malha geometrica
    // Le o .msh do gmsh; names[dim][nome do grupo fisico] = id de material (dim 2: dominio, dim 1: contornos)
    void ReadGeoMesh(const std::string &file, TPZManVector<std::map<std::string, int>, 4> names);
    // Usa uma malha pronta (ex.: a de um nivel de refinamento); a classe passa a ser dona dela
    void SetGeoMesh(TPZGeoMesh *gmesh) { fGMesh = gmesh; }
    // Gemeo com o id newid sobre cada elemento do contorno id, para o contorno ter duas BCs (ex.: u_x = 0 e p).
    // Antes da CreateCompMesh
    void DuplicateBoundary(int id, int newid);

    // ---- material e condicoes de contorno (o problema decide quais)
    // Material do dominio, criado pelo problema com new (lei plastica e parametros); a malha comp vira dona dele
    void SetMaterial(TMat *mat) { fMat = mat; }
    // ex.: base engastada: AddDirichletBC(bcBase, TMat::EDirichletU, 0, {0, 1})
    void AddDirichletBC(int id, int type, int field, const std::vector<int> &comps,
                        const ForcingFunctionBCType<STATE> &f = {})
    {
        fDirichlet.push_back({id, type, field, comps, f});
    }
    // ex.: tracao: AddNaturalBC(bcTopo, TMat::ENeumannU, {0., -10.})
    void AddNaturalBC(int id, int type, const TPZVec<STATE> &val2, const ForcingFunctionBCType<STATE> &f = {})
    {
        fNatural.push_back({id, type, TPZManVector<STATE, 3>(val2), f});
    }

    // ---- malha computacional
    // fields: um item por espaco, na ordem do material: u-p {{2, 2}, {1, 1}}, u-p_w-p_a {{2, 2}, {1, 1}, {1, 1}}.
    // ordens: ordem de u de cada elemento, pelo indice geometrico (refinamento p); vazio: fields[0].order em todos
    void CreateCompMesh(const std::vector<UnsatField> &fields, const std::vector<int> &ordens = {});

    // ---- refinamento
    // Malha geometrica nova e conforme (a analise nao aceita nos pendurados), com os elementos marcados divididos e os
    // de contorno (inclusive os gemeos) acompanhando as arestas. pai[i] = elemento desta malha que contem o elemento i
    // da nova (para levar a ordem p, ou o que mais for preciso, ao proximo nivel)
    TPZGeoMesh *RefinedGeoMesh(const std::set<int64_t> &marked, std::vector<int64_t> &pai) const;

    // ---- pontos de integracao (depois da CreateCompMesh)
    // f(x, mem) em cada ponto de integracao do dominio: x = coordenadas do ponto, mem = o item de memoria dele no
    // material (TPZElastoPlasticMem no u-p, TPZUnsaturatedMem no nao saturado). Serve para dar a cada ponto as suas
    // propriedades (campos aleatorios): mem.m_elastoplastic_state.fmatprop = {c, phi, psi} (phi em rad) e, no nao
    // saturado, mem.fKMult (k do ponto = k do material x fKMult). A regra e a do elemento multifisico na montagem, e a
    // memoria segue a ordem dos pontos dela
    template <class F>
    void ForEachIntegrationPoint(F f)
    {
        for (TPZCompEl *cel : fCMesh->ElementVec()) {
            if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != fMat->Id()) continue;
            TPZGeoEl *gel = cel->Reference();
            auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mf) DebugStop();
            TPZManVector<int, 3> ordens(mf->NMeshes(), 0); // a maior ordem de cada espaco no elemento
            for (int i = 0; i < mf->NMeshes(); i++)
                if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
            const int o = fMat->IntegrationRuleOrder(ordens);
            std::unique_ptr<TPZIntPoints> regra(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
            TPZManVector<int, 3> ordem(gel->Dimension(), o);
            regra->SetOrder(ordem);
            TPZManVector<int64_t> mem;
            cel->GetMemoryIndices(mem);
            if (mem.size() != regra->NPoints()) DebugStop(); // a memoria nao corresponde a regra
            for (int ip = 0; ip < regra->NPoints(); ip++) {
                TPZManVector<REAL, 3> qsi(gel->Dimension()), x(3);
                REAL w;
                regra->Point(ip, qsi, w);
                gel->X(qsi, x);
                f(x, fMat->MemItem(mem[ip]));
            }
        }
    }

    // ---- acesso
    TPZGeoMesh *GeoMesh() { return fGMesh; }
    TPZMultiphysicsCompMesh *CompMesh() { return fCMesh; }
    TMat *Material() { return fMat; }
    const std::vector<DirichletBC> &DirichletBCs() const { return fDirichlet; }
    const std::vector<NaturalBC> &NaturalBCs() const { return fNatural; }

private:
    // Malha H1 de um espaco, com material nulo (TPZNullMaterial) no dominio e nos contornos, para a multifisica
    TPZCompMesh *CreateAtomicCompMesh(const UnsatField &field, const std::vector<int> &ordens);
    // Ids de contorno com alguma BC (de fDirichlet e de fNatural), sem repetir
    std::set<int> BCIds() const;

    TPZGeoMesh *fGMesh = nullptr;
    TPZMultiphysicsCompMesh *fCMesh = nullptr;
    TMat *fMat = nullptr; // depois da CreateCompMesh, mora dentro da fCMesh (ela apaga)
    std::vector<DirichletBC> fDirichlet;
    std::vector<NaturalBC> fNatural;
};
