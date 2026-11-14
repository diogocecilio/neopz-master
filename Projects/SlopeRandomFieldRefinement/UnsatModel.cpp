#include "UnsatModel.h"

#include "TPZGmshReader.h"
#include "TPZNullMaterial.h"
#include "pzgeoelbc.h"
#include "pzintel.h" // TPZInterpolatedElement (PRefine)
#include "pzerror.h"

#include <algorithm>
#include <array>
#include <functional>
#include <map>

// materiais das instancias do fim do arquivo
#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepModifiedCamClay.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "HPMaterial.h"

// A malha multifisica e dona do material e das BCs (apaga junto), mas nao das atomicas; a malha geo vai por ultimo,
// porque os elementos computacionais apontam para os geometricos
template <class TMat>
UnsatModel<TMat>::~UnsatModel()
{
    if (fCMesh) {
        TPZManVector<TPZCompMesh *, 3> atomicas = fCMesh->MeshVector(); // copia: depois do delete o vetor some
        delete fCMesh;
        for (TPZCompMesh *m : atomicas) delete m;
    } else {
        delete fMat; // dado, mas a malha comp nao foi criada: ainda e da classe
    }
    delete fGMesh;
}

template <class TMat>
void UnsatModel<TMat>::ReadGeoMesh(const std::string &file, TPZManVector<std::map<std::string, int>, 4> names)
{
    // TODO (passo A): a sua ReadSlopeMesh do main.cpp, com a tabela names vinda de fora
    // (TPZGmshReader::SetDimNamePhysical e GeometricGmshMesh; parar se algum elemento ficar com id sem nome)
    DebugStop(); // falta implementar
}

template <class TMat>
void UnsatModel<TMat>::DuplicateBoundary(int id, int newid)
{
    // TODO (passo A): TPZGeoElBC(gel, gel->NSides() - 1, newid) em cada elemento com o id. Os gemeos entram no fim
    // da lista de elementos: percorra so os que ja existiam
    DebugStop(); // falta implementar
}

template <class TMat>
std::set<int> UnsatModel<TMat>::BCIds() const
{
    std::set<int> ids;
    for (const DirichletBC &bc : fDirichlet) ids.insert(bc.id);
    for (const NaturalBC &bc : fNatural) ids.insert(bc.id);
    return ids;
}

// Malha de um espaco so (u, p_w ou p_a): o material nulo so diz quantas componentes (nstate) e onde ha elementos.
// A fisica fica toda no fMat, na malha multifisica
template <class TMat>
TPZCompMesh *UnsatModel<TMat>::CreateAtomicCompMesh(const UnsatField &field, const std::vector<int> &ordens)
{
    const int dim = fGMesh->Dimension();
    auto *cmesh = new TPZCompMesh(fGMesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(field.order);
    cmesh->SetAllCreateFunctionsContinuous(); // H1
    auto *mat = new TPZNullMaterial<STATE>(fMat->Id(), dim, field.nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(field.nstate, field.nstate, 0.);
    TPZManVector<STATE, 3> val2(field.nstate, 0.);
    for (int id : BCIds()) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2)); // so para o id existir
    cmesh->AutoBuild(); // um elemento comp para cada elemento geo com id conhecido
    if (!ordens.empty()) { // refinamento p: ordem propria nos elementos do dominio
        for (TPZCompEl *cel : cmesh->ElementVec()) {
            if (!cel || cel->Reference()->MaterialId() != fMat->Id()) continue;
            const int o = ordens[cel->Reference()->Index()];
            if (o != field.order) dynamic_cast<TPZInterpolatedElement *>(cel)->PRefine(o); // ajusta os lados vizinhos
        }
        cmesh->ExpandSolution();
    }
    return cmesh;
}

// Malha multifisica a partir do que o problema deu: fMat (dominio), fDirichlet e fNatural (contornos).
// Uma BC do material por id. Quem impoe Dirichlet e a analise (equacoes eliminadas), entao a BC de Dirichlet no
// material so marca o contorno; ja a natural (tracao, fluxo) e o material que integra, e ela tem prioridade no id
template <class TMat>
void UnsatModel<TMat>::CreateCompMesh(const std::vector<UnsatField> &fields, const std::vector<int> &ordens)
{
    if (!fGMesh || !fMat || fCMesh) DebugStop(); // falta a malha geo ou o material, ou a malha comp ja existe

    // 1. uma malha atomica por espaco; o refinamento p (ordens) so no deslocamento (espaco 0)
    const int nfields = fields.size();
    TPZManVector<TPZCompMesh *, 3> meshvec(nfields);
    for (int i = 0; i < nfields; i++) meshvec[i] = CreateAtomicCompMesh(fields[i], i == 0 ? ordens : std::vector<int>());

    // 2. malha multifisica com o material do dominio e as BCs
    fCMesh = new TPZMultiphysicsCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
    fCMesh->InsertMaterialObject(fMat);
    const int nstate = fMat->NStateVariables(); // componentes do deslocamento: val1 e val2 tem esse tamanho
    std::set<int> criados;
    for (const NaturalBC &bc : fNatural) {
        if (criados.count(bc.id)) DebugStop(); // duas BCs naturais no mesmo id: junte-as numa so
        TPZFMatrix<STATE> val1(nstate, nstate, 0.);
        TPZManVector<STATE, 3> val2(nstate, 0.);
        for (int i = 0; i < nstate && i < bc.val2.size(); i++) val2[i] = bc.val2[i];
        auto *bnd = fMat->CreateBC(fMat, bc.id, bc.type, val1, val2);
        if (bc.f) bnd->SetForcingFunctionBC(bc.f);
        fCMesh->InsertMaterialObject(bnd);
        criados.insert(bc.id);
    }
    for (const DirichletBC &bc : fDirichlet) {
        if (criados.count(bc.id)) continue; // o id ja tem BC no material (a primeira entrada ou uma natural)
        // marca do contorno com valor zero; o valor de verdade (bc.f) a analise le de fDirichlet. No u-p,
        // EDirichletUDirectional prende so as componentes com val1(c, c) = 1 (penalidade, se a analise nao eliminar)
        TPZFMatrix<STATE> val1(nstate, nstate, 0.);
        if (bc.field == 0)
            for (int c : bc.comps) val1(c, c) = 1.;
        TPZManVector<STATE, 3> val2(nstate, 0.);
        fCMesh->InsertMaterialObject(fMat->CreateBC(fMat, bc.id, bc.type, val1, val2));
        criados.insert(bc.id);
    }

    // 3. espacos ativos e elementos multifisicos; so o dominio tem memoria (pontos de integracao)
    TPZManVector<int, 3> active(nfields, 1);
    fCMesh->BuildMultiphysicsSpaceWithMemory(active, meshvec, {fMat->Id()}, criados);
}

// Refinamento h conforme por bisseccao pela aresta mais longa (Rivara), como o RefineMarked do Common:
//  - cada triangulo marcado e dividido duas vezes: o tamanho cai pela metade;
//  - um triangulo so e dividido pela sua aresta mais longa. Se o vizinho por essa aresta tem outra aresta mais longa,
//    o vizinho e dividido antes (recursivamente), ate a aresta ser a mais longa dos dois lados; ai os dois sao
//    divididos juntos, pelo mesmo ponto medio, e nao sobra no pendurado;
//  - os elementos de contorno (1D) de uma aresta dividida viram dois, com o mesmo id e a mesma orientacao.
// A malha nova e plana (sem hierarquia): primeiro os triangulos (os que nao mudaram, na ordem antiga, e depois os
// filhos, na ordem em que nasceram), depois os elementos de contorno, na ordem das arestas. Os nos antigos mantem o
// indice; os novos (pontos medios) vem no fim. pai[i] = elemento da malha atual que contem o elemento i da nova
template <class TMat>
TPZGeoMesh *UnsatModel<TMat>::RefinedGeoMesh(const std::set<int64_t> &marked, std::vector<int64_t> &pai) const
{
    using Edge = std::pair<int64_t, int64_t>; // nos em ordem crescente: a mesma chave para os dois lados da aresta
    auto chave = [](int64_t a, int64_t b) { return Edge(std::min(a, b), std::max(a, b)); };

    // 1. malha plana: os nos, as folhas de dimensao 2 (triangulos) e as de dimensao 1 (contorno)
    std::vector<TPZManVector<REAL, 3>> nos(fGMesh->NNodes(), TPZManVector<REAL, 3>(3, 0.));
    for (int64_t i = 0; i < fGMesh->NNodes(); i++) fGMesh->NodeVec()[i].GetCoordinates(nos[i]);
    std::vector<std::array<int64_t, 3>> tri;
    std::vector<int> triId, nbis; // id de material e bisseccoes que faltam
    std::vector<int64_t> triPai;  // elemento da malha atual de onde o triangulo veio
    std::vector<bool> vivo;
    struct Contorno {
        int id;
        int64_t pai, n0, n1; // n0 -> n1: a orientacao do elemento original
    };
    std::map<Edge, std::vector<Contorno>> contorno; // aresta -> elementos de contorno sobre ela (pode haver mais de um)
    for (TPZGeoEl *gel : fGMesh->ElementVec()) {
        if (!gel || gel->HasSubElement()) continue;
        if (gel->Dimension() == 2) {
            if (gel->Type() != ETriangle) DebugStop(); // a bisseccao e so para triangulos
            tri.push_back({gel->NodeIndex(0), gel->NodeIndex(1), gel->NodeIndex(2)});
            triId.push_back(gel->MaterialId());
            triPai.push_back(gel->Index());
            nbis.push_back(marked.count(gel->Index()) ? 2 : 0);
            vivo.push_back(true);
        } else if (gel->Dimension() == 1) {
            const int64_t a = gel->NodeIndex(0), b = gel->NodeIndex(1);
            contorno[chave(a, b)].push_back({gel->MaterialId(), gel->Index(), a, b});
        }
    }

    // 2. bisseccoes
    std::map<Edge, std::vector<int64_t>> edgeTri; // aresta -> triangulos vivos que a tem (1 no contorno, 2 dentro)
    auto aresta = [&](int64_t t, int i) { return chave(tri[t][i], tri[t][(i + 1) % 3]); };
    for (int64_t t = 0; t < (int64_t)tri.size(); t++)
        for (int i = 0; i < 3; i++) edgeTri[aresta(t, i)].push_back(t);
    auto novoTri = [&](const std::array<int64_t, 3> &v, int id, int64_t p, int n) {
        tri.push_back(v);
        triId.push_back(id);
        triPai.push_back(p);
        nbis.push_back(n);
        vivo.push_back(true);
        for (int i = 0; i < 3; i++) edgeTri[aresta(tri.size() - 1, i)].push_back(tri.size() - 1);
    };
    // aresta mais longa de t; no empate, a de maior chave (assim os dois triangulos de uma aresta escolhem igual)
    auto maisLonga = [&](int64_t t) {
        int best = 0;
        REAL l2best = -1.;
        for (int i = 0; i < 3; i++) {
            const Edge e = aresta(t, i);
            REAL l2 = 0.;
            for (int c = 0; c < 3; c++) l2 += (nos[e.second][c] - nos[e.first][c]) * (nos[e.second][c] - nos[e.first][c]);
            if (l2 > l2best || (l2 == l2best && e > aresta(t, best))) best = i, l2best = l2;
        }
        return best;
    };
    auto vizinho = [&](int64_t t, const Edge &e) -> int64_t { // o outro triangulo da aresta (-1 no contorno)
        for (int64_t s : edgeTri[e])
            if (s != t) return s;
        return -1;
    };
    // divide s pela aresta e no no mid: (a, mid, c) e (mid, b, c) mantem a orientacao; os filhos herdam id e pai
    auto divide = [&](int64_t s, const Edge &e, int64_t mid) {
        const std::array<int64_t, 3> v = tri[s];
        int i = 0;
        while (aresta(s, i) != e) i++;
        for (int j = 0; j < 3; j++) {
            auto &lista = edgeTri[aresta(s, j)];
            lista.erase(std::find(lista.begin(), lista.end(), s));
        }
        vivo[s] = false;
        const int n = std::max(nbis[s] - 1, 0);
        novoTri({v[i], mid, v[(i + 2) % 3]}, triId[s], triPai[s], n);
        novoTri({mid, v[(i + 1) % 3], v[(i + 2) % 3]}, triId[s], triPai[s], n);
    };
    std::function<void(int64_t)> bisecta = [&](int64_t t) {
        Edge e;
        int64_t n;
        for (;;) { // o vizinho pela aresta mais longa tem que ter essa mesma aresta como a mais longa
            e = aresta(t, maisLonga(t));
            n = vizinho(t, e);
            if (n < 0 || aresta(n, maisLonga(n)) == e) break;
            bisecta(n);
            if (!vivo[t]) return; // t ja foi dividido no caminho (os filhos herdaram as bisseccoes que faltavam)
        }
        const int64_t mid = nos.size();
        TPZManVector<REAL, 3> xm(3);
        for (int c = 0; c < 3; c++) xm[c] = 0.5 * (nos[e.first][c] + nos[e.second][c]);
        nos.push_back(xm);
        divide(t, e, mid);
        if (n >= 0) divide(n, e, mid);
        auto it = contorno.find(e);
        if (it != contorno.end()) { // aresta de contorno: cada elemento dela vira dois
            const std::vector<Contorno> lista = it->second;
            contorno.erase(it);
            for (const Contorno &c : lista) {
                contorno[chave(c.n0, mid)].push_back({c.id, c.pai, c.n0, mid});
                contorno[chave(mid, c.n1)].push_back({c.id, c.pai, mid, c.n1});
            }
        }
    };
    for (bool denovo = true; denovo;) { // ate nenhum triangulo vivo ter bisseccao pendente
        denovo = false;
        for (int64_t t = 0; t < (int64_t)tri.size(); t++)
            if (vivo[t] && nbis[t] > 0) bisecta(t), denovo = true;
    }

    // 3. malha nova: triangulos vivos e depois os contornos; o pai de cada elemento
    auto *g = new TPZGeoMesh;
    g->SetDimension(fGMesh->Dimension());
    g->NodeVec().Resize(nos.size());
    for (size_t i = 0; i < nos.size(); i++) g->NodeVec()[i].Initialize(nos[i], *g);
    pai.clear();
    int64_t index;
    for (size_t t = 0; t < tri.size(); t++) {
        if (!vivo[t]) continue;
        TPZManVector<int64_t, 3> v({tri[t][0], tri[t][1], tri[t][2]});
        g->CreateGeoElement(ETriangle, v, triId[t], index);
        pai.push_back(triPai[t]);
    }
    for (const auto &[e, lista] : contorno)
        for (const Contorno &c : lista) {
            TPZManVector<int64_t, 2> v({c.n0, c.n1});
            g->CreateGeoElement(EOned, v, c.id, index);
            pai.push_back(c.pai);
        }
    g->BuildConnectivity();
    return g;
}

// ---- materiais aceitos: uma linha por material. O neopz so compila estes materiais (fim do
// TPZMatPoroElastoPlasticUP.cpp e do TPZMatUnsaturatedUPwPa.cpp). Todas as linhas compilam o mesmo codigo: se uma
// funcao usar algo que so um material tem, ela nao compila aqui, e esse pedaco tem que vir do problema (do main)
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
template class UnsatModel<TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>>;                           // u-p, Mohr-Coulomb
template class UnsatModel<TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>>; // u-p, Cam Clay
template class UnsatModel<TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>>; // u-p_w-p_a, Mohr-Coulomb
// os mesmos com a regra de integracao do hp (HPMaterial.h)
template class UnsatModel<HPMaterial<TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>>>;
template class UnsatModel<HPMaterial<TPZMatPoroElastoPlasticUP<TPZPlasticStepModifiedCamClay, TPZElastoPlasticMem>>>;
template class UnsatModel<HPMaterial<TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>>>;
