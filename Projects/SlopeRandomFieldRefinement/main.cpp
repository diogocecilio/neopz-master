#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "TPZGmshReader.h"  // leitor de arquivos .msh do Gmsh
#include "TPZVTKGeoMesh.h"  // grava a malha geometrica em VTK
#include "pzgeoelbc.h"      // TPZGeoElBC: cria elementos de contorno (os gemeos das laterais)
#include "pzquad.h"         // TPZIntPoints: regras de integracao
#include "pzintel.h"         // TPZInterpolatedElement::PRefine
#include "pzmultiphysicselement.h" // ordens dos espacos de um elemento multifisico
#include "pzgmesh.h"        // TPZGeoMesh: a malha geometrica (nos + elementos)

#include "TPZSpStructMatrix.h"          // matriz esparsa (CSR), resolvida pelo PARDISO do MKL
#include "TPZUnsatAnalysis.h"           // analise do material nao saturado (Newton com equacoes eliminadas)
#include "pzstepsolver.h"               // solver direto (LU)
#include "Plasticity/TPZPoroElastoPlasticUPAnalysis.h" // analise do material u-p

#include "TPZBndCondT.h"             // condicoes de contorno
#include "TPZMultiphysicsCompMesh.h" // malha computacional com varios espacos (u, p_w, ...)
#include "TPZNullMaterial.h"         // material "vazio" das malhas atomicas
#include "pzcmesh.h"                 // TPZCompMesh: malha computacional de um espaco

#include "Plasticity/TPZElasticResponse.h"
#include "Plasticity/TPZMatPoroElastoPlasticUP.h"
#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"

// Lei constitutiva do esqueleto: Mohr-Coulomb (superficie PV2) com resposta elastica linear
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
// Material 1: Biot saturado u-p (deslocamento e poropressao)
using TMatUP = TPZMatPoroElastoPlasticUP<TMC, TPZElastoPlasticMem>;
// Material 2: nao saturado u-p_w-p_a (deslocamento, pressao da agua e pressao do ar)
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;

// Material com a regra de integracao de ordem 2 (p_u - 1), p_u a ordem de u do elemento: 3 pontos em P2 (como antes),
// 6 em P3, 12 em P4 (a base usaria 2 p_u). Classe a parte para nao alterar os materiais
template <class TMat>
class HPRule : public TMat {
public:
    using TMat::TMat;
    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override
    {
        return std::max(2, 2 * (elPMaxOrder[0] - 1));
    }
};

// Ids de material: o numero que o NeoPZ associa a cada grupo fisico do .geo.
// Dominio com id positivo, contornos com id negativo (como no SlopeRandomField).
enum EIds {
    matSolo = 1,       // Physical Surface("Solo")      : o solo (2D)
    bcBase = -1,       // Physical Curve("Base")        : y = 0
    bcDireita = -2,    // Physical Curve("Direita")     : x = 25
    bcSuperficie = -3, // Physical Curve("Superficie")  : terreno do pe, face e crista
    bcEsquerda = -4,   // Physical Curve("Esquerda")    : x = 0
    bcDireitaP = -12,  // gemeos das laterais: poropressao da malha u-p (DuplicateBoundaryBC)
    bcEsquerdaP = -14
};

// Dados do caso deterministico: as medias da Tabela 2 de Vargas Ceron et al. (2025). Unidades kN, m, s, kPa.
struct Data {
    REAL c = 10.;                 // coesao efetiva c' (kPa)
    REAL phi = 30. * M_PI / 180.; // angulo de atrito efetivo phi' (rad); fluxo associado, psi = phi
    REAL gamma = 20.;             // peso especifico do solo (kN/m3)
    REAL E = 2.e4, nu = 0.3;      // elasticidade: o artigo e analise limite (rigido-plastico), estes sao os do relatorio
    REAL kv = 1.e-6;              // permeabilidade k_v do artigo (m4/kN/s): e a mobilidade k/gamma_w
    REAL gammaw = 9.81;           // peso especifico da agua (kN/m3)
    // geometria (a da malha g15) e rebaixamento
    REAL H = 5.;                  // altura do talude (m)
    REAL yCrista = 10.;           // cota da crista na malha (base em y = 0, pe em y = 5)
    REAL hw = 5.;                 // rebaixamento h_w (m); h_w < H exige a carga da agua na superficie (falta)
    REAL Nivel() const { return yCrista - hw; } // cota do nivel d'agua externo depois do rebaixamento
};

// Poropressao imposta no contorno id (Eq. 1 do artigo; laterais como no relatorio)
REAL ContourPressure(int id, REAL y, const Data &d)
{
    switch (id) {
        case bcSuperficie: return d.gammaw * std::max(REAL(0), d.Nivel() - y);
        case bcEsquerda:
        case bcEsquerdaP: return d.gammaw * (d.yCrista - y);
        case bcDireita:
        case bcDireitaP: return d.gammaw * (d.Nivel() - y);
        default: DebugStop();
    }
    return 0.;
}

// Excesso de poropressao u (Eq. 2 do artigo) de uma pressao p na cota y
REAL ExcessPorePressure(REAL p, REAL y, const Data &d) { return p - d.gammaw * (d.yCrista - y); }

TPZGeoMesh *ReadGMesh(const std::string &file)
{
    // names[d] : nome do grupo fisico de dimensao d no .msh -> id de material no NeoPZ
    TPZManVector<std::map<std::string, int>, 4> names(4);
    names[2]["Solo"] = matSolo;
    names[1]["Base"] = bcBase;
    names[1]["Direita"] = bcDireita;
    names[1]["Superficie"] = bcSuperficie;
    names[1]["Esquerda"] = bcEsquerda;

    TPZGmshReader reader;
    reader.SetDimNamePhysical(names);
    TPZGeoMesh *gmesh = reader.GeometricGmshMesh(file);

    // nome que nao esta no .msh: o elemento fica com o numero do Gmsh e o contorno some sem aviso
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        const int id = gel ? gel->MaterialId() : matSolo;
        if (id != matSolo && id != bcBase && id != bcDireita && id != bcSuperficie && id != bcEsquerda) {
            std::cerr << "Id " << id << " sem nome correspondente em ReadGMesh\n";
            DebugStop();
        }
    }

    // vizinhanca entre elementos (lados compartilhados); necessaria antes de qualquer calculo na malha
    gmesh->BuildConnectivity();
    return gmesh;
}

// Imprime o que foi lido: nos, elementos por dimensao e material, e a caixa que contem os nos.
void GMeshSummary(TPZGeoMesh *gmesh)
{
    std::cout << "Dimensao da malha: " << gmesh->Dimension() << "\n";
    std::cout << "Nos: " << gmesh->NNodes() << "\n";

    std::map<int, int> porMaterial;           // id de material -> quantos elementos
    std::map<int, int> porTipo;               // numero de nos do elemento -> quantos (2 = segmento, 3 = triangulo)
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel) continue;
        porMaterial[gel->MaterialId()]++;
        porTipo[gel->NNodes()]++;
    }
    std::cout << "Elementos: " << gmesh->NElements() << "\n";
    for (auto &[nnos, n] : porTipo) std::cout << "  com " << nnos << " nos: " << n << "\n";
    for (auto &[id, n] : porMaterial) std::cout << "  material " << id << ": " << n << "\n";

    // caixa dos nos (no relatorio: 0 a 25 m por 0 a 10 m, o retangulo do Karhunen-Loeve)
    REAL xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    for (int64_t i = 0; i < gmesh->NNodes(); i++) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[i].GetCoordinates(x);
        xmin = std::min(xmin, x[0]); xmax = std::max(xmax, x[0]);
        ymin = std::min(ymin, x[1]); ymax = std::max(ymax, x[1]);
    }
    std::cout << "Caixa dos nos: x de " << xmin << " a " << xmax << ", y de " << ymin << " a " << ymax << "\n";
}


// Gemeos das laterais (ids -12 e -14): o original leva u_x = 0 e o gemeo a poropressao. Chamar antes das malhas.
void DuplicateBoundaryBC(TPZGeoMesh *gmesh)
{
    const int64_t nel = gmesh->NElements(); // os gemeos entram no fim do vetor e nao sao revisitados
    for (int64_t i = 0; i < nel; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == bcDireita || id == bcEsquerda)
            TPZGeoElBC(gel, gel->NSides() - 1, id - 10); // lado NSides()-1: o proprio segmento
    }
}

// Mohr-Coulomb associado (psi = phi) com os dados medios
TMC MohrCoulomb(const Data &d)
{
    TPZElasticResponse ER;
    ER.SetEngineeringData(d.E, d.nu);
    return TMC(TPZYCMohrCoulombPV2(d.phi, d.phi, d.c, ER), ER);
}

// Malha atomica: espaco H1 continuo de ordem 'ordem' com 'nstate' componentes (2 para u, 1 para uma pressao).
// O material nulo e as condicoes de contorno so marcam onde ha elementos; a fisica fica na malha multifisica.
// ordens (opcional): ordem de cada elemento do solo, pelo indice do elemento geometrico (p-refinamento)
TPZCompMesh *CreateAtomicCompMesh(TPZGeoMesh *gmesh, int nstate, int ordem, const std::set<int> &bcids,
                                  const std::vector<int> &ordens = {})
{
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(ordem);
    cmesh->SetAllCreateFunctionsContinuous(); // H1
    auto *mat = new TPZNullMaterial<STATE>(matSolo, 2, nstate);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild(); // cria um elemento computacional para cada elemento geometrico com material inserido
    if (!ordens.empty()) {
        for (TPZCompEl *cel : cmesh->ElementVec()) {
            if (!cel || cel->Reference()->MaterialId() != matSolo) continue;
            const int o = ordens[cel->Reference()->Index()];
            if (o != ordem) dynamic_cast<TPZInterpolatedElement *>(cel)->PRefine(o); // ajusta os lados com os vizinhos
        }
        cmesh->ExpandSolution();
    }
    return cmesh;
}

// Malha u-p (Biot saturado): os Dirichlet sao as condicoes de contorno do proprio material
TPZMultiphysicsCompMesh *CreateUPCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d, const std::vector<int> &ordens = {})
{
    using B = TPZMatPoroElastoPlasticUPBase;
    const std::set<int> bcids = {bcBase, bcDireita, bcSuperficie, bcEsquerda, bcDireitaP, bcEsquerdaP};
    TPZManVector<TPZCompMesh *, 2> meshvec = {CreateAtomicCompMesh(gmesh, 2, porder, bcids, ordens), // [0] u
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids)}; // [1] p_w

    auto *mat = new HPRule<TMatUP>(matSolo, B::EPlaneStrain);
    mat->SetPlasticModel(MohrCoulomb(d));
    mat->SetBiot(1., 0.);       // alpha_B = 1 e 1/M = 0: graos e agua incompressiveis
    mat->SetPermeability(d.kv); // mobilidade k/gamma_w (m4/kN/s), a mesma grandeza do k_v do artigo
    mat->SetBodyForce(TPZManVector<REAL, 3>({0., -d.gamma, 0.}));    // peso do solo saturado
    mat->SetFluidWeight(TPZManVector<REAL, 3>({0., -d.gammaw, 0.})); // peso da agua (gradiente hidrostatico)

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcBase, B::EDirichletU, val1, val2)); // u_x = u_y = 0
    TPZFMatrix<STATE> sel(2, 2, 0.);
    sel(0, 0) = 1.; // so a componente x e prescrita
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcDireita, B::EDirichletUDirectional, sel, val2));  // u_x = 0
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcEsquerda, B::EDirichletUDirectional, sel, val2)); // u_x = 0

    // poropressao: a funcao do contorno da p em cada ponto (a tracao nula da superficie e natural)
    for (int id : {bcSuperficie, bcDireitaP, bcEsquerdaP}) {
        auto *bc = mat->CreateBC(mat, id, B::EDirichletP, val1, val2);
        bc->SetForcingFunctionBC([id, d](const TPZVec<REAL> &x, TPZVec<STATE> &v2, TPZFMatrix<STATE> &) {
            v2[0] = ContourPressure(id, x[1], d);
        });
        mphys->InsertMaterialObject(bc);
    }

    TPZManVector<int, 2> active(2, 1); // os dois espacos sao incognitas
    // memoria nos pontos de integracao so do solo; os contornos nao tem memoria
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matSolo}, bcids);
    return mphys;
}

// Malha nao saturada: u e as pressoes do contorno sao impostos pelo TPZUnsatAnalysis (mapa equacao -> valor)
TPZMultiphysicsCompMesh *CreateUnsatCompMesh(TPZGeoMesh *gmesh, int porder, const Data &d,
                                             const std::vector<int> &ordens = {})
{
    const std::set<int> bcids = {bcBase, bcDireita, bcSuperficie, bcEsquerda};
    TPZManVector<TPZCompMesh *, 3> meshvec = {CreateAtomicCompMesh(gmesh, 2, porder, bcids, ordens), // [0] u
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids),  // [1] p_w
                                              CreateAtomicCompMesh(gmesh, 1, porder - 1, bcids)}; // [2] p_a

    auto *mat = new HPRule<TUnsat>(matSolo);
    mat->SetPlasticModel(MohrCoulomb(d));
    // solo saturado em qualquer pressao: Gardner com a_G = 0 da S_w = 1 e k_rw = 1, tambem com p_w < 0
    TPZUnsaturatedHydraulics h;
    h.fRetention = TPZUnsaturatedHydraulics::EGardner;
    h.fAG = 0.;
    h.fSwr = 0.;
    mat->SetHydraulics(h);
    mat->SetBishop(TUnsat::EChiConstant, 1.); // chi = 1: tensao efetiva de Terzaghi
    mat->SetPorosity(0.45);
    mat->SetBiot(1., 0.);                    // alpha = 1, graos incompressiveis
    mat->SetWaterBulkModulus(1.e30);         // agua incompressivel: sem armazenamento, percolacao permanente
    mat->SetConductivity(d.kv * d.gammaw, d.gammaw, 1.e-3 / 1.8e-5); // k_s = k_v gamma_w (m/s); mu_w / mu_a
    mat->SetAir(1.2 * 9.81e-3, 101.325);     // peso especifico do ar (kN/m3) e pressao atmosferica (kPa)
    mat->SetUnitWeight(d.gamma);
    mat->SetLumpedStorage(true);
    mat->SetHydroMechanicalCoupling(false);  // percolacao sem acoplamento (o esqueleto nao muda p_w), como o Ceron

    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mphys->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 2> val2(2, 0.);
    for (int id : {bcBase, bcDireita, bcEsquerda})
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, TUnsat::EDirichletU, val1, val2)); // so marca o contorno
    mphys->InsertMaterialObject(mat->CreateBC(mat, bcSuperficie, TUnsat::ENeumannU, val1, val2)); // tracao nula

    TPZManVector<int, 3> active(3, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, {matSolo}, bcids);
    return mphys;
}

// Equacoes de cada malha atomica e da multifisica, elementos e pontos de integracao (itens de memoria do solo).
template <class TMat>
void CMeshSummary(TPZMultiphysicsCompMesh *mphys, const std::vector<std::string> &espacos)
{
    TPZVec<TPZCompMesh *> &meshvec = mphys->MeshVector();
    for (int i = 0; i < meshvec.size(); i++)
        std::cout << "  " << espacos[i] << ": " << meshvec[i]->NEquations() << " equacoes (ordem "
        << meshvec[i]->GetDefaultOrder() << ")\n";
    std::cout << "  multifisica: " << mphys->NEquations() << " equacoes, " << mphys->NElements() << " elementos\n";
    auto *mat = dynamic_cast<TMat *>(mphys->FindMaterial(matSolo));
    if (!mat) DebugStop();
    std::cout << "  pontos de integracao com memoria: " << mat->GetMemory()->NElements() << "\n";
}

// Apaga a malha multifisica e as suas malhas atomicas
void EraseMesh(TPZMultiphysicsCompMesh *mphys)
{
    TPZVec<TPZCompMesh *> meshvec = mphys->MeshVector();
    delete mphys;
    for (TPZCompMesh *m : meshvec) delete m;
}

// ===================================================================================================================
// Passo 3: caso seco, Gamma pelo GIM com arc-length e refinamento adaptativo
// ===================================================================================================================

// Malha sem hierarquia: cada nivel vira uma TPZGeoMesh nova e conforme (as analises nao aceitam nos pendurados)
struct FlatMesh {
    std::vector<std::array<REAL, 2>> nodes;
    std::vector<std::array<int64_t, 3>> tri;       // triangulos do solo, na ordem dos elementos da TPZGeoMesh
    std::map<std::pair<int64_t, int64_t>, int> bc; // aresta de contorno (nos em ordem crescente) -> id
    std::vector<int> order;                        // ordem de u de cada triangulo
};

std::pair<int64_t, int64_t> EdgeKey(int64_t a, int64_t b) { return {std::min(a, b), std::max(a, b)}; }

FlatMesh FromGeoMesh(TPZGeoMesh *gmesh)
{
    FlatMesh m;
    for (int64_t i = 0; i < gmesh->NNodes(); i++)
        m.nodes.push_back({gmesh->NodeVec()[i].Coord(0), gmesh->NodeVec()[i].Coord(1)});
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == matSolo) {
            m.tri.push_back({gel->NodeIndex(0), gel->NodeIndex(1), gel->NodeIndex(2)});
            m.order.push_back(2);
        }
        else if (id == bcBase || id == bcDireita || id == bcSuperficie || id == bcEsquerda) // sem os gemeos
            m.bc[EdgeKey(gel->NodeIndex(0), gel->NodeIndex(1))] = id;
    }
    return m;
}

TPZGeoMesh *ToGeoMesh(const FlatMesh &m)
{
    auto *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(m.nodes.size());
    for (size_t i = 0; i < m.nodes.size(); i++) {
        TPZManVector<REAL, 3> x({m.nodes[i][0], m.nodes[i][1], 0.});
        gmesh->NodeVec()[i].Initialize(x, *gmesh);
    }
    int64_t index;
    for (const auto &t : m.tri) { // primeiro os triangulos: o elemento i e m.tri[i]
        TPZManVector<int64_t, 3> nos({t[0], t[1], t[2]});
        gmesh->CreateGeoElement(ETriangle, nos, matSolo, index);
    }
    for (const auto &[e, id] : m.bc) {
        TPZManVector<int64_t, 2> nos({e.first, e.second});
        gmesh->CreateGeoElement(EOned, nos, id, index);
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

// Refinamento uniforme: cada triangulo vira 4 semelhantes, pelos pontos medios das arestas (conforme)
FlatMesh RefineUniform(const FlatMesh &in)
{
    FlatMesh m;
    m.nodes = in.nodes;
    std::map<std::pair<int64_t, int64_t>, int64_t> mid; // aresta -> no do ponto medio
    auto midpoint = [&](int64_t a, int64_t b) {
        const auto e = EdgeKey(a, b);
        auto it = mid.find(e);
        if (it != mid.end()) return it->second;
        m.nodes.push_back({(m.nodes[a][0] + m.nodes[b][0]) / 2., (m.nodes[a][1] + m.nodes[b][1]) / 2.});
        return mid[e] = (int64_t)m.nodes.size() - 1;
    };
    for (size_t i = 0; i < in.tri.size(); i++) {
        const auto &t = in.tri[i];
        const int64_t a = midpoint(t[0], t[1]), b = midpoint(t[1], t[2]), c = midpoint(t[2], t[0]);
        m.tri.push_back({t[0], a, c});
        m.tri.push_back({a, t[1], b});
        m.tri.push_back({c, b, t[2]});
        m.tri.push_back({a, b, c});
        m.order.insert(m.order.end(), 4, in.order[i]); // os filhos herdam a ordem
    }
    for (const auto &[e, id] : in.bc) {
        const int64_t c = midpoint(e.first, e.second);
        m.bc[EdgeKey(e.first, c)] = id;
        m.bc[EdgeKey(c, e.second)] = id;
    }
    return m;
}

// Bisseccao pela aresta mais longa (Rivara): cada triangulo marcado e dividido duas vezes (h/2) e os vizinhos so o
// necessario para a malha continuar conforme
FlatMesh RefineMarked(const FlatMesh &in, const std::set<int64_t> &marked)
{
    using Edge = std::pair<int64_t, int64_t>;
    FlatMesh m = in;
    std::vector<int> nbis(m.tri.size(), 0); // bisseccoes que faltam em cada triangulo
    for (int64_t t : marked) nbis[t] = 2;
    std::vector<bool> alive(m.tri.size(), true);
    std::map<Edge, std::vector<int64_t>> edgeTri; // aresta -> triangulos vivos

    auto edge = [&](int64_t t, int i) { return EdgeKey(m.tri[t][i], m.tri[t][(i + 1) % 3]); };
    auto addTri = [&](const std::array<int64_t, 3> &v, int n, int o) {
        m.tri.push_back(v);
        m.order.push_back(o);
        nbis.push_back(n);
        alive.push_back(true);
        for (int i = 0; i < 3; i++) edgeTri[edge(m.tri.size() - 1, i)].push_back(m.tri.size() - 1);
    };
    for (int64_t t = 0; t < (int64_t)m.tri.size(); t++)
        for (int i = 0; i < 3; i++) edgeTri[edge(t, i)].push_back(t);

    // aresta mais longa; no empate, a de maior chave (os dois triangulos da aresta escolhem igual)
    auto longest = [&](int64_t t) {
        int best = 0;
        REAL l2best = -1.;
        for (int i = 0; i < 3; i++) {
            const Edge e = edge(t, i);
            const REAL dx = m.nodes[e.second][0] - m.nodes[e.first][0], dy = m.nodes[e.second][1] - m.nodes[e.first][1];
            const REAL l2 = dx * dx + dy * dy;
            if (l2 > l2best || (l2 == l2best && e > edge(t, best))) best = i, l2best = l2;
        }
        return best;
    };
    auto neighbor = [&](int64_t t, const Edge &e) -> int64_t {
        for (int64_t s : edgeTri[e]) if (s != t) return s;
        return -1;
    };
    // divide s pela aresta e no no mid: (a, mid, c) e (mid, b, c) mantem a orientacao
    auto split = [&](int64_t s, const Edge &e, int64_t mid) {
        const std::array<int64_t, 3> v = m.tri[s];
        int i = 0;
        while (edge(s, i) != e) i++;
        for (int j = 0; j < 3; j++) {
            auto &lista = edgeTri[edge(s, j)];
            lista.erase(std::find(lista.begin(), lista.end(), s));
        }
        alive[s] = false;
        const int n = std::max(nbis[s] - 1, 0), o = m.order[s]; // os filhos herdam a ordem
        addTri({v[i], mid, v[(i + 2) % 3]}, n, o);
        addTri({mid, v[(i + 1) % 3], v[(i + 2) % 3]}, n, o);
    };
    std::function<void(int64_t)> bisect = [&](int64_t t) {
        Edge e;
        int64_t n;
        for (;;) { // antes, o vizinho pela aresta mais longa tem que ter essa mesma aresta como a mais longa
            e = edge(t, longest(t));
            n = neighbor(t, e);
            if (n < 0 || edge(n, longest(n)) == e) break;
            bisect(n);
        }
        const int64_t mid = m.nodes.size();
        m.nodes.push_back({(m.nodes[e.first][0] + m.nodes[e.second][0]) / 2.,
                           (m.nodes[e.first][1] + m.nodes[e.second][1]) / 2.});
        split(t, e, mid);
        if (n >= 0) split(n, e, mid);
        auto it = m.bc.find(e);
        if (it != m.bc.end()) { // aresta de contorno: os dois pedacos ficam com o mesmo id
            const int id = it->second;
            m.bc.erase(it);
            m.bc[EdgeKey(e.first, mid)] = id;
            m.bc[EdgeKey(mid, e.second)] = id;
        }
    };
    for (bool again = true; again;) {
        again = false;
        for (int64_t t = 0; t < (int64_t)m.tri.size(); t++)
            if (alive[t] && nbis[t] > 0) bisect(t), again = true;
    }

    FlatMesh out;
    out.nodes = m.nodes;
    out.bc = m.bc;
    for (int64_t t = 0; t < (int64_t)m.tri.size(); t++)
        if (alive[t]) out.tri.push_back(m.tri[t]), out.order.push_back(m.order[t]);
    return out;
}

// sqrt(J2) da deformacao plastica de um ponto (cisalhamento de engenharia -> tensor)
REAL SqrtJ2Plastic(const TPZElastoPlasticMem &mem)
{
    TPZTensor<REAL> ep = mem.m_elastoplastic_state.m_eps_p;
    ep.XY() *= 0.5; ep.XZ() *= 0.5; ep.YZ() *= 0.5;
    return std::sqrt(std::max<REAL>(ep.J2(), 0.));
}

// Indicador de SlopeAnalysis: maior sqrt(J2(eps_p))  nos pontos de integracao de cada elemento
template <class TMat>
std::vector<REAL> PlasticIndicator(TPZMultiphysicsCompMesh *cmesh)
{
    auto *mat = dynamic_cast<TMat *>(cmesh->FindMaterial(matSolo));
    if (!mat) DebugStop();
    std::vector<REAL> ind(cmesh->Reference()->NElements(), 0.);
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        REAL &v = ind[cel->Reference()->Index()];
        for (int64_t im : mem)
            if (im >= 0) v = std::max(v, SqrtJ2Plastic(mat->MemItem(im)));
    }
    return ind;
}

// Ordem da regra de integracao de um elemento multifisico (a do material, com a ordem de cada espaco)
int RuleOrder(TPZCompEl *cel)
{
    auto *mf = dynamic_cast<TPZMultiphysicsElement *>(cel);
    auto *mat = dynamic_cast<TPZMatCombinedSpaces *>(cel->Material());
    if (!mf || !mat) DebugStop();
    TPZManVector<int, 3> ordens(mf->NMeshes(), 0);
    for (int64_t i = 0; i < mf->NMeshes(); i++)
        if (auto *sp = dynamic_cast<TPZInterpolationSpace *>(mf->Element(i))) ordens[i] = sp->MaxOrder();
    return mat->IntegrationRuleOrder(ordens);
}

// CSV com x, y e sqrt(J2(eps_p)) de cada ponto de integracao (a regra de cada elemento)
template <class TMat>
void WritePlasticPoints(TPZMultiphysicsCompMesh *cmesh, const std::string &file)
{
    auto *mat = dynamic_cast<TMat *>(cmesh->FindMaterial(matSolo));
    if (!mat) DebugStop();
    std::ofstream out(file);
    out << "x,y,sqrtJ2\n";
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->MaterialId() != matSolo) continue;
        TPZGeoEl *gel = cel->Reference();
        const int o = RuleOrder(cel);
        std::unique_ptr<TPZIntPoints> rule(gel->CreateSideIntegrationRule(gel->NSides() - 1, o));
        TPZManVector<int, 2> ordem(2, o);
        rule->SetOrder(ordem);
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        if (mem.size() != rule->NPoints()) DebugStop(); // a memoria segue a ordem dos pontos da regra
        for (int ip = 0; ip < rule->NPoints(); ip++) {
            TPZManVector<REAL, 3> qsi(2), x(3);
            REAL w;
            rule->Point(ip, qsi, w);
            gel->X(qsi, x);
            out << x[0] << "," << x[1] << "," << SqrtJ2Plastic(mat->MemItem(mem[ip])) << "\n";
        }
    }
}

// Maior f/(2c cos phi) de Mohr-Coulomb nas tensoes da memoria, f = (s1 - s3) + (s1 + s3) sin(phi) - 2c cos(phi):
// 0 = sobre a superficie, > 0 = fora dela
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

// segundos desde t0
using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// Silencia o std::cout no escopo: a montagem serial do NeoPZ imprime '*' a cada 1000 elementos
struct QuietCout {
    std::streambuf *old;
    QuietCout() : old(std::cout.rdbuf(nullptr)) {}
    ~QuietCout() { std::cout.rdbuf(old); }
};

// Linha formatada como printf
template <class... A>
void Print(const char *fmt, A... a)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, fmt, a...);
    std::cout << buf << std::endl;
}

// O que o GIM usa de cada material
struct GIMProblem {
    TPZLinearAnalysis *an = nullptr;
    std::vector<int64_t> free, fixed;       // equacoes livres e fixas (apoios e pressoes)
    TPZFMatrix<STATE> x0;                   // valores das fixas com lambda = 1 (pressoes da percolacao; 0 no seco)
    bool scaleFixed = false;                // fixas = lambda x0 (u-p); senao o material escala as pressoes
    REAL lambda = 1.;
    std::function<void(REAL)> setLambda;    // peso = lambda gamma
    std::function<void(REAL)> setReduction; // c/F e tan(phi)/F (F != 1 so no preditor)
    std::function<void()> accept;           // memoria dos pontos <- solucao carregada na malha
    void SetLambda(REAL lam)
    {
        lambda = lam;
        setLambda(lam);
    }
    REAL Dot(const TPZFMatrix<STATE> &a, const TPZFMatrix<STATE> &b) const
    {
        REAL s = 0.;
        for (int64_t i : free) s += a.GetVal(i, 0) * b.GetVal(i, 0);
        return s;
    }
};

// Filtro da matriz: so as equacoes livres entram no sistema
void SetFreeEquations(TPZLinearAnalysis &an, const std::vector<int64_t> &free)
{
    TPZStack<int64_t> act;
    for (int64_t i : free) act.Push(i);
    an.StructMatrix()->EquationFilter().Reset();
    an.StructMatrix()->EquationFilter().SetActiveEquations(act);
}

TPZFMatrix<STATE> SolveK(TPZLinearAnalysis &an, const TPZFMatrix<STATE> &b);

// Resolve um problema linear nas equacoes livres; as fixas ficam com o valor de x (a percolacao)
bool SolveFree(TPZLinearAnalysis &an, const std::vector<int64_t> &free, TPZFMatrix<STATE> &x)
{
    SetFreeEquations(an, free);
    REAL r0 = 0.;
    for (int it = 0; it < 5; it++) {
        {
            QuietCout q;
            an.LoadSolution(x);
            an.Assemble();
        }
        const TPZFMatrix<STATE> rhs = an.Rhs();
        REAL r = 0.;
        for (int64_t i : free) r += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
        r = std::sqrt(r);
        if (it == 0) r0 = r;
        if (r <= 1.e-10 * r0) return true;
        x += SolveK(an, rhs);
    }
    return false;
}

// -R(x) nas equacoes livres (a matriz tangente fica no solver)
TPZFMatrix<STATE> AssembleAt(GIMProblem &p, const TPZFMatrix<STATE> &x)
{
    TPZFMatrix<STATE> xl(x);
    if (p.scaleFixed)
        for (int64_t i : p.fixed) xl(i, 0) = p.lambda * p.x0.GetVal(i, 0);
    QuietCout q;
    p.an->LoadSolution(xl);
    p.an->Assemble();
    return p.an->Rhs();
}

// K y = b com a fatoracao corrente (Solve do TPZLinearAnalysis carregaria y na malha)
TPZFMatrix<STATE> SolveK(TPZLinearAnalysis &an, const TPZFMatrix<STATE> &b)
{
    TPZEquationFilter &filtro = an.StructMatrix()->EquationFilter();
    TPZFMatrix<STATE> br(filtro.NActiveEquations(), 1, 0.), yr, y(b.Rows(), 1, 0.);
    filtro.Gather(b, br);
    an.MatrixSolver<STATE>().Solve(br, yr);
    filtro.Scatter(yr, y);
    return y;
}

const REAL kTolGIM = 1.e-7; // |R| / |fbar| nas equacoes livres (a tolerancia do TPZUnsatAnalysis)


// Newton com lambda prescrito, a partir de x
bool NewtonGIM(GIMProblem &p, TPZFMatrix<STATE> &x, REAL fref, REAL lambda, int maxit, int &it)
{
    p.SetLambda(lambda);
    for (int k = 1; k <= maxit; k++) {
        it = k;
        const TPZFMatrix<STATE> rhs = AssembleAt(p, x);
        const REAL e = std::sqrt(p.Dot(rhs, rhs)) / fref;
        if (!std::isfinite(e)) return false;
        if (k > 1 && e < kTolGIM) return true;
        try { x += SolveK(*p.an, rhs); } catch (...) { return false; } // tangente singular
    }
    return false;
}

struct ArcResult {
    bool conv = false, vertical = false;
    int iter = 0;
    REAL lambda = 0.;
    TPZFMatrix<STATE> du, x; // incremento e estado final
};

// Um incremento do arc-length cilindrico (de Souza Neto et al. 2008, Box 4.4), como TPZUnsatAnalysis::ArcLengthIterate:
// preditor secante (o incremento anterior com comprimento l) e raiz (4.118) de menor angulo
ArcResult ArcLengthStep(GIMProblem &p, const TPZFMatrix<STATE> &X, const TPZFMatrix<STATE> &fbar, REAL fref,
                        REAL lambdan, REAL l, const TPZFMatrix<STATE> &dun, REAL dlambdan, int maxit)
{
    ArcResult res;
    res.lambda = lambdan;
    res.du.Redim(X.Rows(), 1);
    TPZFMatrix<STATE> x(X), &Du = res.du;
    REAL lambda = lambdan;
    const bool hasPrev = dun.Rows() == X.Rows();
    int kfirst = 1;
    if (hasPrev) {
        const REAL prev2 = p.Dot(dun, dun);
        if (prev2 > 0. && l * l < 1.e4 * prev2) {
            const REAL s = l / std::sqrt(prev2);
            for (int64_t i : p.free) x(i, 0) += s * dun.GetVal(i, 0), Du(i, 0) = s * dun.GetVal(i, 0);
            lambda += s * dlambdan;
            kfirst = 2;
        }
    }
    for (int k = kfirst; k <= maxit; k++) {
        res.iter = k - kfirst + 1;
        if (!(lambda > 0.)) break;
        p.SetLambda(lambda);
        p.setReduction(k == 1 ? 1. + 1.e-6 : 1.); // preditor: superficie um pouco menor, tangente plastica nela
        const TPZFMatrix<STATE> rhs = AssembleAt(p, x);
        const REAL e = std::sqrt(p.Dot(rhs, rhs)) / fref;
        if (!std::isfinite(e) || e > 1.e6) break;
        if (k > 1 && e < kTolGIM) {
            res.conv = true;
            break;
        }
        TPZFMatrix<STATE> dustar, dubar;
        try {
            dustar = SolveK(*p.an, rhs);
            dubar = SolveK(*p.an, fbar);
        } catch (...) {
            break;
        }
        const REAL a = p.Dot(dubar, dubar);
        if (!(a > 0.)) {
            res.vertical = a == 0.;
            break;
        }
        REAL dlambda;
        if (k == 1) { // preditor tangente, com o sinal do criterio secante (4.123)
            dlambda = (hasPrev && p.Dot(dun, dubar) < 0. ? -1. : 1.) * l / std::sqrt(a);
        } else {      // restricao |Du| = l (4.116)-(4.118)
            TPZFMatrix<STATE> w(Du);
            w += dustar;
            const REAL b = 2. * p.Dot(w, dubar), c = p.Dot(w, w) - l * l, disc = b * b - 4. * a * c;
            if (disc < 0.) break;
            const REAL q = -0.5 * (b + (b >= 0. ? 1. : -1.) * std::sqrt(disc));
            const REAL r1 = q / a, r2 = q != 0. ? c / q : r1, dir = p.Dot(dubar, Du);
            dlambda = r1 * dir >= r2 * dir ? r1 : r2;
        }
        lambda += dlambda;
        for (int64_t i : p.free) {
            const REAL di = dustar.GetVal(i, 0) + dlambda * dubar.GetVal(i, 0);
            x(i, 0) += di;
            Du(i, 0) += di;
        }
    }
    p.setReduction(1.);
    res.lambda = lambda;
    res.x = x;
    return res;
}

// GIM com arc-length, parametros do FactorGIMArco do SlopeRandomField: lambda = 0,5 prescrito, depois o caminho com l
// ajustado para ~10 iteracoes, ate d lambda < 0,001 (no l pedido), lambda >= 10, l < 1e-4 |u| ou 100 incrementos.
// Gamma = maior lambda do caminho
struct GIMOut {
    REAL gamma = 0.;
    int steps = 0;
    std::string end;
};

GIMOut GIMArcLength(GIMProblem &p, int64_t neq)
{
    const REAL tol = 0.001, l0 = 0.25, lMaxR = 2., lMinR = 1.e-4, lamLimite = 10.;
    const int alvo = 10, maxIt = 12, maxPassos = 100;
    if (p.x0.Rows() != neq) p.x0.Redim(neq, 1); // seco: fixas em 0
    std::vector<bool> livre(neq, false);
    for (int64_t i : p.free) livre[i] = true;
    p.fixed.clear();
    for (int64_t i = 0; i < neq; i++)
        if (!livre[i]) p.fixed.push_back(i);
    TPZFMatrix<STATE> X(p.x0); // estado aceito: sem tensoes, com as pressoes da percolacao
    p.setReduction(1.);
    p.SetLambda(1.);
    TPZFMatrix<STATE> fbar = AssembleAt(p, X); // carga por unidade de lambda: -R(lambda = 1) + R(lambda = 0)
    p.SetLambda(0.);
    fbar -= AssembleAt(p, X);
    const REAL fref = std::sqrt(p.Dot(fbar, fbar));
    const auto t0 = Clock::now();
    Print("    %-6s %10s %10s %4s %11s %9s", "passo", "lambda", "d lambda", "it", "arco l", "tempo (s)");

    REAL lam = 0.5;
    for (;;) {
        TPZFMatrix<STATE> x(X);
        int it = 0;
        const bool ok = NewtonGIM(p, x, fref, lam, maxIt, it);
        Print("    %-6s %10.5f %10s %4d %11s %9.1f%s", "fixo", lam, "", it, "", Seconds(t0), ok ? "" : "  falhou");
        if (ok) {
            p.accept();
            X = x;
            break;
        }
        lam *= 0.5;
        if (lam < 0.01) return {0., 0, "sem equilibrio inicial"};
    }
    const REAL uRef = std::sqrt(p.Dot(X, X)), lMax = lMaxR * uRef, lMin = lMinR * uRef;
    REAL l = std::min(lMax, std::max(lMin, l0 * uRef)), lPedido = l, lamMax = lam, dlamn = 0.;
    TPZFMatrix<STATE> dun;
    std::string fim = "maxPassos";
    int n = 1;
    for (; n <= maxPassos; n++) {
        const ArcResult r = ArcLengthStep(p, X, fbar, fref, lam, l, dun, dlamn, maxIt);
        char dl[16] = "";
        if (r.conv) std::snprintf(dl, sizeof dl, "%+10.5f", r.lambda - lam);
        Print("    %-6d %10.5f %10s %4d %11.3e %9.1f%s", n, r.lambda, dl, r.iter, l, Seconds(t0), r.conv ? "" : "  falhou");
        if (!r.conv) {
            if (r.vertical) { fim = "vertical"; break; }
            l *= 0.5;
            if (l < lMin) { fim = "lMin"; break; }
            continue;
        }
        p.accept();
        X = r.x;
        dlamn = r.lambda - lam;
        lam = r.lambda;
        dun = r.du;
        lamMax = std::max(lamMax, lam);
        if (dlamn * lPedido / l < tol) { fim = "patamar"; break; }
        if (lam >= lamLimite) { fim = "lambda >= 10"; break; }
        l = std::min(lMax, l * std::min(2., std::max(0.5, std::sqrt(REAL(alvo) / r.iter))));
        lPedido = l;
    }
    p.SetLambda(1.);
    return {lamMax, std::min(n, maxPassos), fim};
}

// Resultado de um nivel: Gamma, equacoes, indicador plastico de cada elemento e p_w de cada no (rebaixamento)
struct LevelResult {
    REAL gamma = 0., fmc = 0.; // fmc: maior f_MC / (2c cos phi) no fim (0: tensoes sobre a superficie)
    int64_t neq = 0, nfree = 0;
    int steps = 0;
    std::string end;
    std::vector<REAL> ind, pw;
};

// Poropressao imposta em cada no da superficie e das laterais (a base fica sem fluxo)
std::map<int64_t, REAL> BoundaryPressures(TPZGeoMesh *gmesh, const Data &d)
{
    std::map<int64_t, REAL> pb;
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel) continue;
        const int id = gel->MaterialId();
        if (id != bcSuperficie && id != bcEsquerda && id != bcDireita) continue;
        for (int i = 0; i < gel->NCornerNodes(); i++)
            pb[gel->NodeIndex(i)] = ContourPressure(id, gel->NodePtr(i)->Coord(1), d);
    }
    return pb;
}

// Cabecalho do nivel: triangulos e ordens (cabecalho) com as equacoes da malha computacional
void PrintLevelHeader(const std::string &cabecalho, int64_t neq, int64_t nfree)
{
    std::cout << cabecalho << " | " << neq << " equacoes (" << nfree << " livres)" << std::endl;
}

void PrintSeepage(const std::vector<REAL> &pw)
{
    Print("    percolacao (esqueleto rigido): p_w de %.1f a %.1f kPa", *std::min_element(pw.begin(), pw.end()),
          *std::max_element(pw.begin(), pw.end()));
}

// Material saturado u-p. Rebaixamento: percolacao com o esqueleto rigido e GIM com essas pressoes (lambda escala o
// peso e a pressao); seco: pressoes 0. Grava os pontos de integracao em pointsFile
LevelResult SolveWithSatMat(TPZGeoMesh *gmesh, const Data &d, const std::string &pointsFile, bool rebaix,
                            const std::vector<int> &ordens, const std::string &cabecalho)
{
    TPZMultiphysicsCompMesh *cmesh = CreateUPCompMesh(gmesh, 2, d, ordens);
    auto *mat = dynamic_cast<TMatUP *>(cmesh->FindMaterial(matSolo));
    if (!mat) DebugStop();
    LevelResult res;
    {
        TPZPoroElastoPlasticUPAnalysis an(cmesh, mat);
        TPZSpStructMatrix<STATE> strmat(cmesh); // PARDISO, nao simetrica; reordena as equacoes sozinho
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        an.SetSolver(solver);

        an.IdentifyEquations(); // equacoes de pressao e com Dirichlet
        const int64_t neq = cmesh->NEquations();
        std::vector<int64_t> livres; // do GIM: deslocamentos sem Dirichlet
        for (int64_t eq = 0; eq < neq; eq++)
            if (!an.PressureEquations()[eq] && !an.ConstrainedEquations()[eq]) livres.push_back(eq);
        PrintLevelHeader(cabecalho, neq, livres.size());
        TPZFMatrix<STATE> x0(neq, 1, 0.);
        if (rebaix) { // percolacao: u = 0 e p do contorno fixos; livres as pressoes do interior e da base
            std::set<int64_t> contorno;
            for (const auto &[no, v] : BoundaryPressures(gmesh, d)) {
                const int64_t eq = an.NodeEquation(no, 1, 0);
                x0(eq, 0) = v;
                contorno.insert(eq);
            }
            std::vector<int64_t> livresP;
            for (int64_t eq = 0; eq < neq; eq++)
                if (an.PressureEquations()[eq] && !contorno.count(eq)) livresP.push_back(eq);
            mat->SetTimeStep(1.); // liga o fluxo; sem armazenamento (1/M = 0) o passo ja e o permanente
            if (!SolveFree(an, livresP, x0)) DebugStop();
            mat->SetTimeStep(0.);
            for (int64_t no = 0; no < gmesh->NNodes(); no++) res.pw.push_back(x0(an.NodeEquation(no, 1, 0), 0));
            PrintSeepage(res.pw);
        }
        GIMProblem p;
        p.an = &an;
        p.x0 = x0;
        p.scaleFixed = true; // pressoes fixas = lambda p
        p.free = livres;
        p.setLambda = [mat, &d](REAL lam) { mat->SetBodyForce(TPZManVector<REAL, 3>({0., -lam * d.gamma, 0.})); };
        p.setReduction = [mat](REAL F) { mat->PlasticModel().SetStrengthReductionFactor(F); };
        p.accept = [&an]() { QuietCout q; an.AcceptSolution(); };
        SetFreeEquations(an, p.free);
        res.neq = cmesh->NEquations();
        res.nfree = p.free.size();
        const GIMOut g = GIMArcLength(p, res.neq);
        res.gamma = g.gamma, res.steps = g.steps, res.end = g.end;
    }
    res.fmc = MaxYield(mat, d);
    res.ind = PlasticIndicator<TMatUP>(cmesh);
    WritePlasticPoints<TMatUP>(cmesh, pointsFile);
    EraseMesh(cmesh);
    return res;
}

// Material nao saturado: como SolveWithSatMat; o material escala as pressoes congeladas com lambda
LevelResult SolveWithUnsatMat(TPZGeoMesh *gmesh, const Data &d, const std::string &pointsFile, bool rebaix,
                              const std::vector<int> &ordens, const std::string &cabecalho)
{
    TPZMultiphysicsCompMesh *cmesh = CreateUnsatCompMesh(gmesh, 2, d, ordens);
    LevelResult res;
    {
        // analise antes da matriz: o construtor renumera as equacoes
        TPZUnsatAnalysis an(cmesh, {{bcBase, {true, true}}, {bcDireita, {true, false}}, {bcEsquerda, {true, false}}});
        TPZSpStructMatrix<STATE> strmat(cmesh); // PARDISO, nao simetrica; reordena as equacoes sozinho
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> solver;
        solver.SetDirect(ELU);
        an.SetSolver(solver);

        std::vector<int64_t> livres; // do GIM: deslocamentos sem apoio (pressoes congeladas)
        const auto fixed = an.FrozenFixed(an.X());
        for (int64_t eq = 0; eq < an.NEq(); eq++)
            if (!fixed.count(eq)) livres.push_back(eq);
        PrintLevelHeader(cabecalho, an.NEq(), livres.size());
        TPZFMatrix<STATE> x0(an.NEq(), 1, 0.);
        if (rebaix) { // percolacao: u = 0, p_a = 0 e p_w do contorno fixos
            auto fixas = an.AllDisplacements(an.X());
            for (const auto &[no, eq] : an.AirEq()) fixas[eq] = 0.;
            for (const auto &[no, v] : BoundaryPressures(gmesh, d)) fixas[an.WaterEq().at(no)] = v;
            {
                QuietCout q;
                if (!an.Solve(1., fixas, {}, false).conv) DebugStop(); // sem armazenamento: o permanente
                an.Accept();
            }
            x0 = an.X();
            an.ResetMemory(); // o GIM parte do estado sem tensoes
            for (int64_t no = 0; no < gmesh->NNodes(); no++) res.pw.push_back(an.Nodal(x0, no, 1));
            PrintSeepage(res.pw);
        }
        an.SetGravityScalesSuction(true); // lambda escala o peso e toda a pressao (o Gamma do artigo)
        for (auto *m : an.Materials()) {
            m->SetTimeStep(0.);
            m->SetFrozenPressures(true); // so o equilibrio
        }
        GIMProblem p;
        p.an = &an;
        p.x0 = x0;
        p.scaleFixed = false;
        p.free = livres;
        p.setLambda = [&an](REAL lam) { an.SetGravityFactor(lam); };
        p.setReduction = [&an](REAL F) { an.SetStrengthReductionFactor(F); };
        p.accept = [&an]() { QuietCout q; an.Accept(); };
        SetFreeEquations(an, p.free);
        res.neq = an.NEq();
        res.nfree = p.free.size();
        const GIMOut g = GIMArcLength(p, res.neq);
        res.gamma = g.gamma, res.steps = g.steps, res.end = g.end;
    }
    res.fmc = MaxYield(dynamic_cast<TUnsat *>(cmesh->FindMaterial(matSolo)), d);
    res.ind = PlasticIndicator<TUnsat>(cmesh);
    WritePlasticPoints<TUnsat>(cmesh, pointsFile);
    EraseMesh(cmesh);
    return res;
}

// VTK de um nivel: triangulos, o indicador plastico de cada um e, no rebaixamento, p_w e o excesso u nos nos
void WriteLevelVTK(const FlatMesh &m, const std::vector<REAL> &ind, const std::vector<REAL> &pw, const Data &d,
                   const std::string &file)
{
    std::ofstream out(file);
    out << "# vtk DataFile Version 3.0\nGIM\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << m.nodes.size() << " double\n";
    for (const auto &x : m.nodes) out << x[0] << " " << x[1] << " 0\n";
    out << "CELLS " << m.tri.size() << " " << 4 * m.tri.size() << "\n";
    for (const auto &t : m.tri) out << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    out << "CELL_TYPES " << m.tri.size() << "\n";
    for (size_t i = 0; i < m.tri.size(); i++) out << "5\n";
    out << "CELL_DATA " << m.tri.size() << "\nSCALARS sqrtJ2_eps_p double 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < m.tri.size(); i++) out << ind[i] << "\n";
    out << "SCALARS ordem_u int 1\nLOOKUP_TABLE default\n";
    for (int o : m.order) out << o << "\n";
    if (pw.empty()) return;
    out << "POINT_DATA " << m.nodes.size() << "\nSCALARS p_w double 1\nLOOKUP_TABLE default\n";
    for (REAL v : pw) out << v << "\n";
    out << "SCALARS u_excesso double 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < m.nodes.size(); i++) out << ExcessPorePressure(pw[i], m.nodes[i][1], d) << "\n";
}

// Como FactorOfSafety de main(2).cpp: Gamma, marca sqrt(J2(eps_p)) >= 0,1 do maximo, refina; nref vezes. Os marcados
// sao divididos (h) e a ordem de u sobe 1 neles, ate pmax (p).
// Arquivos de cada nivel k: <prefixo>_<mat>_nivel<k>.vtk (malha), _pontos.csv (pontos de integracao); resumo em
// <prefixo>_<mat>.csv, com o alvo (para figuras_gim.py)
// Linha da tabela final
struct LevelSummary {
    int level = 0;
    int64_t ntri = 0, neq = 0, nfree = 0;
    std::map<int, int> orders; // ordem de u -> triangulos
    REAL gamma = 0.;
    int steps = 0;
    double seconds = 0.;
};

std::vector<LevelSummary> AdaptiveGIM(FlatMesh m, int nref, const Data &d, bool sat, const std::string &prefixo,
                                      REAL alvo, bool rebaix, int pmax)
{
    const std::string base = prefixo + (sat ? "_up" : "_unsat");
    std::ofstream resumo(base + ".csv");
    resumo << "nivel,triangulos,equacoes,livres,gamma,tempo_s,alvo\n";
    const auto t0 = Clock::now();
    std::vector<LevelSummary> niveis;
    for (int k = 0;; k++) {
        const auto tk = Clock::now();
        LevelSummary ls;
        ls.level = k;
        ls.ntri = m.tri.size();
        for (int o : m.order) ls.orders[o]++;
        std::string ordens;
        for (const auto &[o, n] : ls.orders) ordens += " P" + std::to_string(o) + " " + std::to_string(n);
        const std::string cabecalho =
            "\n  nivel " + std::to_string(k) + " | " + std::to_string(m.tri.size()) + " triangulos (u:" + ordens + ")";
        const std::string nivel = base + "_nivel" + std::to_string(k);
        TPZGeoMesh *gmesh = ToGeoMesh(m);
        DuplicateBoundaryBC(gmesh);
        const LevelResult r = sat ? SolveWithSatMat(gmesh, d, nivel + "_pontos.csv", rebaix, m.order, cabecalho)
                                  : SolveWithUnsatMat(gmesh, d, nivel + "_pontos.csv", rebaix, m.order, cabecalho);
        delete gmesh;
        ls.neq = r.neq, ls.nfree = r.nfree, ls.gamma = r.gamma, ls.steps = r.steps, ls.seconds = Seconds(tk);
        niveis.push_back(ls);
        WriteLevelVTK(m, r.ind, r.pw, d, nivel + ".vtk");
        resumo << k << "," << m.tri.size() << "," << r.neq << "," << r.nfree << "," << r.gamma << "," << Seconds(tk)
               << "," << alvo << std::endl;
        Print("    => Gamma = %.5f (erro %+.2f%%) | %d passos, fim: %s | f_MC max %.0e | %.1f s (total %.1f s)",
              r.gamma, 100. * (r.gamma / alvo - 1.), r.steps, r.end.c_str(), r.fmc, Seconds(tk), Seconds(t0));
        if (k >= nref) break;
        const REAL vmax = *std::max_element(r.ind.begin(), r.ind.begin() + m.tri.size());
        std::set<int64_t> marked;
        for (int64_t t = 0; t < (int64_t)m.tri.size(); t++)
            if (vmax > 0. && r.ind[t] >= 0.1 * vmax) marked.insert(t);
        for (int64_t t : marked) m.order[t] = std::min(m.order[t] + 1, pmax); // p: os filhos herdam
        m = RefineMarked(m, marked);                                         // h
    }
    return niveis;
}

// Tabela final no terminal: os niveis dos dois materiais
void PrintSummary(const std::string &titulo, REAL alvo, const std::vector<LevelSummary> &up,
                  const std::vector<LevelSummary> &unsat, double total)
{
    int pmaxv = 2;
    for (const auto *v : {&up, &unsat})
        for (const auto &l : *v) pmaxv = std::max(pmaxv, l.orders.rbegin()->first);
    std::string traco(88 + 7 * (pmaxv - 1), '-'), duplo(traco.size(), '=');
    std::cout << "\n" << duplo << "\n RESUMO | " << titulo << " | alvo Gamma = " << alvo << "\n" << traco << "\n";
    std::cout << " material  nivel  triang.  equacoes   livres";
    for (int o = 2; o <= pmaxv; o++) std::cout << "     P" << o;
    std::cout << "      Gamma      erro  passos  tempo (s)\n";
    for (const auto &[nome, v] : {std::pair{"u-p", &up}, std::pair{"unsat", &unsat}})
        for (const auto &l : *v) {
            char buf[64];
            std::snprintf(buf, sizeof buf, " %-8s %6d %8lld %9lld %8lld", nome, l.level, (long long)l.ntri,
                          (long long)l.neq, (long long)l.nfree);
            std::cout << buf;
            for (int o = 2; o <= pmaxv; o++) {
                auto it = l.orders.find(o);
                std::snprintf(buf, sizeof buf, " %6d", it == l.orders.end() ? 0 : it->second);
                std::cout << buf;
            }
            Print(" %10.5f %+8.2f%% %7d %10.1f", l.gamma, 100. * (l.gamma / alvo - 1.), l.steps, l.seconds);
        }
    const REAL gu = up.back().gamma, gn = unsat.back().gamma;
    std::cout << traco << "\n";
    Print(" Gamma final: u-p %.5f (%+.2f%%) | nao saturado %.5f (%+.2f%%) | diferenca %.1e", gu,
          100. * (gu / alvo - 1.), gn, 100. * (gn / alvo - 1.), std::fabs(gu - gn));
    Print(" tempo total: %.1f s", total);
    std::cout << duplo << std::endl;
}

// Casos (argumento do programa): malha, geometria e alvo do GIM seco. Alvo: o 1,777 da Secao 5.3.2, que e do talude de
// H = 10 m da Fig. 19 (o texto diz 5 m); como Gamma = N c / (gamma H), com H = 5 m o alvo e 2 x 1,777 = 3,554
// Rebaixamento total (h_w = H): Gamma = 1,336 (Tabelas 5 e 6); no talude10 so 1,336 / 2 (Gamma ~ 1/H), indicativo,
// porque o dominio e outro e a percolacao muda
struct Case {
    std::string mesh;
    REAL H, yCrista, targetDry, targetDrawdown;
};
const std::map<std::string, Case> kCases = {
    {"ceron", {std::string(MESH_DIR) + "/ceron_h5_a45_g15.msh", 5., 10., 3.554, 1.336}}, // talude do relatorio
    {"talude10", {std::string(PROJECT_DIR) + "/malhas/talude_h10_grossa.msh", 10., 40., 1.777, 0.668}}, // Fig. 8.30
    {"ceron10H", {std::string(PROJECT_DIR) + "/malhas/ceron_h5_a45_10H.msh", 5., 55., 3.554, 1.336}}}; // dominio 10H

int main(int argc, char *argv[])
{
    const std::string caso = argc > 1 ? argv[1] : "ceron10H";
    const std::string modo = argc > 2 ? argv[2] : "rebaix";
    const int nunif = argc > 3 ? std::atoi(argv[3]) : 0; // refinamentos uniformes da malha inicial
    const int nref = argc > 4 ? std::atoi(argv[4]) : 4;  // niveis de refinamento adaptativo
    const int pmax = argc > 5 ? std::atoi(argv[5]) : 3;  // ordem maxima de u no refinamento p
    if (!kCases.count(caso) || (modo != "seco" && modo != "rebaix") || nunif < 0 || nref < 0 || pmax < 2) {
        std::cerr << "Uso: SlopeRandomFieldRefinement [ceron | talude10 | ceron10H] [seco | rebaix] [nunif = 1] [nref = 3] [pmax = 3]\n";
        return 1;
    }
    const Case &cs = kCases.at(caso);
    const bool rebaix = modo == "rebaix";
    const REAL alvo = rebaix ? cs.targetDrawdown : cs.targetDry;

    // ---- passo 1: malha geometrica
    const std::string arquivo = cs.mesh;
    const std::string nome = arquivo.substr(arquivo.rfind('/') + 1);
    std::cout << "Caso " << caso << ": lendo " << arquivo << "\n";

    TPZGeoMesh *gmesh = ReadGMesh(arquivo);
    GMeshSummary(gmesh);

    // VTK: cada elemento leva o seu id de material (matColor = true), para colorir no ParaView
    const std::string vtk = "malha_" + nome.substr(0, nome.rfind('.')) + ".vtk";
    std::ofstream out(vtk);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    std::cout << "Gravado " << vtk << " (abra no ParaView)\n";

    DuplicateBoundaryBC(gmesh); // gemeos -12 e -14, antes das malhas computacionais

    // ---- passo 2: malhas computacionais com os dois materiais
    Data d;
    d.H = cs.H;
    d.yCrista = cs.yCrista;
    d.hw = cs.H; // rebaixamento total
    const int porder = 2; // u quadratico, pressoes lineares

    std::cout << "\nMalha u-p (TPZMatPoroElastoPlasticUP):\n";
    TPZMultiphysicsCompMesh *malhaUP = CreateUPCompMesh(gmesh, porder, d);
    CMeshSummary<TMatUP>(malhaUP, {"u", "p_w"});

    std::cout << "\nMalha u-p_w-p_a (TPZMatUnsaturatedUPwPa):\n";
    TPZMultiphysicsCompMesh *malhaUnsat = CreateUnsatCompMesh(gmesh, porder, d);
    CMeshSummary<TUnsat>(malhaUnsat, {"u", "p_w", "p_a"});

    EraseMesh(malhaUP);
    EraseMesh(malhaUnsat);

    // ---- passo 3 (seco) e 4 (rebaixamento): GIM com arc-length e refinamento adaptativo
    FlatMesh m0 = FromGeoMesh(gmesh);
    for (int i = 0; i < nunif; i++) m0 = RefineUniform(m0);
    std::cout << "\nMalha inicial: " << m0.tri.size() << " triangulos (" << nunif << " refinamento(s) uniforme(s)), "
              << nref << " niveis adaptativos" << std::endl;
    const std::string prefixo = "gim_" + modo + "_" + caso + "_u" + std::to_string(nunif) + "h" + std::to_string(nref) +
                                "p" + std::to_string(pmax);
    const std::string titulo = "caso " + caso + " | " + (rebaix ? "rebaixamento (h_w = H)" : "seco") + " | " +
                               std::to_string(nunif) + " unif. + " + std::to_string(nref) + " adapt., u ate P" +
                               std::to_string(pmax);
    const std::string linha(100, '=');
    const auto t0 = Clock::now();
    std::cout << "\n" << linha << "\n GIM | " << titulo << " | alvo Gamma = " << alvo << "\n" << linha << std::endl;
    std::cout << "\n--- material u-p (TPZMatPoroElastoPlasticUP) ---" << std::endl;
    const auto nivUP = AdaptiveGIM(m0, nref, d, true, prefixo, alvo, rebaix, pmax);
    std::cout << "\n--- material nao saturado (TPZMatUnsaturatedUPwPa) ---" << std::endl;
    const auto nivUnsat = AdaptiveGIM(m0, nref, d, false, prefixo, alvo, rebaix, pmax);
    PrintSummary(titulo, alvo, nivUP, nivUnsat, Seconds(t0));

    delete gmesh;
    return 0;
}
