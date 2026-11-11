#include "SlopeMesh.h"

#include "TPZGmshReader.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>

namespace slope {

TPZGeoMesh *ReadGMesh(const std::string &file)
{
    TPZManVector<std::map<std::string, int>, 4> names(4); // nome do grupo fisico -> id
    names[2]["Solo"] = matSolo;
    names[1]["Base"] = bcBase;
    names[1]["Direita"] = bcDireita;
    names[1]["Superficie"] = bcSuperficie;
    names[1]["Esquerda"] = bcEsquerda;
    TPZGmshReader reader;
    reader.SetDimNamePhysical(names);
    TPZGeoMesh *gmesh = reader.GeometricGmshMesh(file);
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        const int id = gel ? gel->MaterialId() : matSolo;
        if (id != matSolo && id != bcBase && id != bcDireita && id != bcSuperficie && id != bcEsquerda) {
            std::cerr << "Id " << id << " sem nome correspondente em ReadGMesh\n";
            DebugStop();
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

void DuplicateBoundaryBC(TPZGeoMesh *gmesh)
{
    const int64_t nel = gmesh->NElements(); // os gemeos entram no fim e nao sao revisitados
    for (int64_t i = 0; i < nel; i++) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == bcDireita || id == bcEsquerda) TPZGeoElBC(gel, gel->NSides() - 1, id - 10);
    }
}

static std::pair<int64_t, int64_t> EdgeKey(int64_t a, int64_t b) { return {std::min(a, b), std::max(a, b)}; }

FlatMesh FromGeoMesh(TPZGeoMesh *gmesh, int order)
{
    FlatMesh m;
    for (int64_t i = 0; i < gmesh->NNodes(); i++)
        m.nodes.push_back({gmesh->NodeVec()[i].Coord(0), gmesh->NodeVec()[i].Coord(1)});
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == matSolo) {
            m.tri.push_back({gel->NodeIndex(0), gel->NodeIndex(1), gel->NodeIndex(2)});
            m.order.push_back(order);
        } else if (id == bcBase || id == bcDireita || id == bcSuperficie || id == bcEsquerda) { // sem os gemeos
            m.bc[EdgeKey(gel->NodeIndex(0), gel->NodeIndex(1))] = id;
        }
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
    for (const auto &t : m.tri) {
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
        m.order.insert(m.order.end(), 4, in.order[i]);
    }
    for (const auto &[e, id] : in.bc) {
        const int64_t c = midpoint(e.first, e.second);
        m.bc[EdgeKey(e.first, c)] = id;
        m.bc[EdgeKey(c, e.second)] = id;
    }
    return m;
}

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
        for (int64_t s : edgeTri[e])
            if (s != t) return s;
        return -1;
    };
    // divide s pela aresta e no no mid: (a, mid, c) e (mid, b, c) mantem a orientacao; os filhos herdam a ordem
    auto split = [&](int64_t s, const Edge &e, int64_t mid) {
        const std::array<int64_t, 3> v = m.tri[s];
        int i = 0;
        while (edge(s, i) != e) i++;
        for (int j = 0; j < 3; j++) {
            auto &lista = edgeTri[edge(s, j)];
            lista.erase(std::find(lista.begin(), lista.end(), s));
        }
        alive[s] = false;
        const int n = std::max(nbis[s] - 1, 0), o = m.order[s];
        addTri({v[i], mid, v[(i + 2) % 3]}, n, o);
        addTri({mid, v[(i + 1) % 3], v[(i + 2) % 3]}, n, o);
    };
    std::function<void(int64_t)> bisect = [&](int64_t t) {
        Edge e;
        int64_t n;
        for (;;) { // o vizinho pela aresta mais longa tem que ter essa mesma aresta como a mais longa
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
        if (it != m.bc.end()) { // aresta de contorno: os dois pedacos com o mesmo id
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

void WriteFlatVTK(const FlatMesh &m, const std::map<std::string, std::vector<REAL>> &cellData,
                  const std::map<std::string, std::vector<REAL>> &pointData, const std::string &file)
{
    std::ofstream out(file);
    out << "# vtk DataFile Version 3.0\ntalude\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << m.nodes.size() << " double\n";
    for (const auto &x : m.nodes) out << x[0] << " " << x[1] << " 0\n";
    out << "CELLS " << m.tri.size() << " " << 4 * m.tri.size() << "\n";
    for (const auto &t : m.tri) out << "3 " << t[0] << " " << t[1] << " " << t[2] << "\n";
    out << "CELL_TYPES " << m.tri.size() << "\n";
    for (size_t i = 0; i < m.tri.size(); i++) out << "5\n";
    if (!cellData.empty()) out << "CELL_DATA " << m.tri.size() << "\n";
    for (const auto &[nome, v] : cellData) {
        out << "SCALARS " << nome << " double 1\nLOOKUP_TABLE default\n";
        for (size_t i = 0; i < m.tri.size(); i++) out << v[i] << "\n";
    }
    if (!pointData.empty()) out << "POINT_DATA " << m.nodes.size() << "\n";
    for (const auto &[nome, v] : pointData) {
        out << "SCALARS " << nome << " double 1\nLOOKUP_TABLE default\n";
        for (size_t i = 0; i < m.nodes.size(); i++) out << v[i] << "\n";
    }
}

} // namespace slope
