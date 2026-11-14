#include "LayerGMesh.h"

#include "pzerror.h"
#include "pzgeoel.h"

#include <algorithm>
#include <cmath>

std::vector<REAL> GeometricSteps(REAL L, REAL h, REAL r, REAL hmax)
{
    if (L <= 0. || h <= 0. || r < 1. || hmax < h) DebugStop();
    std::vector<REAL> p;
    REAL soma = 0., s = h;
    while (soma < L * (1. - 1.e-12)) {
        p.push_back(std::min(s, hmax));
        soma += p.back();
        s *= r;
    }
    // sem o ultimo passo, a soma fica abaixo de L: fica com o ultimo se ele mexer menos nos passos
    const REAL semUltimo = soma - p.back();
    if (p.size() > 1 && std::abs(L / semUltimo - 1.) < std::abs(L / soma - 1.)) p.pop_back(), soma = semUltimo;
    for (REAL &x : p) x *= L / soma;
    return p;
}

namespace {
// altura do terreno em x (linear em cada trecho) e o trecho que contem x
REAL Superficie(const LayerGMeshData &d, REAL x, int *trecho = nullptr)
{
    const auto &s = d.superficie;
    int k = 0;
    while (k + 2 < int(s.size()) && x > s[k + 1][0]) k++;
    if (trecho) *trecho = k;
    const REAL t = (x - s[k][0]) / (s[k + 1][0] - s[k][0]);
    return s[k][1] + t * (s[k + 1][1] - s[k][1]);
}
} // namespace

TPZGeoMesh *LayerGMesh(const LayerGMeshData &d)
{
    const auto &s = d.superficie;
    const int nc = int(d.xs.size()) - 1, nl = int(d.eta.size()) - 1; // colunas e camadas
    if (s.size() < 2 || d.idSuperficie.size() + 1 != s.size() || nc < 1 || nl < 1) DebugStop();
    const REAL tol = 1.e-9 * (s.back()[0] - s.front()[0]);
    if (std::abs(d.xs.front() - s.front()[0]) > tol || std::abs(d.xs.back() - s.back()[0]) > tol) DebugStop();
    if (std::abs(d.eta.front()) > 1.e-12 || std::abs(d.eta.back() - 1.) > 1.e-12) DebugStop();
    for (int i = 0; i < nc; i++)
        if (d.xs[i + 1] <= d.xs[i]) DebugStop(); // colunas crescentes
    for (int j = 0; j < nl; j++)
        if (d.eta[j + 1] <= d.eta[j]) DebugStop(); // camadas crescentes
    for (const auto &v : s) // um vertice do terreno fora das colunas cortaria uma celula ao meio
        if (std::none_of(d.xs.begin(), d.xs.end(), [&](REAL x) { return std::abs(x - v[0]) <= tol; })) DebugStop();

    auto *g = new TPZGeoMesh;
    g->SetDimension(2);
    auto no = [nl](int i, int j) { return int64_t(i) * (nl + 1) + j; };
    g->NodeVec().Resize(int64_t(nc + 1) * (nl + 1));
    for (int i = 0; i <= nc; i++) {
        const REAL ys = Superficie(d, d.xs[i]);
        for (int j = 0; j <= nl; j++) {
            TPZManVector<REAL, 3> x = {d.xs[i], ys * d.eta[j], 0.};
            g->NodeVec()[no(i, j)] = TPZGeoNode(no(i, j), x, *g);
        }
    }
    auto dist = [g](int64_t a, int64_t b) {
        const TPZGeoNode &p = g->NodeVec()[a], &q = g->NodeVec()[b];
        return std::hypot(p.Coord(0) - q.Coord(0), p.Coord(1) - q.Coord(1));
    };

    // celulas: vertices a, b, c, e no sentido anti-horario
    int64_t index;
    for (int i = 0; i < nc; i++)
        for (int j = 0; j < nl; j++) {
            const int64_t a = no(i, j), b = no(i + 1, j), c = no(i + 1, j + 1), e = no(i, j + 1);
            if (d.quadrilateros) {
                TPZManVector<int64_t, 4> q = {a, b, c, e};
                g->CreateGeoElement(EQuadrilateral, q, d.idMaterial, index);
                continue;
            }
            const REAL dac = dist(a, c), dbe = dist(b, e);
            const bool porAC = std::abs(dac - dbe) > 1.e-9 * (dac + dbe) ? dac < dbe : (i + j) % 2 == 0;
            TPZManVector<int64_t, 3> t1(3), t2(3);
            if (porAC) t1 = {a, b, c}, t2 = {a, c, e};
            else t1 = {a, b, e}, t2 = {b, c, e};
            g->CreateGeoElement(ETriangle, t1, d.idMaterial, index);
            g->CreateGeoElement(ETriangle, t2, d.idMaterial, index);
        }

    // contornos
    auto linha = [&](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> l = {a, b};
        g->CreateGeoElement(EOned, l, id, index);
    };
    for (int i = 0; i < nc; i++) linha(no(i, 0), no(i + 1, 0), d.idBase);
    for (int j = 0; j < nl; j++) linha(no(nc, j), no(nc, j + 1), d.idDireita);
    for (int i = nc - 1; i >= 0; i--) {
        int k;
        Superficie(d, 0.5 * (d.xs[i] + d.xs[i + 1]), &k);
        linha(no(i + 1, nl), no(i, nl), d.idSuperficie[k]);
    }
    for (int j = nl - 1; j >= 0; j--) linha(no(0, j + 1), no(0, j), d.idEsquerda);
    g->BuildConnectivity();
    return g;
}
