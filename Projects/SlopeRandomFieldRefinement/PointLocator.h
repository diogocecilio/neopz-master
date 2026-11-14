// Interpolacao de um campo nodal de uma malha em pontos quaisquer: acha o elemento que contem o ponto (triangulo P1 ou
// quadrilatero Q1) e combina os valores dos vertices com as funcoes de forma (coordenadas baricentricas no triangulo;
// bilineares no quadrilatero, com a inversa do mapa geometrico por Newton). Serve para levar a p_w da malha do transiente
// para os nos de outra malha (a do fator de seguranca). Usa as folhas de dimensao 2 da malha e uma grade de baldes (cada
// balde guarda os elementos que o tocam), para nao testar todos os elementos em cada ponto
#pragma once

#include "pzerror.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

#include <algorithm>
#include <cmath>
#include <vector>

class PointLocator {
public:
    explicit PointLocator(TPZGeoMesh *g, int nx = 64, int ny = 64) : fnx(nx), fny(ny)
    {
        REAL x0 = 1.e300, x1 = -1.e300, y0 = 1.e300, y1 = -1.e300;
        for (TPZGeoEl *gel : g->ElementVec()) {
            if (!gel || gel->HasSubElement() || gel->Dimension() != 2) continue;
            El e;
            e.gel = gel;
            e.nv = gel->NCornerNodes();
            if (e.nv != 3 && e.nv != 4) DebugStop();
            for (int i = 0; i < e.nv; i++) {
                e.n[i] = gel->NodeIndex(i);
                e.x[i] = gel->NodePtr(i)->Coord(0), e.y[i] = gel->NodePtr(i)->Coord(1);
                x0 = std::min(x0, e.x[i]), x1 = std::max(x1, e.x[i]), y0 = std::min(y0, e.y[i]), y1 = std::max(y1, e.y[i]);
            }
            fEl.push_back(e);
        }
        fx0 = x0, fy0 = y0, fdx = (x1 - x0) / nx, fdy = (y1 - y0) / ny;
        fBaldes.assign(nx * ny, {});
        for (size_t k = 0; k < fEl.size(); k++) { // o elemento vai para os baldes da sua caixa envolvente
            const El &e = fEl[k];
            const int i0 = Coluna(*std::min_element(e.x, e.x + e.nv)), i1 = Coluna(*std::max_element(e.x, e.x + e.nv));
            const int j0 = Linha(*std::min_element(e.y, e.y + e.nv)), j1 = Linha(*std::max_element(e.y, e.y + e.nv));
            for (int i = i0; i <= i1; i++)
                for (int j = j0; j <= j1; j++) fBaldes[j * fnx + i].push_back(k);
        }
    }

    // Valor de v (um por no da malha) no ponto (x, y). Um ponto um pouco fora de todos os elementos (arredondamento no
    // contorno) usa o elemento do balde em que ele fica menos fora
    REAL Interpolate(REAL x, REAL y, const std::vector<REAL> &v) const
    {
        REAL N[4];
        const int k = Locate(x, y, N);
        if (k < 0) return 0.; // ponto fora da caixa da malha
        const El &e = fEl[k];
        REAL val = 0.;
        for (int i = 0; i < e.nv; i++) val += N[i] * v[e.n[i]];
        return val;
    }

    // Elemento que contem o ponto (x, y) (ou o menos fora dele, como na Interpolate) e as coordenadas qsi do ponto no
    // elemento de referencia do neopz: triangulo (N1, N2), quadrilatero (xi, eta). nullptr: fora da caixa da malha
    TPZGeoEl *Find(REAL x, REAL y, TPZVec<REAL> &qsi) const
    {
        REAL N[4];
        const int k = Locate(x, y, N);
        if (k < 0) return nullptr;
        const El &e = fEl[k];
        qsi.Resize(2);
        if (e.nv == 3) qsi[0] = N[1], qsi[1] = N[2];
        else qsi[0] = N[1] + N[2] - N[0] - N[3], qsi[1] = N[2] + N[3] - N[0] - N[1]; // das funcoes bilineares
        return e.gel;
    }

private:
    struct El {
        TPZGeoEl *gel = nullptr;
        int nv = 0; // 3 ou 4 vertices
        int64_t n[4];
        REAL x[4], y[4];
    };
    // Elemento do balde do ponto que o contem (ou, se nenhum, o menos fora) e as funcoes de forma N dos vertices no
    // ponto; -1: ponto fora da caixa da malha
    int Locate(REAL x, REAL y, REAL N[4]) const
    {
        const std::vector<int> &balde = fBaldes[Linha(y) * fnx + Coluna(x)];
        int melhor = -1;
        REAL fora = -1.e300;
        std::fill(N, N + 4, REAL(0));
        for (int k : balde) {
            REAL Nk[4];
            const REAL m = Forma(fEl[k], x, y, Nk); // >= 0: dentro
            if (m > fora) {
                fora = m, melhor = k;
                std::copy(Nk, Nk + 4, N);
            }
            if (m >= 0.) break;
        }
        return melhor;
    }
    // Funcoes de forma dos vertices no ponto e o quanto ele esta dentro do elemento: a menor coordenada baricentrica
    // (triangulo) ou 1 - max(|xi|, |eta|) (quadrilatero); >= 0 dentro
    static REAL Forma(const El &e, REAL x, REAL y, REAL N[4])
    {
        if (e.nv == 3) {
            const REAL det = (e.y[1] - e.y[2]) * (e.x[0] - e.x[2]) + (e.x[2] - e.x[1]) * (e.y[0] - e.y[2]);
            N[0] = ((e.y[1] - e.y[2]) * (x - e.x[2]) + (e.x[2] - e.x[1]) * (y - e.y[2])) / det;
            N[1] = ((e.y[2] - e.y[0]) * (x - e.x[2]) + (e.x[0] - e.x[2]) * (y - e.y[2])) / det;
            N[2] = 1. - N[0] - N[1];
            N[3] = 0.;
            return std::min({N[0], N[1], N[2]});
        }
        // quadrilatero: vertices em (-1,-1), (1,-1), (1,1), (-1,1); Newton para (xi, eta) com x(xi, eta) = (x, y)
        const REAL sx[4] = {-1., 1., 1., -1.}, se[4] = {-1., -1., 1., 1.};
        REAL xi = 0., et = 0.;
        for (int it = 0; it < 30; it++) {
            REAL fx = -x, fy = -y, j00 = 0., j01 = 0., j10 = 0., j11 = 0.;
            for (int i = 0; i < 4; i++) {
                const REAL Ni = 0.25 * (1. + sx[i] * xi) * (1. + se[i] * et);
                const REAL dxi = 0.25 * sx[i] * (1. + se[i] * et), det = 0.25 * se[i] * (1. + sx[i] * xi);
                fx += Ni * e.x[i], fy += Ni * e.y[i];
                j00 += dxi * e.x[i], j01 += det * e.x[i], j10 += dxi * e.y[i], j11 += det * e.y[i];
            }
            const REAL D = j00 * j11 - j01 * j10;
            const REAL dxi = (j11 * fx - j01 * fy) / D, det = (-j10 * fx + j00 * fy) / D;
            xi = std::clamp(xi - dxi, REAL(-5.), REAL(5.)), et = std::clamp(et - det, REAL(-5.), REAL(5.));
            if (std::abs(dxi) + std::abs(det) < 1.e-13) break;
        }
        for (int i = 0; i < 4; i++) N[i] = 0.25 * (1. + sx[i] * xi) * (1. + se[i] * et);
        return 1. - std::max(std::abs(xi), std::abs(et));
    }
    int Coluna(REAL x) const { return std::clamp(int((x - fx0) / fdx), 0, fnx - 1); }
    int Linha(REAL y) const { return std::clamp(int((y - fy0) / fdy), 0, fny - 1); }

    std::vector<El> fEl;
    REAL fx0 = 0., fy0 = 0., fdx = 1., fdy = 1.;
    int fnx, fny;
    std::vector<std::vector<int>> fBaldes;
};
