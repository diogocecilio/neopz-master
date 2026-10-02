// TPZPoroCamClayUP.cpp — ver TPZPoroCamClayUP.h

#include "TPZPoroCamClayUP.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzquad.h"
#include "pzskylnsymmat.h"
#include "tpzcube.h"
#include "tpzquadrilateral.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticquad.h"
#include "TPZCutHillMcKee.h"

namespace {

const REAL kCorner[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
// arestas do hexaedro na ordem dos nós 8-19 do TPZQuadraticCube (lados 8-19 do TPZCube)
const int kEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                           {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
// VTK_QUADRATIC_HEXAHEDRON: vértices, arestas da face inferior, da superior e verticais
const int kVTK20[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};

/// coordenadas paramétricas do nó a do hexaedro (ordem do NeoPZ)
void NodeXi(int a, std::array<REAL, 3> &xi) {
    if (a < 8) {
        for (int i = 0; i < 3; i++) xi[i] = kCorner[a][i];
    } else {
        for (int i = 0; i < 3; i++) xi[i] = 0.5 * (kCorner[kEdges[a - 8][0]][i] + kCorner[kEdges[a - 8][1]][i]);
    }
}

REAL Det3(const REAL J[3][3]) {
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
           J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

void Inv3(const REAL J[3][3], REAL Ji[3][3]) {
    const REAL d = Det3(J);
    Ji[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / d;
    Ji[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / d;
    Ji[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / d;
    Ji[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / d;
    Ji[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / d;
    Ji[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / d;
    Ji[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / d;
    Ji[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / d;
    Ji[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / d;
}

REAL Norm2(const std::vector<REAL> &v) {
    REAL s = 0.;
    for (REAL x : v) s += x * x;
    return std::sqrt(s);
}

/// Componentes não nulas da coluna (3a + i) de B (Voigt {xx, xy, xz, yy, yz, zz}, distorções de engenharia)
inline void BColumn(int i, const REAL dN[3], int rows[3], REAL vals[3]) {
    switch (i) {
        case 0:
            rows[0] = _XX_; vals[0] = dN[0];
            rows[1] = _XY_; vals[1] = dN[1];
            rows[2] = _XZ_; vals[2] = dN[2];
            break;
        case 1:
            rows[0] = _YY_; vals[0] = dN[1];
            rows[1] = _XY_; vals[1] = dN[0];
            rows[2] = _YZ_; vals[2] = dN[2];
            break;
        default:
            rows[0] = _ZZ_; vals[0] = dN[2];
            rows[1] = _XZ_; vals[1] = dN[0];
            rows[2] = _YZ_; vals[2] = dN[1];
            break;
    }
}

} // namespace

// =====================================================================================================
TPZPoroCamClayUP::TPZPoroCamClayUP(TPZGeoMesh *gmesh, EElement type, const TParFunc &parOfX,
                                   const TCoupling &coupling, const TScalarFunc &p0OfX, bool initAtCentroid)
    : fType(type), fCoupling(coupling) {
    fNu = (type == EHex8) ? 8 : 20;
    const int nfn = (type == EHex8) ? 4 : 8;
    fNN = gmesh->NNodes();
    fX.resize(fNN);
    for (int64_t i = 0; i < fNN; i++) {
        TPZManVector<REAL, 3> c(3);
        gmesh->NodeVec()[i].GetCoordinates(c);
        fX[i] = {c[0], c[1], c[2]};
    }
    for (int64_t iel = 0; iel < gmesh->NElements(); iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel || gel->HasSubElement()) continue;
        if (gel->Dimension() == 3) {
            if (gel->NNodes() != fNu) throw std::invalid_argument("TPZPoroCamClayUP: elemento 3D incompatível");
            TElem e;
            e.nodes.resize(fNu);
            for (int a = 0; a < fNu; a++) e.nodes[a] = gel->NodeIndex(a);
            fEl.push_back(std::move(e));
        } else if (gel->Dimension() == 2) {
            if (gel->NNodes() != nfn) throw std::invalid_argument("TPZPoroCamClayUP: face incompatível");
            TFace f;
            f.marker = gel->MaterialId();
            for (int a = 0; a < nfn; a++) f.nodes.push_back(gel->NodeIndex(a));
            fFaces.push_back(f);
        }
    }
    // graus de liberdade de pressão: vértices dos elementos (em ordem crescente de nó)
    std::set<int64_t> pnodes;
    for (auto &e : fEl)
        for (int a = 0; a < 8; a++) pnodes.insert(e.nodes[a]);
    fPMap.assign(fNN, -1);
    for (int64_t v : pnodes) {
        fPMap[v] = int64_t(fPNode.size());
        fPNode.push_back(v);
    }
    fNP = int64_t(fPNode.size());

    // regra de integração de volume (Gauss-Legendre n x n x n)
    const int ng = (type == EHex20) ? 3 : 2;
    TPZIntCube3D rule(2 * ng - 1);
    if (rule.NPoints() != ng * ng * ng) throw std::logic_error("TPZPoroCamClayUP: regra de integração inesperada");
    {
        std::set<REAL> g;
        TPZManVector<REAL, 3> xi(3);
        REAL w;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            rule.Point(ip, xi, w);
            bool novo = true;
            for (REAL v : g)
                if (std::fabs(v - xi[0]) < 1.e-12) novo = false;
            if (novo) g.insert(xi[0]);
        }
        fGauss1D.assign(g.begin(), g.end());
    }

    const REAL alpha = fCoupling.alpha, invM = fCoupling.invBiotModulus, perm = fCoupling.perm;
    fFg.assign(fNP, 0.);
    fFb.assign(3 * fNN, 0.);
    fState.resize(fEl.size());
    for (size_t iel = 0; iel < fEl.size(); iel++) {
        TElem &e = fEl[iel];
        e.Qe.Redim(3 * fNu, 8);
        e.Se.Redim(8, 8);
        e.He.Redim(8, 8);
        int centroidPar = -1;
        if (initAtCentroid) {
            TPZManVector<REAL, 3> c(3, 0.);
            for (int a = 0; a < 8; a++)
                for (int i = 0; i < 3; i++) c[i] += fX[e.nodes[a]][i] / 8.;
            fPar.push_back(parOfX(c));
            centroidPar = int(fPar.size()) - 1;
        }
        TPZFNMatrix<20, REAL> phi(fNu, 1);
        TPZFNMatrix<60, REAL> dphi(3, fNu);
        TPZFNMatrix<8, REAL> phip(8, 1);
        TPZFNMatrix<24, REAL> dphip(3, 8);
        TPZManVector<REAL, 3> xi(3);
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            REAL w;
            rule.Point(ip, xi, w);
            ShapeU(xi, phi, dphi);
            pztopology::TPZCube::Shape(xi, phip, dphip);
            REAL J[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}, Ji[3][3];
            for (int a = 0; a < fNu; a++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) J[i][j] += dphi(i, a) * fX[e.nodes[a]][j];
            const REAL detJ = Det3(J);
            if (detJ <= 0.) throw std::invalid_argument("TPZPoroCamClayUP: jacobiano não positivo");
            Inv3(J, Ji);
            TGP gp;
            gp.Nu.resize(fNu);
            gp.dNu.Redim(3, fNu);
            gp.x = {0., 0., 0.};
            gp.xi = {xi[0], xi[1], xi[2]};
            for (int a = 0; a < fNu; a++) {
                gp.Nu[a] = phi(a, 0);
                for (int i = 0; i < 3; i++) {
                    REAL s = 0.;
                    for (int k = 0; k < 3; k++) s += Ji[i][k] * dphi(k, a);
                    gp.dNu(i, a) = s;
                    gp.x[i] += phi(a, 0) * fX[e.nodes[a]][i];
                }
            }
            for (int b = 0; b < 8; b++) {
                gp.Np[b] = phip(b, 0);
                for (int i = 0; i < 3; i++) {
                    REAL s = 0.;
                    for (int k = 0; k < 3; k++) s += Ji[i][k] * dphip(k, b);
                    gp.dNp[i][b] = s;
                }
            }
            gp.wdJ = w * detJ;
            // parâmetros e estado inicial
            if (initAtCentroid) {
                e.par.push_back(centroidPar);
            } else {
                TPZManVector<REAL, 3> x = {gp.x[0], gp.x[1], gp.x[2]};
                fPar.push_back(parOfX(x));
                e.par.push_back(int(fPar.size()) - 1);
            }
            TGPState st;
            st.sig = fPar[e.par.back()].InitialStress();
            fState[iel].push_back(st);
            // acoplamento, fluxo, força de corpo
            for (int a = 0; a < fNu; a++)
                for (int c = 0; c < 3; c++)
                    for (int b = 0; b < 8; b++) e.Qe(3 * a + c, b) += alpha * gp.dNu(c, a) * gp.Np[b] * gp.wdJ;
            for (int a = 0; a < 8; a++) {
                for (int b = 0; b < 8; b++) {
                    e.Se(a, b) += invM * gp.Np[a] * gp.Np[b] * gp.wdJ;
                    e.He(a, b) += perm * (gp.dNp[0][a] * gp.dNp[0][b] + gp.dNp[1][a] * gp.dNp[1][b] +
                                          gp.dNp[2][a] * gp.dNp[2][b]) * gp.wdJ;
                }
                fFg[fPMap[e.nodes[a]]] += perm * gp.dNp[2][a] * (-fCoupling.gammaW) * gp.wdJ;
            }
            for (int a = 0; a < fNu; a++)
                for (int c = 0; c < 3; c++) fFb[3 * e.nodes[a] + c] += gp.Nu[a] * fCoupling.body[c] * gp.wdJ;
            e.gps.push_back(std::move(gp));
        }
    }
    fU.assign(3 * fNN, 0.);
    fP.assign(fNP, 0.);
    for (int64_t k = 0; k < fNP; k++) {
        const auto &c = fX[fPNode[k]];
        TPZManVector<REAL, 3> x = {c[0], c[1], c[2]};
        fP[k] = p0OfX ? p0OfX(x) : 0.;
    }
    fFint = InternalForces();  // σ = σ'0 em todos os pontos

    // renumeração dos nós (Cuthill-McKee reverso) para reduzir o perfil do sistema
    TPZManVector<int64_t> elgraph(int64_t(fEl.size()) * fNu), elgraphindex(fEl.size() + 1);
    elgraphindex[0] = 0;
    for (size_t iel = 0; iel < fEl.size(); iel++) {
        for (int a = 0; a < fNu; a++) elgraph[iel * fNu + a] = fEl[iel].nodes[a];
        elgraphindex[iel + 1] = int64_t(iel + 1) * fNu;
    }
    TPZCutHillMcKee renum(int64_t(fEl.size()), fNN, true);
    renum.fVerbose = false;
    renum.SetElementGraph(elgraph, elgraphindex);
    TPZManVector<int64_t> nodeperm, inodeperm;
    renum.Resequence(nodeperm, inodeperm);
    fNodePerm.resize(fNN);
    for (int64_t i = 0; i < fNN; i++) fNodePerm[i] = nodeperm[i];
}

void TPZPoroCamClayUP::ShapeU(TPZVec<REAL> &xi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) const {
    if (fType == EHex8) {
        pztopology::TPZCube::Shape(xi, phi, dphi);
    } else {
        pzgeom::TPZQuadraticCube::Shape(xi, phi, dphi);
    }
}

std::vector<int64_t> TPZPoroCamClayUP::NodesWhere(const std::function<bool(const std::array<REAL, 3> &)> &crit) const {
    std::vector<int64_t> out;
    for (int64_t i = 0; i < fNN; i++)
        if (crit(fX[i])) out.push_back(i);
    return out;
}

std::vector<int64_t> TPZPoroCamClayUP::NodesOnMarker(int marker) const {
    std::set<int64_t> s;
    for (auto &f : fFaces)
        if (f.marker == marker) s.insert(f.nodes.begin(), f.nodes.end());
    return std::vector<int64_t>(s.begin(), s.end());
}

// ------------------------------------------------------------------------------------------- cargas
std::vector<REAL> TPZPoroCamClayUP::FaceLoad(const TFaceSelect &select, const std::array<REAL, 3> &t) const {
    std::vector<REAL> F(3 * fNN, 0.);
    const int ng = (fType == EHex8) ? 2 : 3;
    TPZIntQuad rule(2 * ng - 1);
    const int nfn = (fType == EHex8) ? 4 : 8;
    TPZFNMatrix<8, REAL> phi(nfn, 1);
    TPZFNMatrix<16, REAL> dphi(2, nfn);
    TPZManVector<REAL, 2> xi(2);
    for (auto &f : fFaces) {
        if (!select(f)) continue;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            REAL w;
            rule.Point(ip, xi, w);
            if (nfn == 4) pztopology::TPZQuadrilateral::Shape(xi, phi, dphi);
            else pzgeom::TPZQuadraticQuad::Shape(xi, phi, dphi);
            REAL J[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
            for (int a = 0; a < nfn; a++)
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 3; j++) J[i][j] += dphi(i, a) * fX[f.nodes[a]][j];
            const REAL n[3] = {J[0][1] * J[1][2] - J[0][2] * J[1][1], J[0][2] * J[1][0] - J[0][0] * J[1][2],
                               J[0][0] * J[1][1] - J[0][1] * J[1][0]};
            const REAL dA = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            for (int a = 0; a < nfn; a++)
                for (int c = 0; c < 3; c++) F[3 * f.nodes[a] + c] += phi(a, 0) * t[c] * w * dA;
        }
    }
    return F;
}

std::vector<REAL> TPZPoroCamClayUP::FacePressure(const TFaceSelect &select, REAL pressure,
                                                 const TVecFunc &outward) const {
    std::vector<REAL> F(3 * fNN, 0.);
    const int ng = (fType == EHex8) ? 2 : 3;
    TPZIntQuad rule(2 * ng - 1);
    const int nfn = (fType == EHex8) ? 4 : 8;
    TPZFNMatrix<8, REAL> phi(nfn, 1);
    TPZFNMatrix<16, REAL> dphi(2, nfn);
    TPZManVector<REAL, 2> xi(2);
    for (auto &f : fFaces) {
        if (!select(f)) continue;
        for (int ip = 0; ip < rule.NPoints(); ip++) {
            REAL w;
            rule.Point(ip, xi, w);
            if (nfn == 4) pztopology::TPZQuadrilateral::Shape(xi, phi, dphi);
            else pzgeom::TPZQuadraticQuad::Shape(xi, phi, dphi);
            REAL J[2][3] = {{0., 0., 0.}, {0., 0., 0.}};
            TPZManVector<REAL, 3> x(3, 0.);
            for (int a = 0; a < nfn; a++)
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 2; i++) J[i][j] += dphi(i, a) * fX[f.nodes[a]][j];
                    x[j] += phi(a, 0) * fX[f.nodes[a]][j];
                }
            REAL n[3] = {J[0][1] * J[1][2] - J[0][2] * J[1][1], J[0][2] * J[1][0] - J[0][0] * J[1][2],
                         J[0][0] * J[1][1] - J[0][1] * J[1][0]};
            const auto o = outward(x);
            if (n[0] * o[0] + n[1] * o[1] + n[2] * o[2] < 0.)
                for (REAL &v : n) v = -v;
            for (int a = 0; a < nfn; a++)
                for (int c = 0; c < 3; c++) F[3 * f.nodes[a] + c] -= pressure * phi(a, 0) * n[c] * w;
        }
    }
    return F;
}

// ------------------------------------------------------------------------------------------- solução
void TPZPoroCamClayUP::BuildEquations(const std::map<int64_t, REAL> &fixedU, const std::map<int64_t, REAL> &fixedP) {
    std::vector<int64_t> keyU, keyP;
    for (auto &kv : fixedU) keyU.push_back(kv.first);
    for (auto &kv : fixedP) keyP.push_back(kv.first);
    if (!fEqU.empty() && keyU == fFixedKeyU && keyP == fFixedKeyP) return;
    fFixedKeyU = keyU;
    fFixedKeyP = keyP;
    fEqU.assign(3 * fNN, -1);
    fEqP.assign(fNP, -1);
    std::vector<int64_t> order(fNN);
    for (int64_t i = 0; i < fNN; i++) order[fNodePerm[i]] = i;
    fNEq = 0;
    for (int64_t node : order) {
        for (int c = 0; c < 3; c++)
            if (!fixedU.count(3 * node + c)) fEqU[3 * node + c] = fNEq++;
        const int64_t pd = fPMap[node];
        if (pd >= 0 && !fixedP.count(pd)) fEqP[pd] = fNEq++;
    }
    fSkyline.resize(fNEq);
    for (int64_t i = 0; i < fNEq; i++) fSkyline[i] = i;
    std::vector<int64_t> eqs;
    for (auto &e : fEl) {
        eqs.clear();
        for (int a = 0; a < fNu; a++)
            for (int c = 0; c < 3; c++)
                if (fEqU[3 * e.nodes[a] + c] >= 0) eqs.push_back(fEqU[3 * e.nodes[a] + c]);
        for (int b = 0; b < 8; b++)
            if (fEqP[fPMap[e.nodes[b]]] >= 0) eqs.push_back(fEqP[fPMap[e.nodes[b]]]);
        if (eqs.empty()) continue;
        const int64_t mn = *std::min_element(eqs.begin(), eqs.end());
        for (int64_t q : eqs) fSkyline[q] = std::min(fSkyline[q], mn);
    }
}

void TPZPoroCamClayUP::Assemble(const std::vector<REAL> &U, std::vector<REAL> &F,
                                std::vector<std::vector<TGPState>> &trial, std::vector<TPZFMatrix<REAL>> *Ke) const {
    F.assign(3 * fNN, 0.);
    trial.resize(fEl.size());
    if (Ke) Ke->resize(fEl.size());
    const int nd = 3 * fNu;
    std::vector<REAL> ue(nd);
    TPZModifiedCamClay::TResult r;
    for (size_t iel = 0; iel < fEl.size(); iel++) {
        const TElem &e = fEl[iel];
        for (int a = 0; a < fNu; a++)
            for (int c = 0; c < 3; c++) ue[3 * a + c] = U[3 * e.nodes[a] + c];
        if (Ke) {
            (*Ke)[iel].Redim(nd, nd);
        }
        trial[iel].resize(e.gps.size());
        for (size_t ip = 0; ip < e.gps.size(); ip++) {
            const TGP &g = e.gps[ip];
            TPZTensor<REAL> eps;
            for (int a = 0; a < fNu; a++) {
                const REAL d0 = g.dNu.GetVal(0, a), d1 = g.dNu.GetVal(1, a), d2 = g.dNu.GetVal(2, a);
                const REAL ux = ue[3 * a], uy = ue[3 * a + 1], uz = ue[3 * a + 2];
                eps[_XX_] += d0 * ux;
                eps[_YY_] += d1 * uy;
                eps[_ZZ_] += d2 * uz;
                eps[_XY_] += d1 * ux + d0 * uy;
                eps[_XZ_] += d2 * ux + d0 * uz;
                eps[_YZ_] += d2 * uy + d1 * uz;
            }
            const TGPState &sn = fState[iel][ip];
            fPar[e.par[ip]].ReturnMapping(eps, sn.epsp, sn.alpha, r, &sn.sig, &sn.eps);
            TGPState &st = trial[iel][ip];
            st.epsp = r.plastic_strain;
            st.sig = r.stress;
            st.eps = eps;
            st.alpha = r.alpha;
            st.plastic = r.plastic;
            const TPZTensor<REAL> &s = r.stress;
            for (int a = 0; a < fNu; a++) {
                const REAL d0 = g.dNu.GetVal(0, a), d1 = g.dNu.GetVal(1, a), d2 = g.dNu.GetVal(2, a);
                const int64_t n = e.nodes[a];
                F[3 * n] += (s[_XX_] * d0 + s[_XY_] * d1 + s[_XZ_] * d2) * g.wdJ;
                F[3 * n + 1] += (s[_XY_] * d0 + s[_YY_] * d1 + s[_YZ_] * d2) * g.wdJ;
                F[3 * n + 2] += (s[_XZ_] * d0 + s[_YZ_] * d1 + s[_ZZ_] * d2) * g.wdJ;
            }
            if (!Ke) continue;
            // K += Bᵀ D B wdJ (B tem 3 componentes não nulas por coluna)
            TPZFMatrix<REAL> &K = (*Ke)[iel];
            std::vector<std::array<REAL, 6>> DB(nd);
            std::vector<std::array<int, 3>> brow(nd);
            std::vector<std::array<REAL, 3>> bval(nd);
            for (int a = 0; a < fNu; a++) {
                const REAL dN[3] = {g.dNu.GetVal(0, a), g.dNu.GetVal(1, a), g.dNu.GetVal(2, a)};
                for (int i = 0; i < 3; i++) {
                    const int col = 3 * a + i;
                    BColumn(i, dN, brow[col].data(), bval[col].data());
                    for (int rr = 0; rr < 6; rr++) {
                        REAL v = 0.;
                        for (int k = 0; k < 3; k++) v += r.Dep(rr, brow[col][k]) * bval[col][k];
                        DB[col][rr] = v;
                    }
                }
            }
            for (int row = 0; row < nd; row++)
                for (int col = 0; col < nd; col++) {
                    REAL v = 0.;
                    for (int k = 0; k < 3; k++) v += bval[row][k] * DB[col][brow[row][k]];
                    K(row, col) += v * g.wdJ;
                }
        }
    }
}

std::vector<REAL> TPZPoroCamClayUP::InternalForces() const {
    std::vector<REAL> F;
    std::vector<std::vector<TGPState>> trial;
    Assemble(fU, F, trial, nullptr);
    return F;
}

int TPZPoroCamClayUP::Step(const std::vector<REAL> &Fext, const std::map<int64_t, REAL> &fixedU,
                           const std::map<int64_t, REAL> &fixedP, REAL dt, bool flow, REAL tol, int maxit,
                           const std::vector<REAL> *Uguess, bool verbose) {
    BuildEquations(fixedU, fixedP);
    std::vector<REAL> U = Uguess ? *Uguess : fU;
    std::vector<REAL> P = fP;
    for (auto &kv : fixedU) U[kv.first] = kv.second;
    for (auto &kv : fixedP) P[kv.first] = kv.second;
    const REAL ref = std::max(Norm2(Fext), REAL(1.));
    const REAL hfac = flow ? dt : 0.;
    std::vector<REAL> F, Ru, Rp;
    std::vector<std::vector<TGPState>> trial;
    std::vector<TPZFMatrix<REAL>> Ke;
    bool conv = false;
    REAL nr = 0.;
    int it;
    for (it = 1; it <= maxit; it++) {
        Assemble(U, F, trial, &Ke);
        // resíduos
        Ru.assign(3 * fNN, 0.);
        Rp.assign(fNP, 0.);
        for (int64_t i = 0; i < 3 * fNN; i++) Ru[i] = F[i] - Fext[i];
        for (auto &e : fEl) {
            REAL pe[8], dpe[8];
            for (int b = 0; b < 8; b++) {
                const int64_t k = fPMap[e.nodes[b]];
                pe[b] = P[k];
                dpe[b] = P[k] - fP[k];
            }
            for (int a = 0; a < fNu; a++)
                for (int c = 0; c < 3; c++) {
                    const int64_t gi = 3 * e.nodes[a] + c;
                    const REAL du = U[gi] - fU[gi];
                    for (int b = 0; b < 8; b++) {
                        Ru[gi] -= e.Qe.GetVal(3 * a + c, b) * pe[b];
                        Rp[fPMap[e.nodes[b]]] += e.Qe.GetVal(3 * a + c, b) * du;
                    }
                }
            for (int b = 0; b < 8; b++) {
                REAL v = 0.;
                for (int c = 0; c < 8; c++) v += e.Se.GetVal(b, c) * dpe[c] + hfac * e.He.GetVal(b, c) * pe[c];
                Rp[fPMap[e.nodes[b]]] += v;
            }
        }
        if (flow)
            for (int64_t k = 0; k < fNP; k++) Rp[k] -= dt * fFg[k];
        REAL s = 0.;
        for (int64_t i = 0; i < 3 * fNN; i++)
            if (fEqU[i] >= 0) s += Ru[i] * Ru[i];
        for (int64_t k = 0; k < fNP; k++)
            if (fEqP[k] >= 0) s += Rp[k] * Rp[k];
        nr = std::sqrt(s) / ref;
        if (verbose) std::cout << "      it " << it << ": |R|/|F| = " << nr << "\n";
        if (nr < tol) {
            conv = true;
            break;
        }
        // jacobiana [[K, -Q], [Qᵀ, S + Δt H]] nos gdl livres
        TPZManVector<int64_t> sky(fNEq);
        for (int64_t i = 0; i < fNEq; i++) sky[i] = fSkyline[i];
        TPZSkylNSymMatrix<STATE> J(fNEq, sky);
        std::vector<int64_t> eu(3 * fNu), ep(8);
        for (size_t iel = 0; iel < fEl.size(); iel++) {
            const TElem &e = fEl[iel];
            for (int a = 0; a < fNu; a++)
                for (int c = 0; c < 3; c++) eu[3 * a + c] = fEqU[3 * e.nodes[a] + c];
            for (int b = 0; b < 8; b++) ep[b] = fEqP[fPMap[e.nodes[b]]];
            const TPZFMatrix<REAL> &K = Ke[iel];
            for (int i = 0; i < 3 * fNu; i++) {
                if (eu[i] < 0) continue;
                for (int j = 0; j < 3 * fNu; j++)
                    if (eu[j] >= 0) J(eu[i], eu[j]) += K.GetVal(i, j);
                for (int b = 0; b < 8; b++)
                    if (ep[b] >= 0) {
                        J(eu[i], ep[b]) -= e.Qe.GetVal(i, b);
                        J(ep[b], eu[i]) += e.Qe.GetVal(i, b);
                    }
            }
            for (int b = 0; b < 8; b++) {
                if (ep[b] < 0) continue;
                for (int c = 0; c < 8; c++)
                    if (ep[c] >= 0) J(ep[b], ep[c]) += e.Se.GetVal(b, c) + hfac * e.He.GetVal(b, c);
            }
        }
        TPZFMatrix<STATE> rhs(fNEq, 1, 0.);
        for (int64_t i = 0; i < 3 * fNN; i++)
            if (fEqU[i] >= 0) rhs(fEqU[i], 0) = -Ru[i];
        for (int64_t k = 0; k < fNP; k++)
            if (fEqP[k] >= 0) rhs(fEqP[k], 0) = -Rp[k];
        if (!J.Decompose_LU()) throw std::runtime_error("TPZPoroCamClayUP::Step: fatoração LU falhou");
        J.Subst_LForward(&rhs);
        J.Subst_Backward(&rhs);
        for (int64_t i = 0; i < 3 * fNN; i++)
            if (fEqU[i] >= 0) U[i] += rhs(fEqU[i], 0);
        for (int64_t k = 0; k < fNP; k++)
            if (fEqP[k] >= 0) P[k] += rhs(fEqP[k], 0);
    }
    if (!conv) {
        std::stringstream sout;
        sout << "TPZPoroCamClayUP::Step: Newton não convergiu (|R|/|F| = " << nr << ")";
        throw std::runtime_error(sout.str());
    }
    fU = U;
    fP = P;
    fState = trial;
    fFint = F;
    return it;
}

std::vector<REAL> TPZPoroCamClayUP::Reactions(const std::vector<REAL> &Fext) const {
    std::vector<REAL> R(3 * fNN);
    for (int64_t i = 0; i < 3 * fNN; i++) R[i] = fFint[i] - Fext[i];
    for (auto &e : fEl)
        for (int a = 0; a < fNu; a++)
            for (int c = 0; c < 3; c++)
                for (int b = 0; b < 8; b++)
                    R[3 * e.nodes[a] + c] -= e.Qe.GetVal(3 * a + c, b) * fP[fPMap[e.nodes[b]]];
    return R;
}

// ------------------------------------------------------------------------------------------- pós-processamento
bool TPZPoroCamClayUP::Locate(const std::array<REAL, 3> &x, int64_t &elout, std::array<REAL, 3> &xiout) const {
    TPZFNMatrix<20, REAL> phi(fNu, 1);
    TPZFNMatrix<60, REAL> dphi(3, fNu);
    TPZManVector<REAL, 3> xi(3);
    for (size_t iel = 0; iel < fEl.size(); iel++) {
        const TElem &e = fEl[iel];
        bool fora = false;
        for (int i = 0; i < 3 && !fora; i++) {
            REAL mn = 1.e300, mx = -1.e300;
            for (int a = 0; a < fNu; a++) {
                mn = std::min(mn, fX[e.nodes[a]][i]);
                mx = std::max(mx, fX[e.nodes[a]][i]);
            }
            if (x[i] < mn - 1.e-12 || x[i] > mx + 1.e-12) fora = true;
        }
        if (fora) continue;
        xi[0] = xi[1] = xi[2] = 0.;
        for (int k = 0; k < 30; k++) {
            ShapeU(xi, phi, dphi);
            REAL r[3] = {-x[0], -x[1], -x[2]}, J[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}, Ji[3][3];
            for (int a = 0; a < fNu; a++)
                for (int j = 0; j < 3; j++) {
                    r[j] += phi(a, 0) * fX[e.nodes[a]][j];
                    for (int i = 0; i < 3; i++) J[j][i] += dphi(i, a) * fX[e.nodes[a]][j];  // ∂x_j/∂ξ_i
                }
            if (std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < 1.e-14) break;
            Inv3(J, Ji);
            for (int i = 0; i < 3; i++) xi[i] -= Ji[i][0] * r[0] + Ji[i][1] * r[1] + Ji[i][2] * r[2];
        }
        if (std::fabs(xi[0]) <= 1. + 1.e-7 && std::fabs(xi[1]) <= 1. + 1.e-7 && std::fabs(xi[2]) <= 1. + 1.e-7) {
            elout = int64_t(iel);
            xiout = {xi[0], xi[1], xi[2]};
            return true;
        }
    }
    return false;
}

std::vector<REAL> TPZPoroCamClayUP::GaussWeights(int64_t el, const std::array<REAL, 3> &xi) const {
    const std::vector<REAL> &g = fGauss1D;
    const int n = int(g.size());
    auto lag = [&](REAL t, int l) {
        REAL v = 1.;
        for (int m = 0; m < n; m++)
            if (m != l) v *= (t - g[m]) / (g[l] - g[m]);
        return v;
    };
    auto idx = [&](REAL t) {
        int best = 0;
        for (int m = 1; m < n; m++)
            if (std::fabs(g[m] - t) < std::fabs(g[best] - t)) best = m;
        return best;
    };
    const TElem &e = fEl[el];
    std::vector<REAL> w(e.gps.size());
    for (size_t ip = 0; ip < e.gps.size(); ip++) {
        const auto &p = e.gps[ip].xi;
        w[ip] = lag(xi[0], idx(p[0])) * lag(xi[1], idx(p[1])) * lag(xi[2], idx(p[2]));
    }
    return w;
}

TPZTensor<REAL> TPZPoroCamClayUP::StressAt(int64_t el, const std::array<REAL, 3> &xi) const {
    const std::vector<REAL> w = GaussWeights(el, xi);
    TPZTensor<REAL> s;
    for (size_t ip = 0; ip < w.size(); ip++)
        for (int i = 0; i < 6; i++) s[i] += w[ip] * fState[el][ip].sig[i];
    return s;
}

REAL TPZPoroCamClayUP::ZonePressure(int64_t el) const {
    REAL s = 0.;
    for (int b = 0; b < 8; b++) s += fP[fPMap[fEl[el].nodes[b]]];
    return s / 8.;
}

std::vector<REAL> TPZPoroCamClayUP::NodalPressure() const {
    std::vector<REAL> p(fNN, 0.);
    for (int64_t k = 0; k < fNP; k++) p[fPNode[k]] = fP[k];
    if (fNu == 20)
        for (auto &e : fEl)
            for (int a = 8; a < 20; a++)
                p[e.nodes[a]] = 0.5 * (fP[fPMap[e.nodes[kEdges[a - 8][0]]]] + fP[fPMap[e.nodes[kEdges[a - 8][1]]]]);
    return p;
}

void TPZPoroCamClayUP::WriteVTK(const std::string &file, const std::string &title,
                                const std::vector<std::pair<std::string, REAL>> &fieldData,
                                const TScalarFunc &hydro) const {
    std::ofstream out(file);
    out << std::setprecision(10);
    out << "# vtk DataFile Version 3.0\n" << title << "\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    if (!fieldData.empty()) {
        out << "FIELD FieldData " << fieldData.size() << "\n";
        for (auto &fd : fieldData) out << fd.first << " 1 1 double\n" << fd.second << "\n";
    }
    out << "POINTS " << fNN << " double\n";
    for (auto &x : fX) out << x[0] << " " << x[1] << " " << x[2] << "\n";
    const int64_t nel = int64_t(fEl.size());
    out << "CELLS " << nel << " " << nel * (fNu + 1) << "\n";
    for (auto &e : fEl) {
        out << fNu;
        for (int a = 0; a < fNu; a++) out << " " << e.nodes[fNu == 20 ? kVTK20[a] : a];
        out << "\n";
    }
    out << "CELL_TYPES " << nel << "\n";
    for (int64_t i = 0; i < nel; i++) out << (fNu == 20 ? 25 : 12) << "\n";

    // ---- campos nodais
    const std::vector<REAL> pn = NodalPressure();
    std::vector<TPZTensor<REAL>> sn(fNN);
    std::vector<int> cnt(fNN, 0);
    for (int64_t iel = 0; iel < nel; iel++) {
        const TElem &e = fEl[iel];
        for (int a = 0; a < fNu; a++) {
            std::array<REAL, 3> xi;
            NodeXi(a, xi);
            const TPZTensor<REAL> s = StressAt(iel, xi);  // extrapolação dos pontos de Gauss
            sn[e.nodes[a]] += s;
            cnt[e.nodes[a]]++;
        }
    }
    for (int64_t i = 0; i < fNN; i++)
        if (cnt[i]) sn[i] *= 1. / cnt[i];
    auto tensor = [&out](const TPZTensor<REAL> &s) {
        out << s[_XX_] << " " << s[_XY_] << " " << s[_XZ_] << "\n"
            << s[_XY_] << " " << s[_YY_] << " " << s[_YZ_] << "\n"
            << s[_XZ_] << " " << s[_YZ_] << " " << s[_ZZ_] << "\n\n";
    };
    out << "POINT_DATA " << fNN << "\n";
    out << "VECTORS deslocamento double\n";
    for (int64_t i = 0; i < fNN; i++) out << fU[3 * i] << " " << fU[3 * i + 1] << " " << fU[3 * i + 2] << "\n";
    out << "SCALARS poropressao double 1\nLOOKUP_TABLE default\n";
    for (int64_t i = 0; i < fNN; i++) out << pn[i] << "\n";
    if (hydro) {
        out << "SCALARS excesso_poropressao double 1\nLOOKUP_TABLE default\n";
        for (int64_t i = 0; i < fNN; i++) {
            TPZManVector<REAL, 3> x = {fX[i][0], fX[i][1], fX[i][2]};
            out << pn[i] - hydro(x) << "\n";
        }
    }
    out << "TENSORS tensao_efetiva double\n";
    for (int64_t i = 0; i < fNN; i++) tensor(sn[i]);
    out << "TENSORS tensao_total double\n";
    for (int64_t i = 0; i < fNN; i++) {
        TPZTensor<REAL> st(sn[i]);
        for (int c : {_XX_, _YY_, _ZZ_}) st[c] -= fCoupling.alpha * pn[i];
        tensor(st);
    }
    out << "SCALARS p_efetiva double 1\nLOOKUP_TABLE default\n";
    for (int64_t i = 0; i < fNN; i++) {
        REAL p, q;
        TPZModifiedCamClay::Invariants(sn[i], p, q);
        out << -p << "\n";
    }
    out << "SCALARS q double 1\nLOOKUP_TABLE default\n";
    for (int64_t i = 0; i < fNN; i++) {
        REAL p, q;
        TPZModifiedCamClay::Invariants(sn[i], p, q);
        out << q << "\n";
    }

    // ---- campos por elemento (médias dos pontos de Gauss)
    out << "CELL_DATA " << nel << "\n";
    std::vector<TPZTensor<REAL>> se(nel);
    std::vector<REAL> al(nel, 0.), pc(nel, 0.), fpl(nel, 0.), fyi(nel, 0.);
    for (int64_t iel = 0; iel < nel; iel++) {
        const int ngp = int(fState[iel].size());
        for (int ip = 0; ip < ngp; ip++) {
            const TGPState &st = fState[iel][ip];
            se[iel] += st.sig;
            al[iel] += st.alpha / ngp;
            pc[iel] += fPar[fEl[iel].par[ip]].Pc(st.alpha) / ngp;
            fpl[iel] += (st.plastic ? 1. : 0.) / ngp;
            fyi[iel] += (st.alpha != 0. ? 1. : 0.) / ngp;
        }
        se[iel] *= 1. / ngp;
    }
    out << "TENSORS tensao_efetiva_media double\n";
    for (int64_t i = 0; i < nel; i++) tensor(se[i]);
    std::vector<std::pair<std::string, std::vector<REAL>>> cs;
    std::vector<REAL> pm(nel), qm(nel), zp(nel);
    for (int64_t i = 0; i < nel; i++) {
        REAL p, q;
        TPZModifiedCamClay::Invariants(se[i], p, q);
        pm[i] = -p;
        qm[i] = q;
        zp[i] = ZonePressure(i);
    }
    cs.push_back({"p_efetiva_media", pm});
    cs.push_back({"q_medio", qm});
    cs.push_back({"alpha_medio", al});
    cs.push_back({"pc_medio", pc});
    cs.push_back({"fracao_plastica_passo", fpl});
    cs.push_back({"fracao_plastificada", fyi});
    cs.push_back({"poropressao_zona", zp});
    for (auto &c : cs) {
        out << "SCALARS " << c.first << " double 1\nLOOKUP_TABLE default\n";
        for (REAL v : c.second) out << v << "\n";
    }
}
