// A mesh, its multiphysics space and the analysis, with the helpers of the column tests (nodal sets, lumped
// masses, stored volumes of water and air).
#ifndef UNSATPROBLEM_H
#define UNSATPROBLEM_H

#include "UnsatAnalysis.h"

#include <memory>

struct UnsatProblem {
    TPZGeoMesh *gmesh = nullptr;
    TPZMultiphysicsCompMesh *mphys = nullptr;
    std::unique_ptr<UnsatAnalysis> an;
    std::map<int64_t, REAL> mass; ///< int N_j dOmega of each vertex (P1)

    UnsatProblem(TPZGeoMesh *g, const std::map<int, UnsatParams> &params,
                 const std::map<int, std::array<bool, 2>> &fixU, const std::vector<Traction> &tractions = {})
        : gmesh(g) {
        mphys = UnsatMesh(gmesh, params, tractions);
        an = std::make_unique<UnsatAnalysis>(mphys, fixU);
        for (int64_t i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *gel = gmesh->Element(i);
            if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
            const TPZGeoNode &a = *gel->NodePtr(0), &b = *gel->NodePtr(1), &c = *gel->NodePtr(2);
            const REAL area = 0.5 * std::fabs((b.Coord(0) - a.Coord(0)) * (c.Coord(1) - a.Coord(1)) -
                                              (c.Coord(0) - a.Coord(0)) * (b.Coord(1) - a.Coord(1)));
            for (int k = 0; k < 3; ++k) mass[gel->NodeIndex(k)] += area / 3.;
        }
    }
    UnsatProblem(const UnsatProblem &) = delete;
    ~UnsatProblem() {
        an.reset();
        if (mphys) {
            TPZVec<TPZCompMesh *> mv = mphys->MeshVector();
            delete mphys;
            for (auto *m : mv) delete m;
        }
        delete gmesh;
    }

    UnsatAnalysis &A() { return *an; }
    REAL Y(int64_t node) const { return gmesh->NodeVec()[node].Coord(1); }
    REAL X(int64_t node) const { return gmesh->NodeVec()[node].Coord(0); }

    /// lumped volume of water n S_w and of air n S_a (per unit thickness), with the vertex pressures of x
    void StoredVolumes(const TPZFMatrix<STATE> &x, const TPZUnsaturatedHydraulics &h, REAL n, REAL &Vw, REAL &Va,
                       REAL *airMass = nullptr, REAL patm = 101.325) const {
        Vw = Va = 0.;
        if (airMass) *airMass = 0.;
        for (auto &m : mass) {
            const REAL pw = an->Nodal(x, m.first, 1), pa = an->Nodal(x, m.first, 2);
            REAL S, G;
            h.Retention(pa - pw, S, G);
            Vw += m.second * n * S;
            Va += m.second * n * (1. - S);
            if (airMass) *airMass += m.second * n * (1. - S) * (pa + patm);
        }
    }

    /// equations of p_w (space 1) or p_a (space 2) of a set of nodes
    std::vector<int64_t> Eqs(const std::map<int64_t, REAL> &nodes, int space) const {
        std::vector<int64_t> e;
        for (auto &n : nodes) e.push_back((space == 1 ? an->WaterEq() : an->AirEq()).at(n.first));
        return e;
    }
};

#endif
