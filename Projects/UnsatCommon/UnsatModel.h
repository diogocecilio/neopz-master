// Unsaturated soil (u - p_w - p_a) with the closed-form Mohr-Coulomb of Lira Cecilio (J. Eng. Math. 157:10, 2026):
// parameters, material, structured meshes and the multiphysics mesh shared by the test projects.
// Formulation: docs/biot_nao_saturado/formulacao/formulacao.tex.
#ifndef UNSATMODEL_H
#define UNSATMODEL_H

#include "Plasticity/TPZMatUnsaturatedUPwPa.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzgeoel.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <vector>

/// Associative Mohr-Coulomb (psi = phi) in rotated Haigh-Westergaard space with consistent tangent
using TMC = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TUnsat = TPZMatUnsaturatedUPwPa<TMC, TPZUnsaturatedMem>;

/// Boundary ids of the structured meshes
enum EBnd { EBase = -1, ERight = -2, ETop = -3, ELeft = -4 };

/// Parameters of a soil (units kN, m, s, kPa); defaults: chi = S_w (notes), van Genuchten - Mualem
struct UnsatParams {
    REAL E = 2.e4, nu = 0.3, phi = 30. * M_PI / 180., c = 10.;
    REAL gamma = 20.;            ///< constant unit weight (kN/m3)
    bool mixture = false;        ///< unit weight of the mixture (1-n) gamma_s + n S_w gamma_w
    REAL gammas = 26.;           ///< unit weight of the grains (mixture)
    TUnsat::EBishop bishop = TUnsat::EChiSaturation;
    REAL chi = 1.;               ///< constant chi (bishop = EChiConstant)
    REAL n = 0.45, Kw = 2.2e6, gw = 9.81, ga = 1.2 * 9.81e-3, muRatio = 1.e-3 / 1.8e-5, patm = 101.325;
    REAL ks = 1.e-5;             ///< saturated hydraulic conductivity (m/s)
    REAL alpha = 1., Ks = 0.;    ///< Biot coefficient and bulk modulus of the grains (<= 0: incompressible)
    bool coupled = true;         ///< hydro-mechanical coupling of the mass balances
    bool lumped = true;          ///< lumped storage
    std::array<REAL, 3> g = {0., -1., 0.}; ///< unit vector of gravity
    TPZUnsaturatedHydraulics hydro;
    UnsatParams() {
        hydro.fRetention = TPZUnsaturatedHydraulics::EVanGenuchten;
        hydro.fKrw = TPZUnsaturatedHydraulics::EMualem;
        hydro.fKra = TPZUnsaturatedHydraulics::EMualemAir;
    }
};

inline TMC MohrCoulomb(const UnsatParams &p) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(p.E, p.nu);
    return TMC(TPZYCMohrCoulombPV2(p.phi, p.phi, p.c, ER), ER);
}

inline TUnsat *UnsatMaterial(int id, const UnsatParams &p) {
    auto *mat = new TUnsat(id);
    mat->SetPlasticModel(MohrCoulomb(p));
    mat->SetHydraulics(p.hydro);
    mat->SetPorosity(p.n);
    mat->SetBiot(p.alpha, p.Ks);
    mat->SetWaterBulkModulus(p.Kw);
    mat->SetConductivity(p.ks, p.gw, p.muRatio);
    mat->SetAir(p.ga, p.patm);
    if (p.mixture) mat->SetUnitWeightMixture(p.gammas);
    else mat->SetUnitWeight(p.gamma);
    mat->SetGravityDirection(TPZManVector<REAL, 3>({p.g[0], p.g[1], p.g[2]}));
    mat->SetBishop(p.bishop, p.chi);
    mat->SetLumpedStorage(p.lumped);
    mat->SetHydroMechanicalCoupling(p.coupled);
    mat->SetIntegrationOrder(2); // 3-point rule
    return mat;
}

/// Graded coordinates: n intervals on [a, b], ratio of the last to the first interval = r
inline std::vector<REAL> Graded(REAL a, REAL b, int n, REAL r = 1.) {
    std::vector<REAL> x(n + 1);
    const REAL q = n > 1 ? std::pow(r, 1. / (n - 1)) : 1.;
    REAL sum = 0., h = 1.;
    for (int i = 0; i < n; ++i) sum += std::pow(q, i);
    x[0] = a;
    for (int i = 0; i < n; ++i) {
        x[i + 1] = x[i] + (b - a) * std::pow(q, i) / sum;
        (void)h;
    }
    x[n] = b;
    return x;
}

/// Joins coordinate lists (the first point of each list after the first is dropped)
inline std::vector<REAL> Join(std::vector<std::vector<REAL>> parts) {
    std::vector<REAL> x;
    for (auto &p : parts)
        for (size_t i = (x.empty() ? 0 : 1); i < p.size(); ++i) x.push_back(p[i]);
    return x;
}

/**
 * Structured triangular mesh: vertical lines at xs; on each line, nodes at y = f * ytop(x) for the fractions fr
 * (0 = base, 1 = top). Each quadrilateral is split in two triangles. Domain id by the centroid (default 1);
 * boundary ids: base -1, right -2, top -3, left -4.
 */
inline TPZGeoMesh *StructuredMesh(const std::vector<REAL> &xs, const std::function<REAL(REAL)> &ytop,
                                  const std::vector<REAL> &fr,
                                  const std::function<int(REAL, REAL)> &idOf = nullptr) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    const int nx = int(xs.size()) - 1, ny = int(fr.size()) - 1;
    gmesh->NodeVec().Resize((nx + 1) * (ny + 1));
    auto node = [&](int i, int j) { return int64_t(i) * (ny + 1) + j; };
    for (int i = 0; i <= nx; ++i)
        for (int j = 0; j <= ny; ++j) {
            TPZManVector<REAL, 3> x = {xs[i], fr[j] * ytop(xs[i]), 0.};
            gmesh->NodeVec()[node(i, j)] = TPZGeoNode(node(i, j), x, *gmesh);
        }
    auto coord = [&](int64_t k, int d) { return gmesh->NodeVec()[k].Coord(d); };
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            const int64_t a = node(i, j), b = node(i + 1, j), c = node(i + 1, j + 1), d = node(i, j + 1);
            for (auto t : {std::array<int64_t, 3>{a, b, c}, std::array<int64_t, 3>{a, c, d}}) {
                const REAL xc = (coord(t[0], 0) + coord(t[1], 0) + coord(t[2], 0)) / 3.;
                const REAL yc = (coord(t[0], 1) + coord(t[1], 1) + coord(t[2], 1)) / 3.;
                TPZManVector<int64_t, 3> nn = {t[0], t[1], t[2]};
                new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nn, idOf ? idOf(xc, yc) : 1, *gmesh);
            }
        }
    auto line = [&](int64_t a, int64_t b, int id) {
        TPZManVector<int64_t, 2> nn = {a, b};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nn, id, *gmesh);
    };
    for (int i = 0; i < nx; ++i) line(node(i, 0), node(i + 1, 0), EBase);
    for (int j = 0; j < ny; ++j) line(node(nx, j), node(nx, j + 1), ERight);
    for (int i = nx; i > 0; --i) line(node(i, ny), node(i - 1, ny), ETop);
    for (int j = ny; j > 0; --j) line(node(0, j), node(0, j - 1), ELeft);
    gmesh->BuildConnectivity();
    return gmesh;
}

/// Rectangle [0, W] x [0, H] with nx x ny cells (fractions in y from fr if given)
inline TPZGeoMesh *RectangleMesh(REAL W, REAL H, int nx, int ny, std::vector<REAL> fr = {},
                                 const std::function<int(REAL, REAL)> &idOf = nullptr) {
    if (fr.empty()) fr = Graded(0., 1., ny);
    return StructuredMesh(Graded(0., W, nx), [H](REAL) { return H; }, fr, idOf);
}

/**
 * Slope of SlopeMohrCoulomb (Lira Cecilio 2026 implementation, docs/slope_mohr_coulomb): domain 70 m wide, crest at
 * y = 40 m for x <= 30 m, face at 45 degrees down to the toe at y = 30 m for x >= 40 m. Graded towards the face.
 */
inline TPZGeoMesh *SlopeMesh(int level, REAL surface = 0.6) {
    const int k = std::max(level, 1);
    const std::vector<REAL> xs = Join({Graded(0., 20., 3 * k, 0.4), Graded(20., 30., 4 * k), Graded(30., 40., 5 * k),
                                       Graded(40., 52., 5 * k), Graded(52., 70., 3 * k, 2.5)});
    // rows: graded towards the surface (surface = last / first interval of the upper half)
    const std::vector<REAL> fr = Join({Graded(0., 0.5, 3 * k, 0.5), Graded(0.5, 1., 7 * k, surface)});
    auto ytop = [](REAL x) { return x <= 30. ? 40. : (x >= 40. ? 30. : 40. - (x - 30.)); };
    return StructuredMesh(xs, ytop, fr);
}

/**
 * Same slope with the elements concentrated in the zone of the failure mechanisms: nf uniform columns between
 * x = 26 and 46 m and nf/2 rows in the upper 40% of each column, graded towards the surface.
 */
inline TPZGeoMesh *SlopeMeshFace(int nf) {
    const std::vector<REAL> xs = Join({Graded(0., 26., std::max(nf / 8, 3), 0.25), Graded(26., 46., nf),
                                       Graded(46., 70., std::max(nf / 6, 3), 4.)});
    const std::vector<REAL> fr = Join({Graded(0., 0.6, std::max(nf / 8, 3), 0.3), Graded(0.6, 1., nf / 2, 0.25)});
    auto ytop = [](REAL x) { return x <= 30. ? 40. : (x >= 40. ? 30. : 40. - (x - 30.)); };
    return StructuredMesh(xs, ytop, fr);
}

/// H1 atomic space with null materials (domain ids and all the boundary ids of the mesh)
inline TPZCompMesh *AtomicMesh(TPZGeoMesh *gmesh, int nstate, int order) {
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    std::set<int> dom, bnd;
    for (int64_t i = 0; i < gmesh->NElements(); ++i) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel) continue;
        (gel->Dimension() == 2 ? dom : bnd).insert(gel->MaterialId());
    }
    TPZNullMaterial<STATE> *mat = nullptr;
    for (int id : dom) {
        mat = new TPZNullMaterial<STATE>(id, 2, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    TPZFMatrix<STATE> val1(nstate, nstate, 0.);
    TPZManVector<STATE, 2> val2(nstate, 0.);
    for (int id : bnd) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

/// Traction on a boundary (kPa), assembled by the material (ENeumannU)
struct Traction {
    int id;
    REAL tx, ty;
};

/// Multiphysics mesh u (P2) - p_w (P1) - p_a (P1): material per domain id (params[id], or params.begin())
inline TPZMultiphysicsCompMesh *UnsatMesh(TPZGeoMesh *gmesh, const std::map<int, UnsatParams> &params,
                                          const std::vector<Traction> &tractions = {}) {
    TPZManVector<TPZCompMesh *, 3> meshvec = {AtomicMesh(gmesh, 2, 2), AtomicMesh(gmesh, 1, 1), AtomicMesh(gmesh, 1, 1)};
    auto *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    std::set<int> dom, bnd;
    for (int64_t i = 0; i < gmesh->NElements(); ++i) {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel) continue;
        (gel->Dimension() == 2 ? dom : bnd).insert(gel->MaterialId());
    }
    TUnsat *mat = nullptr;
    for (int id : dom) {
        auto it = params.find(id);
        mat = UnsatMaterial(id, it == params.end() ? params.begin()->second : it->second);
        mphys->InsertMaterialObject(mat);
    }
    TPZFMatrix<STATE> val1(2, 2, 0.);
    for (int id : bnd) {
        TPZManVector<STATE, 2> val2(2, 0.);
        int type = TUnsat::EDirichletU; // no contribution: the Dirichlet conditions are imposed by the analysis
        for (auto &t : tractions)
            if (t.id == id) type = TUnsat::ENeumannU, val2[0] = t.tx, val2[1] = t.ty;
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, type, val1, val2));
    }
    TPZManVector<int, 3> active(3, 1);
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, dom, bnd);
    return mphys;
}

#endif
