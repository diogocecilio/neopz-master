// 3D version of ../SlopeMohrCoulomb: the same slope (70 x 40 m, height 10 m at 45 degrees) extruded
// 40 m along z, associative perfectly plastic Mohr-Coulomb, factor of safety by gravity increase
// and by shear strength reduction. The front and back faces (z = 0, z = Lz) are on rollers (u_z = 0),
// so the exact solution is the plane-strain one of the 2D project.
//
// Mesh: the 2D mesh of the original project at one level of uniform refinement (TriGMesh(1),
// elements of 5 m) extruded in layers of the same size:
//  - tetra (default): each triangular prism split into 3 tetrahedra (conforming diagonals);
//  - hexa: a structured quadrilateral version of the 2D mesh extruded into hexahedra (prepared,
//    same boundary ids).
//
// Usage: SlopeMohrCoulomb3D [hexa] [pv] [mesh] [nref=<n>] [ref=<n>] [nz=<n>] [lz=<m>] [p=<order>] [nu=<poisson>]
//   mesh   : only writes the geometric mesh (VTK) and the number of equations, no analysis
//   ref    : uniform refinements of the 2D base mesh (default 1, as in the 2D project)
//   nref   : adaptive refinement cycles of the plastic zone (default 0; tetrahedra become pyramids)
#include "SlopeAnalysis3D.h"
#include "Plasticity/TPZMatElastoPlastic_impl.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoCube.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoelbc.h"
#include "pzgeotetrahedra.h"
#include "tpzgeoelrefpattern.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using TMCVoigt = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMCPV = TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;

// the library instantiates the 3D material for TPZPlasticStepPV<TPZYCMohrCoulombPV> only
template class TPZMatElastoPlastic<TMCVoigt, TPZElastoPlasticMem>;

/// TPZMatElastoPlastic ignores the body force of SetBodyForce (it only uses a forcing function times
/// the bulk density); this adds m_force per unit volume, as TPZMatElastoPlastic2D does.
template <class TPlastic>
class TPZMatElastoPlasticGravity : public TPZMatElastoPlastic<TPlastic, TPZElastoPlasticMem> {
    using TBase = TPZMatElastoPlastic<TPlastic, TPZElastoPlasticMem>;

public:
    using TBase::TBase;
    TPZMaterial *NewMaterial() const override { return new TPZMatElastoPlasticGravity(*this); }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ek,
                    TPZFMatrix<REAL> &ef) override {
        TBase::Contribute(data, weight, ek, ef);
        AddBodyForce(data, weight, ef);
    }

    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ef) override {
        TBase::Contribute(data, weight, ef);
        AddBodyForce(data, weight, ef);
    }

private:
    void AddBodyForce(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ef) const {
        if (this->fForcingFunction) return; // the base class already applies it
        for (int in = 0; in < data.phi.Rows(); in++)
            for (int k = 0; k < 3; k++) ef(3 * in + k, 0) += weight * data.phi(in, 0) * this->m_force[k];
    }
};

struct Soil {
    REAL E = 20000., nu = 0.49, c = 10., phi = 30. * M_PI / 180., gamma = 20.;
};

TMCVoigt ModelVoigt(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    return TMCVoigt(TPZYCMohrCoulombPV2(s.phi, s.phi, s.c, ER), ER);
}

TMCPV ModelPV(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    TMCPV pv;
    pv.fYC.SetUp(s.phi, s.phi, s.c, ER);
    pv.fER = ER;
    return pv;
}

/// 2D mesh in the (x, y) plane: nodes and cells (3 nodes: triangles, 4 nodes: quadrilaterals, CCW)
struct Mesh2D {
    std::vector<std::array<REAL, 2>> x;
    std::vector<std::vector<int64_t>> cells;
};

/// Coarse triangles of TriGMesh (2D project): slope 70 x 40 m, height 10 m at 45 degrees
Mesh2D CoarseTriangles() {
    Mesh2D m;
    m.x = {{0, 0},   {10, 0},  {20, 0},  {30, 0},  {40, 0},  {50, 0},  {60, 0},  {70, 0},  {0, 10},
           {10, 10}, {20, 10}, {30, 10}, {40, 10}, {50, 10}, {60, 10}, {70, 10}, {0, 20},  {10, 20},
           {20, 20}, {30, 20}, {40, 20}, {50, 20}, {60, 20}, {70, 20}, {0, 30},  {10, 30}, {20, 30},
           {30, 30}, {40, 30}, {50, 30}, {60, 30}, {70, 30}, {0, 40},  {10, 40}, {20, 40}, {30, 40}};
    m.cells = {
        {0, 1, 8},    {1, 9, 8},    {1, 2, 9},    {2, 10, 9},   {2, 3, 10},   {3, 11, 10},  {3, 4, 11},
        {4, 12, 11},  {4, 5, 12},   {5, 13, 12},  {5, 6, 13},   {6, 14, 13},  {6, 7, 14},   {7, 15, 14},
        {8, 9, 16},   {9, 17, 16},  {9, 10, 17},  {10, 18, 17}, {10, 11, 18}, {11, 19, 18}, {11, 12, 19},
        {12, 20, 19}, {12, 13, 20}, {13, 21, 20}, {13, 14, 21}, {14, 22, 21}, {14, 15, 22}, {15, 23, 22},
        {16, 17, 24}, {17, 25, 24}, {17, 18, 25}, {18, 26, 25}, {18, 19, 26}, {19, 27, 26}, {19, 20, 27},
        {20, 28, 27}, {20, 21, 28}, {21, 29, 28}, {21, 22, 29}, {22, 30, 29}, {22, 23, 30}, {23, 31, 30},
        {24, 25, 32}, {25, 33, 32}, {25, 26, 33}, {26, 34, 33}, {26, 27, 34}, {27, 35, 34}, {27, 28, 35}};
    return m;
}

/// Same domain with quadrilaterals: 10 m squares below y = 30 and a mapped band between y = 30
/// and the crest (x from 0 to the slope face, 4 elements), conforming at y = 30
Mesh2D CoarseQuads() {
    Mesh2D m;
    for (int j = 0; j < 4; j++)
        for (int i = 0; i < 8; i++) m.x.push_back({10. * i, 10. * j}); // nodes 0..31
    for (int i = 0; i < 5; i++) m.x.push_back({7.5 * i, 40.});        // nodes 32..36 (crest, x <= 30)
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 7; i++) m.cells.push_back({8 * j + i, 8 * j + i + 1, 8 * (j + 1) + i + 1, 8 * (j + 1) + i});
    for (int i = 0; i < 4; i++) m.cells.push_back({24 + i, 25 + i, 33 + i, 32 + i});
    return m;
}

/// Uniform refinement by edge midpoints (and centre for quadrilaterals), as the uniform refinement
/// patterns of NeoPZ for linear triangles and quadrilaterals
Mesh2D Refine2D(const Mesh2D &m) {
    Mesh2D r;
    r.x = m.x;
    std::map<std::pair<int64_t, int64_t>, int64_t> mid;
    auto node = [&](std::array<REAL, 2> p) { r.x.push_back(p); return int64_t(r.x.size() - 1); };
    auto edge = [&](int64_t a, int64_t b) {
        const auto key = std::minmax(a, b);
        auto it = mid.find(key);
        if (it != mid.end()) return it->second;
        const int64_t n = node({0.5 * (m.x[a][0] + m.x[b][0]), 0.5 * (m.x[a][1] + m.x[b][1])});
        mid[key] = n;
        return n;
    };
    for (auto &c : m.cells) {
        if (c.size() == 3) {
            const int64_t ab = edge(c[0], c[1]), bc = edge(c[1], c[2]), ca = edge(c[2], c[0]);
            r.cells.push_back({c[0], ab, ca});
            r.cells.push_back({ab, c[1], bc});
            r.cells.push_back({ca, bc, c[2]});
            r.cells.push_back({ab, bc, ca});
        } else {
            const int64_t ab = edge(c[0], c[1]), bc = edge(c[1], c[2]), cd = edge(c[2], c[3]), da = edge(c[3], c[0]);
            std::array<REAL, 2> xc = {0., 0.};
            for (int k = 0; k < 4; k++) { xc[0] += 0.25 * m.x[c[k]][0]; xc[1] += 0.25 * m.x[c[k]][1]; }
            const int64_t mc = node(xc);
            r.cells.push_back({c[0], ab, mc, da});
            r.cells.push_back({ab, c[1], bc, mc});
            r.cells.push_back({mc, bc, c[2], cd});
            r.cells.push_back({da, mc, cd, c[3]});
        }
    }
    return r;
}

/// Extrusion of the 2D mesh along z (nz layers, thickness lz). Triangles: each prism is split into
/// 3 tetrahedra with the face diagonals chosen by the global node numbers (conforming between
/// neighbours); quadrilaterals: hexahedra. Boundary faces (ids as in the 2D project):
/// -1 base y = 0, -2 right x = 70, -5 left x = 0, -7 front z = 0, -8 back z = lz; free otherwise.
TPZGeoMesh *ExtrudeGMesh(const Mesh2D &m, REAL lz, int nz) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(3);
    const int64_t n2 = m.x.size();
    gmesh->NodeVec().Resize(n2 * (nz + 1));
    for (int l = 0; l <= nz; l++)
        for (int64_t i = 0; i < n2; i++) {
            TPZManVector<REAL, 3> x = {m.x[i][0], m.x[i][1], lz * l / nz};
            gmesh->NodeVec()[l * n2 + i] = TPZGeoNode(l * n2 + i, x, *gmesh);
        }
    auto coord = [&](int64_t n, int k) { return gmesh->NodeVec()[n].Coord(k); };
    auto tetra = [&](TPZManVector<int64_t, 4> nodes) {
        REAL a[3][3]; // positive volume: (x1-x0) . ((x2-x0) x (x3-x0)) > 0
        for (int i = 0; i < 3; i++)
            for (int k = 0; k < 3; k++) a[i][k] = coord(nodes[i + 1], k) - coord(nodes[0], k);
        const REAL vol = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                         a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                         a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        if (vol < 0.) std::swap(nodes[1], nodes[2]);
        new TPZGeoElRefPattern<pzgeom::TPZGeoTetrahedra>(nodes, 1, *gmesh);
    };
    for (int l = 0; l < nz; l++) {
        const int64_t bot = l * n2, top = (l + 1) * n2;
        for (auto &c : m.cells) {
            if (c.size() == 4) {
                TPZManVector<int64_t, 8> nodes(8);
                for (int k = 0; k < 4; k++) { nodes[k] = bot + c[k]; nodes[k + 4] = top + c[k]; }
                new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(nodes, 1, *gmesh);
            } else {
                std::vector<int64_t> v(c); // v0 < v1 < v2: quad face (vi, vj), i < j, cut along vi' - vj
                std::sort(v.begin(), v.end());
                tetra({bot + v[0], bot + v[1], bot + v[2], top + v[0]});
                tetra({bot + v[1], bot + v[2], top + v[0], top + v[1]});
                tetra({bot + v[2], top + v[0], top + v[1], top + v[2]});
            }
        }
    }
    gmesh->BuildConnectivity();
    const REAL tol = 1.e-6, xmax = 70.;
    const int64_t nvol = gmesh->NElements();
    for (int64_t el = 0; el < nvol; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        for (int s = 0; s < gel->NSides(); s++) {
            TPZGeoElSide side(gel, s);
            if (side.Dimension() != 2 || side.Neighbour() != side) continue; // interior face
            TPZManVector<REAL, 3> xc(3);
            side.CenterX(xc);
            int id = 0;
            if (std::fabs(xc[1]) < tol) id = -1;
            else if (std::fabs(xc[0] - xmax) < tol) id = -2;
            else if (std::fabs(xc[0]) < tol) id = -5;
            else if (std::fabs(xc[2]) < tol) id = -7;
            else if (std::fabs(xc[2] - lz) < tol) id = -8;
            if (id) TPZGeoElBC bc(side, id);
        }
    }
    return gmesh;
}

/// H1 mesh with memory; base fixed, lateral sides and front/back on rollers
/// (BC type 3: val2 = constrained directions)
template <class TPlastic>
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int porder, const TPlastic &model, const Soil &s) {
    auto *mat = new TPZMatElastoPlasticGravity<TPlastic>(1);
    TPlastic m(model);
    mat->SetPlasticityModel(m);
    mat->SetBodyForce({0., -s.gamma, 0.});
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetDimModel(3);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(3, 3, 0.);
    const std::vector<std::pair<int, TPZManVector<STATE, 3>>> bcs = {
        {-1, {1., 1., 1.}}, {-2, {1., 0., 0.}}, {-5, {1., 0., 0.}}, {-7, {0., 0., 1.}}, {-8, {0., 0., 1.}}};
    for (auto &bc : bcs) cmesh->InsertMaterialObject(mat->CreateBC(mat, bc.first, 3, val1, bc.second));
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

struct Options {
    bool hexa = false, pv = false, meshonly = false;
    int nref = 0, ref = 1, nz = -1, porder = 2;
    REAL lz = 40.;
};

TPZGeoMesh *GMesh(const Options &o) {
    Mesh2D m = o.hexa ? CoarseQuads() : CoarseTriangles();
    for (int k = 0; k < o.ref; k++) m = Refine2D(m);
    const REAL h = 10. / (1 << o.ref); // in-plane element size
    const int nz = o.nz > 0 ? o.nz : std::max(1, int(std::lround(o.lz / h)));
    return ExtrudeGMesh(m, o.lz, nz);
}

template <class TPlastic>
void Run(const TPlastic &model, const Soil &s, const Options &o, const std::string &tag) {
    TPZGeoMesh *gmesh = GMesh(o);
    {
        std::ofstream vtk(tag + "_gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtk, true);
    }
    TPZCompMesh *cmesh = CreateCMesh(gmesh, o.porder, model, s);
    std::cout << tag << ": " << gmesh->NElements() << " geometric elements, " << cmesh->NEquations() << " equations\n";
    if (o.meshonly) {
        delete cmesh;
        delete gmesh;
        return;
    }
    if (!o.hexa && o.nref > 0)
        std::cout << "warning: refining tetrahedra creates pyramids (uniform pattern: 4 tetrahedra + 2 pyramids)\n";
    SlopeAnalysis3D<TPlastic> slope(cmesh);
    std::vector<std::vector<REAL>> table;
    for (int k = 0;; k++) {
        const int64_t neq = cmesh->NEquations();
        const REAL fsGI = slope.GravityIncrease();
        slope.PostProcess(tag + "_GI_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1); // failure mechanism at collapse
        const REAL fsSRM = slope.StrengthReduction();
        slope.PostProcess(tag + "_SRM_ref" + std::to_string(k) + ".vtk");
        slope.MarkPlasticZone(0.1);
        table.push_back({REAL(k), REAL(neq), fsGI, fsSRM});
        if (k == o.nref) break;
        slope.Refine();
    }
    std::ofstream vtk(tag + "_mesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtk, true);
    std::cout << "\n" << tag << ": refinement  equations  FS(gravity increase)  FS(strength reduction)\n";
    for (auto &r : table) std::cout << "  " << r[0] << "  " << r[1] << "  " << r[2] << "  " << r[3] << "\n";
    delete cmesh;
    delete gmesh;
}

int main(int argc, char *argv[]) {
    Options o;
    Soil s;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "hexa")) o.hexa = true;
        else if (!strcmp(argv[i], "pv")) o.pv = true;
        else if (!strcmp(argv[i], "mesh")) o.meshonly = true;
        else if (!strncmp(argv[i], "nref=", 5)) o.nref = atoi(argv[i] + 5);
        else if (!strncmp(argv[i], "ref=", 4)) o.ref = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "nz=", 3)) o.nz = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "lz=", 3)) o.lz = atof(argv[i] + 3);
        else if (!strncmp(argv[i], "p=", 2)) o.porder = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "nu=", 3)) s.nu = atof(argv[i] + 3);
    }
    const std::string tag = std::string("slope3d_") + (o.hexa ? "hexa" : "tetra") + (o.pv ? "_pv" : "_rhw");
    if (o.pv) Run(ModelPV(s), s, o, tag);
    else Run(ModelVoigt(s), s, o, tag);
    return 0;
}
