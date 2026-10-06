/**
 * @file MCCPaperTools.h
 * @brief Utilities shared by the examples of the article "Return mapping for Modified Cam-Clay
 * plasticity in rotated Haigh-Westergaard space with consistent tangent operator and coupled u-p
 * consolidation": structured meshes (Q8-Q4, Hex20-Hex8: rectangle, box, block, unit cube, slab and quarter
 * cylinder), CSV files of the geometric mesh for plotting (WriteMeshCSV), atomic and multiphysics
 * computational meshes, post-processing at the integration points, VTK file series of every converged state
 * for ParaView (TVTKSeries), material point drivers and the closed-form solutions of Appendix B.
 *
 * The functions mirror the routines of the Python transcription of the Wolfram Language packages
 * (fe_user.py, camclay_hw.py): SubdivideQuadMesh, BoxMesh3D, QuarterCylinderMesh, LocatePoint,
 * GaussPointWeights, TriaxialPointCC, TriaxialDrainedClosedCC and the Terzaghi series.
 */

#ifndef MCCPAPERTOOLS_H
#define MCCPAPERTOOLS_H

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticquad.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcube.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzintel.h"
#include "pzquad.h"
#include "TPZTensor.h"
#include "TPZPlasticStepModifiedCamClay.h"
#include "TPZMatPoroElastoPlasticUP.h"
#include "TPZPoroElastoPlasticUPAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"

#include <vector>
#include <map>
#include <set>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <numeric>
#include <memory>
#include <string>
#include <filesystem>

/**
 * @ingroup mccpaper
 * @brief Utilities of the examples of the Modified Cam-Clay u-p article
 */
namespace mcc {

/** @brief The plastic model of all the examples */
typedef TPZPlasticStepModifiedCamClay TPlastic;
/** @brief The u-p material of all the finite element examples */
typedef TPZMatPoroElastoPlasticUP<TPlastic, TPZElastoPlasticMem> TPoroMaterial;
/** @brief The incremental Newton driver */
typedef TPZPoroElastoPlasticUPAnalysis TAnalysis;

/** @name Stress invariants (soil mechanics signs: compression positive) */
/** @{ */
/** @brief Mean effective stress \f$p'=-I_1/3\f$ */
inline REAL MeanEffectiveStress(const TPZTensor<REAL> &sig) { return -sig.I1() / 3.; }
/** @brief Deviatoric stress \f$q=\sqrt{3J_2}\f$ */
inline REAL DeviatoricStress(const TPZTensor<REAL> &sig) { return std::sqrt(3. * std::max(REAL(0.), sig.J2())); }
/** @brief Isotropic stress tensor \f$pI\f$ */
inline TPZTensor<REAL> IsotropicTensor(REAL p) {
    TPZTensor<REAL> t;
    t.XX() = t.YY() = t.ZZ() = p;
    return t;
}
/** @brief Specific volume on the normal compression line: \f$v_0=v_\lambda-\lambda\ln p'_{c0}+\kappa\ln(p'_{c0}/p'_0)\f$ */
inline REAL SpecificVolumeNCL(REAL vlambda, REAL lambda, REAL kappa, REAL pc0, REAL p0) {
    return vlambda - lambda * std::log(pc0) + kappa * std::log(pc0 / p0);
}
/** @} */

/** @name Geometric meshes */
/** @{ */

/**
 * @brief Marker function of the boundary edges of a 2D structured mesh
 * @param side 0 bottom, 1 right, 2 top, 3 left
 * @param xmid midpoint of the boundary edge
 * @return material ids of the boundary elements created on the edge (several ids create coincident elements)
 */
typedef std::function<std::vector<int>(int side, const TPZVec<REAL> &xmid)> TMarker2D;

/**
 * @brief Structured mesh of nx x ny quadrilaterals on the rectangle [x0,x1] x [y0,y1] (SubdivideQuadMesh);
 * nodes numbered row by row, boundary line elements created with the ids returned by marker
 */
inline TPZGeoMesh *CreateRectangleMesh(REAL x0, REAL y0, REAL x1, REAL y1, int nx, int ny, int matid,
                                       const TMarker2D &marker) {
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    auto idx = [nx](int i, int j) { return int64_t(j) * (nx + 1) + i; };
    gmesh->NodeVec().Resize((nx + 1) * (ny + 1));
    TPZManVector<REAL, 3> co(3, 0.);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            co[0] = x0 + (x1 - x0) * i / nx;
            co[1] = y0 + (y1 - y0) * j / ny;
            gmesh->NodeVec()[idx(i, j)].Initialize(co, *gmesh);
        }
    }
    TPZManVector<int64_t, 4> quad(4);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            quad[0] = idx(i, j);
            quad[1] = idx(i + 1, j);
            quad[2] = idx(i + 1, j + 1);
            quad[3] = idx(i, j + 1);
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(quad, matid, *gmesh);
        }
    }
    auto addline = [&](int side, int64_t a, int64_t b) {
        TPZManVector<REAL, 3> xa(3), xb(3), xm(3);
        gmesh->NodeVec()[a].GetCoordinates(xa);
        gmesh->NodeVec()[b].GetCoordinates(xb);
        for (int k = 0; k < 3; ++k) xm[k] = 0.5 * (xa[k] + xb[k]);
        TPZManVector<int64_t, 2> line(2);
        line[0] = a;
        line[1] = b;
        for (int id : marker(side, xm)) new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(line, id, *gmesh);
    };
    for (int i = 0; i < nx; ++i) addline(0, idx(i, 0), idx(i + 1, 0));
    for (int j = 0; j < ny; ++j) addline(1, idx(nx, j), idx(nx, j + 1));
    for (int i = 0; i < nx; ++i) addline(2, idx(i + 1, ny), idx(i, ny));
    for (int j = 0; j < ny; ++j) addline(3, idx(0, j + 1), idx(0, j));
    gmesh->BuildConnectivity();
    return gmesh;
}

/**
 * @brief Marker function of the boundary faces of a 3D mesh
 * @param X coordinates of the four vertices of the face (4 x 3)
 * @return material ids of the boundary elements created on the face
 */
typedef std::function<std::vector<int>(const std::array<std::array<REAL, 3>, 4> &X)> TMarker3D;

/**
 * @brief Helper of the marker functions: true if the four vertices of a face lie on the plane
 * \f$x_{axis}=value\f$
 * @param X coordinates of the four vertices of the face (argument of TMarker3D)
 * @param axis 0 (x), 1 (y) or 2 (z)
 * @param value coordinate of the plane
 * @param tol absolute tolerance on the coordinate
 */
inline bool FaceOnPlane(const std::array<std::array<REAL, 3>, 4> &X, int axis, REAL value, REAL tol = 1.e-9) {
    for (auto &p : X)
        if (std::fabs(p[axis] - value) > tol) return false;
    return true;
}

namespace internal {
/**
 * @brief Builds a hexahedral mesh from the vertices and the cells (8 vertices, PZ/Python ordering);
 * order 1 creates trilinear elements (TPZGeoCube), order 2 creates 20-node serendipity elements
 * (TPZQuadraticCube) whose mid-edge nodes are placed by midnode (default: midpoint).
 * Boundary faces (faces of a single cell) receive the ids returned by marker.
 */
inline TPZGeoMesh *HexMesh(const std::vector<std::array<REAL, 3>> &vertices,
                           const std::vector<std::array<int64_t, 8>> &cells, int matid, const TMarker3D &marker,
                           int geoorder,
                           const std::function<void(const std::vector<int> &markers, std::array<REAL, 3> &x)> &movenode =
                               nullptr) {
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    std::vector<std::array<REAL, 3>> coords(vertices);
    std::map<std::pair<int64_t, int64_t>, int64_t> edgenode;
    auto midnode = [&](int64_t a, int64_t b) {
        std::pair<int64_t, int64_t> key(std::min(a, b), std::max(a, b));
        auto it = edgenode.find(key);
        if (it != edgenode.end()) return it->second;
        std::array<REAL, 3> x;
        for (int k = 0; k < 3; ++k) x[k] = 0.5 * (coords[a][k] + coords[b][k]);
        coords.push_back(x);
        edgenode[key] = int64_t(coords.size()) - 1;
        return int64_t(coords.size()) - 1;
    };
    // element node lists (vertices, then the 12 edges in the order of the sides 8-19 of TPZCube)
    std::vector<std::vector<int64_t>> elnodes;
    for (auto &c : cells) {
        std::vector<int64_t> nodes(c.begin(), c.end());
        if (geoorder == 2) {
            for (int side = 8; side < 20; ++side) {
                const int a = pztopology::TPZCube::SideNodeLocId(side, 0);
                const int b = pztopology::TPZCube::SideNodeLocId(side, 1);
                nodes.push_back(midnode(c[a], c[b]));
            }
        }
        elnodes.push_back(nodes);
    }
    // boundary faces: faces (sides 20-25) that belong to a single cell
    std::map<std::array<int64_t, 4>, std::vector<std::pair<int, int>>> faces;
    for (int ic = 0; ic < (int)cells.size(); ++ic) {
        for (int side = 20; side < 26; ++side) {
            std::array<int64_t, 4> key;
            for (int k = 0; k < 4; ++k) key[k] = cells[ic][pztopology::TPZCube::SideNodeLocId(side, k)];
            std::sort(key.begin(), key.end());
            faces[key].push_back({ic, side});
        }
    }
    struct TFace {
        std::vector<int64_t> nodes;
        std::vector<int> ids;
    };
    std::vector<TFace> bcfaces;
    for (auto &f : faces) {
        if (f.second.size() != 1) continue;
        const int ic = f.second[0].first, side = f.second[0].second;
        std::array<std::array<REAL, 3>, 4> X;
        std::vector<int64_t> nodes(4);
        for (int k = 0; k < 4; ++k) {
            nodes[k] = cells[ic][pztopology::TPZCube::SideNodeLocId(side, k)];
            X[k] = vertices[nodes[k]];
        }
        TFace face;
        face.ids = marker(X);
        if (geoorder == 2) {
            for (int k = 0; k < 4; ++k) face.nodes.push_back(nodes[k]);
            for (int k = 0; k < 4; ++k) face.nodes.push_back(midnode(nodes[k], nodes[(k + 1) % 4]));
        } else {
            face.nodes = nodes;
        }
        bcfaces.push_back(face);
    }
    // optional relocation of the nodes of the boundary faces (curved boundaries)
    if (movenode) {
        std::map<int64_t, std::vector<int>> nodemarkers;
        for (auto &f : bcfaces)
            for (auto n : f.nodes)
                for (int id : f.ids) nodemarkers[n].push_back(id);
        for (auto &nm : nodemarkers) movenode(nm.second, coords[nm.first]);
    }
    gmesh->NodeVec().Resize(coords.size());
    for (size_t i = 0; i < coords.size(); ++i) {
        TPZManVector<REAL, 3> co(3);
        for (int k = 0; k < 3; ++k) co[k] = coords[i][k];
        gmesh->NodeVec()[i].Initialize(co, *gmesh);
    }
    for (auto &nodes : elnodes) {
        TPZManVector<int64_t, 20> topo(nodes.size());
        for (size_t k = 0; k < nodes.size(); ++k) topo[k] = nodes[k];
        if (geoorder == 2) new TPZGeoElRefPattern<pzgeom::TPZQuadraticCube>(topo, matid, *gmesh);
        else new TPZGeoElRefPattern<pzgeom::TPZGeoCube>(topo, matid, *gmesh);
    }
    for (auto &f : bcfaces) {
        TPZManVector<int64_t, 8> topo(f.nodes.size());
        for (size_t k = 0; k < f.nodes.size(); ++k) topo[k] = f.nodes[k];
        for (int id : f.ids) {
            if (geoorder == 2) new TPZGeoElRefPattern<pzgeom::TPZQuadraticQuad>(topo, id, *gmesh);
            else new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(topo, id, *gmesh);
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}
} // namespace internal

/**
 * @brief Box [0,Lx] x [0,Ly] x [0,Lz] divided into nx x ny x nz trilinear hexahedra (BoxMesh3D)
 */
inline TPZGeoMesh *CreateBoxMesh(REAL Lx, REAL Ly, REAL Lz, int nx, int ny, int nz, int matid,
                                 const TMarker3D &marker) {
    auto idx = [nx, ny](int i, int j, int k) { return int64_t(k) * (ny + 1) * (nx + 1) + int64_t(j) * (nx + 1) + i; };
    std::vector<std::array<REAL, 3>> vertices;
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) vertices.push_back({Lx * i / nx, Ly * j / ny, Lz * k / nz});
    std::vector<std::array<int64_t, 8>> cells;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                cells.push_back({idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                                 idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)});
    return internal::HexMesh(vertices, cells, matid, marker, 1);
}

/**
 * @brief Quarter of a cylinder of radius R and height H (QuarterCylinderMesh): central square
 * [0, a R]^2 with nc x nc cells and two outer blocks with nc cells along the arc and nr radially,
 * extruded in nz layers. The elements are 20-node quadratic hexahedra; the nodes of the lateral faces
 * (id lateralid) are moved radially to r = R, as in the isoparametric Python model.
 */
inline TPZGeoMesh *CreateQuarterCylinderMesh(REAL R, REAL H, int nc, int nr, int nz, int matid,
                                             const TMarker3D &marker, int lateralid, REAL aratio = 0.5) {
    const REAL a = aratio * R;
    std::map<std::pair<int64_t, int64_t>, int64_t> keys;
    std::vector<std::array<REAL, 2>> pts;
    auto node = [&](REAL x, REAL y) {
        std::pair<int64_t, int64_t> key(std::llround(x / R * 1e9), std::llround(y / R * 1e9));
        auto it = keys.find(key);
        if (it != keys.end()) return it->second;
        pts.push_back({x, y});
        keys[key] = int64_t(pts.size()) - 1;
        return int64_t(pts.size()) - 1;
    };
    std::vector<std::array<int64_t, 4>> quads;
    std::vector<std::vector<int64_t>> S(nc + 1, std::vector<int64_t>(nc + 1));
    for (int i = 0; i <= nc; ++i)
        for (int j = 0; j <= nc; ++j) S[i][j] = node(a * i / nc, a * j / nc);
    for (int i = 0; i < nc; ++i)
        for (int j = 0; j < nc; ++j) quads.push_back({S[i][j], S[i + 1][j], S[i + 1][j + 1], S[i][j + 1]});
    for (int block = 1; block <= 2; ++block) {
        std::vector<std::vector<int64_t>> Gb(nc + 1, std::vector<int64_t>(nr + 1));
        for (int j = 0; j <= nc; ++j) {
            REAL ix, iy, th;
            if (block == 1) {
                ix = a;
                iy = a * j / nc;
                th = M_PI / 4. * j / nc;
            } else {
                ix = a * j / nc;
                iy = a;
                th = M_PI / 2. - M_PI / 4. * j / nc;
            }
            const REAL ox = R * std::cos(th), oy = R * std::sin(th);
            for (int k = 0; k <= nr; ++k)
                Gb[j][k] = node(ix + (REAL(k) / nr) * (ox - ix), iy + (REAL(k) / nr) * (oy - iy));
        }
        for (int j = 0; j < nc; ++j)
            for (int k = 0; k < nr; ++k) quads.push_back({Gb[j][k], Gb[j][k + 1], Gb[j + 1][k + 1], Gb[j + 1][k]});
    }
    // counter-clockwise orientation
    for (auto &q : quads) {
        REAL area = 0.;
        for (int k = 0; k < 4; ++k) {
            const auto &p = pts[q[k]], &n = pts[q[(k + 1) % 4]];
            area += p[0] * n[1] - n[0] * p[1];
        }
        if (area < 0.) q = {q[0], q[3], q[2], q[1]};
    }
    const int64_t n2 = pts.size();
    std::vector<std::array<REAL, 3>> vertices;
    for (int k = 0; k <= nz; ++k)
        for (auto &p : pts) vertices.push_back({p[0], p[1], H * k / nz});
    std::vector<std::array<int64_t, 8>> cells;
    for (int k = 0; k < nz; ++k)
        for (auto &q : quads)
            cells.push_back({k * n2 + q[0], k * n2 + q[1], k * n2 + q[2], k * n2 + q[3], (k + 1) * n2 + q[0],
                             (k + 1) * n2 + q[1], (k + 1) * n2 + q[2], (k + 1) * n2 + q[3]});
    auto move = [R, lateralid](const std::vector<int> &markers, std::array<REAL, 3> &x) {
        if (std::find(markers.begin(), markers.end(), lateralid) == markers.end()) return;
        const REAL r = std::hypot(x[0], x[1]);
        x[0] *= R / r;
        x[1] *= R / r;
    };
    return internal::HexMesh(vertices, cells, matid, marker, 2, move);
}

/**
 * @brief Block \f$[x_0,x_1]\times[y_0,y_1]\times[z_0,z_1]\f$ divided into nx x ny x nz trilinear hexahedra
 * (TPZGeoCube), the generalization of CreateBoxMesh to any origin
 * @param x0 lower corner (x, y, z)
 * @param x1 upper corner (x, y, z)
 * @param nx,ny,nz number of elements in each direction
 * @param matid material id of the hexahedra
 * @param marker material ids of the boundary quadrilaterals (TPZGeoQuad) created on each boundary face
 * (several ids create coincident elements, none creates no element); FaceOnPlane helps to write it
 *
 * The vertices are numbered with x fastest, then y, then z: node (i, j, k) has index
 * \f$k(n_y+1)(n_x+1)+j(n_x+1)+i\f$ and coordinates \f$x_0+(x_1-x_0)i/n_x\f$, etc.; the elements follow the
 * same order, and the vertices of each element follow TPZCube (see WriteMeshCSV). Geometry is trilinear:
 * the 20-node serendipity displacement space (CreateDisplacementMesh) is built on it without geometric
 * mid-edge nodes.
 */
inline TPZGeoMesh *CreateBlockMesh(const std::array<REAL, 3> &x0, const std::array<REAL, 3> &x1, int nx, int ny,
                                   int nz, int matid, const TMarker3D &marker) {
    auto idx = [nx, ny](int i, int j, int k) { return int64_t(k) * (ny + 1) * (nx + 1) + int64_t(j) * (nx + 1) + i; };
    std::vector<std::array<REAL, 3>> vertices;
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i)
                vertices.push_back({x0[0] + (x1[0] - x0[0]) * i / nx, x0[1] + (x1[1] - x0[1]) * j / ny,
                                    x0[2] + (x1[2] - x0[2]) * k / nz});
    std::vector<std::array<int64_t, 8>> cells;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                cells.push_back({idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                                 idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)});
    return internal::HexMesh(vertices, cells, matid, marker, 1);
}

/**
 * @brief A single trilinear hexahedron on the unit cube \f$[0,1]^3\f$ (element tests with a homogeneous
 * state, e.g. the triaxial tests of RS2 and FLAC3D in 3D)
 * @param matid material id of the hexahedron
 * @param marker material ids of the boundary quadrilaterals of each of the six faces (see CreateBlockMesh)
 *
 * Nodes 0-7 are the vertices of TPZCube: (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1),
 * (0,1,1).
 */
inline TPZGeoMesh *CreateUnitCubeMesh(int matid, const TMarker3D &marker) {
    return CreateBlockMesh({0., 0., 0.}, {1., 1., 1.}, 1, 1, 1, matid, marker);
}

/**
 * @brief Slab \f$[x_0,x_1]\times[y_0,y_1]\times[0,t]\f$ with nx x ny x 1 trilinear hexahedra: the 3D
 * counterpart of the 2D structured mesh CreateRectangleMesh (plane strain is imposed with
 * \f$u_z=0\f$ on the faces z = 0 and z = t)
 * @param x0,y0 lower corner of the rectangle
 * @param x1,y1 upper corner of the rectangle
 * @param thickness thickness t of the slab (z direction)
 * @param nx,ny number of elements in x and y
 * @param matid material id of the hexahedra
 * @param marker material ids of the boundary quadrilaterals of each boundary face (see CreateBlockMesh)
 *
 * The vertices of the plane z = 0 have the indices of the nodes of CreateRectangleMesh (row by row,
 * \f$j(n_x+1)+i\f$) and those of the plane z = t follow with the offset \f$(n_x+1)(n_y+1)\f$.
 */
inline TPZGeoMesh *CreateSlabMesh(REAL x0, REAL y0, REAL x1, REAL y1, REAL thickness, int nx, int ny, int matid,
                                  const TMarker3D &marker) {
    return CreateBlockMesh({x0, y0, 0.}, {x1, y1, thickness}, nx, ny, 1, matid, marker);
}

/**
 * @brief Writes the geometric mesh as CSV files for plotting (e.g. the finite element meshes of the
 * figures of the article drawn with matplotlib)
 * @param gmesh geometric mesh (only the leaf elements, without subelements, are written)
 * @param prefix path prefix of the files
 *
 * Files (one header line, comma separated, integers for the indices, coordinates with 16 digits):
 *  - \<prefix\>_nodes.csv: @c node, @c x, @c y, @c z for every node of the mesh; @c node is the index
 *    of the node in TPZGeoMesh::NodeVec (the row order), referred to by the other files;
 *  - \<prefix\>_elements.csv: @c element, @c matid, @c nnodes, @c n0, ..., @c n{N-1} for every element of
 *    the dimension of the mesh (volume elements in 3D): element index in the geometric mesh, material id,
 *    number of nodes and node indices in the NeoPZ order of the element, padded with -1 up to the largest
 *    number of nodes N;
 *  - \<prefix\>_faces.csv: same columns for the elements of dimension one less (the boundary condition
 *    elements: quadrilaterals in 3D, lines in 2D). A face with several boundary conditions appears once per
 *    material id (coincident elements); the orientation of the faces is not normalized;
 *  - \<prefix\>_edges.csv: @c n0, @c n1, @c nmid, @c boundary for each distinct edge of the elements of
 *    the mesh dimension: end nodes (n0 < n1), mid-edge node of a quadratic geometry (-1 for a straight
 *    edge) and 1 if the edge belongs to an element of the faces file (0 otherwise), for wireframe plots of
 *    the whole mesh or of its boundary only.
 *
 * Node order of the hexahedra (TPZCube; TPZGeoCube has 8 nodes, TPZQuadraticCube 20). Vertices 0-7 at the
 * parametric points \f$(\xi,\eta,\zeta)\f$ = (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1), (-1,-1,1),
 * (1,-1,1), (1,1,1), (-1,1,1); nodes 8-19 at the middle of the edges (sides 8-19 of TPZCube)
 * 0-1, 1-2, 2-3, 3-0, 0-4, 1-5, 2-6, 3-7, 4-5, 5-6, 6-7, 7-4. Quadrilaterals (TPZQuadrilateral;
 * TPZGeoQuad 4 nodes, TPZQuadraticQuad 8): vertices 0-3 at (-1,-1), (1,-1), (1,1), (-1,1), nodes 4-7 at
 * the middle of the edges 0-1, 1-2, 2-3, 3-0. In general the quadratic maps of NeoPZ number the mid-edge
 * node of edge side s as node s, which is how the edges file is built. A curved edge of a
 * TPZQuadraticCube (e.g. the lateral face of CreateQuarterCylinderMesh) is the parabola through its end
 * nodes and its mid-edge node: x(t) = x0 (1-t)(1-2t) + 4 xm t(1-t) + x1 t(2t-1), t in [0, 1].
 * Elements of other dimensions (points, lines in 3D) are not written.
 */
inline void WriteMeshCSV(TPZGeoMesh *gmesh, const std::string &prefix) {
    const int dim = gmesh->Dimension();
    {
        std::ofstream out(prefix + "_nodes.csv");
        out << "node,x,y,z\n" << std::setprecision(16);
        TPZManVector<REAL, 3> x(3, 0.);
        for (int64_t i = 0; i < gmesh->NNodes(); ++i) {
            gmesh->NodeVec()[i].GetCoordinates(x);
            out << i << "," << x[0] << "," << x[1] << "," << x[2] << "\n";
        }
    }
    auto writeElements = [gmesh](const std::string &file, int eldim) {
        std::vector<TPZGeoEl *> els;
        int nmax = 0;
        for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel || gel->HasSubElement() || gel->Dimension() != eldim) continue;
            els.push_back(gel);
            nmax = std::max(nmax, gel->NNodes());
        }
        std::ofstream out(file);
        out << "element,matid,nnodes";
        for (int k = 0; k < nmax; ++k) out << ",n" << k;
        out << "\n";
        for (TPZGeoEl *gel : els) {
            out << gel->Index() << "," << gel->MaterialId() << "," << gel->NNodes();
            for (int k = 0; k < nmax; ++k) out << "," << (k < gel->NNodes() ? gel->NodeIndex(k) : int64_t(-1));
            out << "\n";
        }
    };
    writeElements(prefix + "_elements.csv", dim);
    writeElements(prefix + "_faces.csv", dim - 1);
    {
        std::map<std::pair<int64_t, int64_t>, int64_t> edges;
        std::set<std::pair<int64_t, int64_t>> boundary;
        for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel || gel->HasSubElement() || gel->Dimension() < dim - 1 || gel->Dimension() > dim) continue;
            const bool quadratic = gel->NNodes() > gel->NCornerNodes();
            for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side) {
                if (gel->SideDimension(side) != 1) continue;
                const int64_t a = gel->SideNodeIndex(side, 0), b = gel->SideNodeIndex(side, 1);
                const std::pair<int64_t, int64_t> key(std::min(a, b), std::max(a, b));
                if (gel->Dimension() == dim - 1) {
                    boundary.insert(key);
                    continue;
                }
                const int64_t mid = (quadratic && side < gel->NNodes()) ? gel->NodeIndex(side) : int64_t(-1);
                edges.insert({key, mid});
            }
        }
        std::ofstream out(prefix + "_edges.csv");
        out << "n0,n1,nmid,boundary\n";
        for (auto &e : edges)
            out << e.first.first << "," << e.first.second << "," << e.second << "," << int(boundary.count(e.first))
                << "\n";
    }
}
/** @} */

/** @name Computational meshes */
/** @{ */

/**
 * @brief H1 displacement mesh of order 2 on the volume material and its boundary materials.
 * With serendipity = true the face and volume connects of quadrilaterals and hexahedra are reduced
 * to order 1 (no internal shape functions): the space is then exactly the serendipity Q8/Hex20 space
 * of the article (the edge functions of NeoPZ are made linear across the element by adding the bubbles).
 */
inline TPZCompMesh *CreateDisplacementMesh(TPZGeoMesh *gmesh, int dim, int volid, const std::set<int> &bcids,
                                           bool serendipity = true, int order = 2) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("Displacement");
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(volid, dim, dim);
    cmesh->InsertMaterialObject(mat);
    TPZFNMatrix<9, STATE> val1(dim, dim, 0.);
    TPZManVector<STATE, 3> val2(dim, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    if (serendipity) {
        for (int64_t iel = 0; iel < cmesh->NElements(); ++iel) {
            auto *intel = dynamic_cast<TPZInterpolatedElement *>(cmesh->Element(iel));
            if (!intel) continue;
            TPZGeoEl *gel = intel->Reference();
            if (gel->Type() != EQuadrilateral && gel->Type() != ECube) continue;
            for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side)
                if (gel->SideDimension(side) >= 2) intel->SetSideOrder(side, 1);
        }
        cmesh->ExpandSolution();
    }
    return cmesh;
}

/** @brief H1 pore pressure mesh of order 1 on the volume material and its boundary materials */
inline TPZCompMesh *CreatePressureMesh(TPZGeoMesh *gmesh, int dim, int volid, const std::set<int> &bcids,
                                       int order = 1) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("PorePressure");
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    auto *mat = new TPZNullMaterial<STATE>(volid, dim, 1);
    cmesh->InsertMaterialObject(mat);
    TPZFNMatrix<1, STATE> val1(1, 1, 0.);
    TPZManVector<STATE, 1> val2(1, 0.);
    for (int id : bcids) cmesh->InsertMaterialObject(mat->CreateBC(mat, id, 0, val1, val2));
    cmesh->AutoBuild();
    return cmesh;
}

/**
 * @brief Builds the multiphysics space (displacement, pore pressure) with memory in the volume
 * material; the materials must already be inserted in mphys
 */
inline void BuildMultiphysics(TPZMultiphysicsCompMesh *mphys, TPZCompMesh *cmeshU, TPZCompMesh *cmeshP, int volid,
                              const std::set<int> &bcids) {
    TPZManVector<TPZCompMesh *, 2> meshvec(2);
    meshvec[0] = cmeshU;
    meshvec[1] = cmeshP;
    TPZManVector<int, 2> active(2, 1);
    std::set<int> withmem = {volid};
    mphys->BuildMultiphysicsSpaceWithMemory(active, meshvec, withmem, bcids);
}

/**
 * @brief Sets the nodal values of the pore pressure mesh (vertex connects) with f(x) and loads them in the
 * multiphysics mesh (initial pore pressure field)
 */
inline void SetInitialPressure(TPZMultiphysicsCompMesh *mphys, const std::function<REAL(const TPZVec<REAL> &x)> &f) {
    TPZCompMesh *cmeshP = mphys->MeshVector()[1];
    TPZFMatrix<STATE> &sol = cmeshP->Solution();
    for (int64_t iel = 0; iel < cmeshP->NElements(); ++iel) {
        TPZCompEl *cel = cmeshP->Element(iel);
        if (!cel || !cel->Reference()) continue;
        TPZGeoEl *gel = cel->Reference();
        for (int i = 0; i < gel->NCornerNodes() && i < cel->NConnects(); ++i) {
            TPZConnect &c = cel->Connect(i);
            const int64_t seq = c.SequenceNumber();
            if (seq < 0 || cmeshP->Block().Size(seq) == 0) continue;
            TPZManVector<REAL, 3> x(3);
            gel->NodePtr(i)->GetCoordinates(x);
            sol(cmeshP->Block().Position(seq), 0) = f(x);
        }
    }
    mphys->LoadSolutionFromMeshes();
}
/** @} */

/** @name Post-processing */
/** @{ */

/** @brief Converged state of an integration point */
struct TGaussPoint {
    int64_t fElement;
    int fIndex;
    TPZManVector<REAL, 3> fX;
    TPZManVector<REAL, 3> fQsi;
    TPZTensor<REAL> fSigma; ///< effective stress
    TPZTensor<REAL> fEps;   ///< total strain (engineering shear components)
    REAL fPc;
    REAL fV0;
    int fType;
};

/** @brief Collects the converged state of all the integration points of the u-p material */
inline std::vector<TGaussPoint> GaussPoints(TPoroMaterial *mat, TPZCompMesh *mesh) {
    std::vector<TGaussPoint> out;
    mat->ForEachIntegrationPoint(mesh, [&](TPZCompEl *cel, int ip, const TPZVec<REAL> &x, const TPZVec<REAL> &qsi,
                                           REAL, TPZElastoPlasticMem &mem) {
        TGaussPoint g;
        g.fElement = cel->Index();
        g.fIndex = ip;
        g.fX = x;
        g.fQsi = qsi;
        g.fSigma = mem.m_sigma;
        g.fEps = mem.m_elastoplastic_state.m_eps_t;
        g.fPc = mem.m_elastoplastic_state.m_hardening;
        g.fV0 = mem.m_elastoplastic_state.fmatprop.size() ? mem.m_elastoplastic_state.fmatprop[0] : 0.;
        g.fType = mem.m_elastoplastic_state.m_m_type;
        out.push_back(g);
    });
    return out;
}

/** @brief Number of integration points in the elastic (0), subcritical (1) and supercritical (2) states */
inline std::array<int, 3> CountTypes(TPoroMaterial *mat, TPZCompMesh *mesh) {
    std::array<int, 3> n = {0, 0, 0};
    for (auto &g : GaussPoints(mat, mesh))
        if (g.fType >= 0 && g.fType < 3) n[g.fType]++;
    return n;
}

/**
 * @brief Writes the integration points as a VTK point cloud (legacy ASCII) with the effective stress,
 * p', q, preconsolidation pressure, specific volume and type of response (0 elastic, 1 subcritical,
 * 2 supercritical)
 */
inline void WriteGaussPointsVTK(TPoroMaterial *mat, TPZCompMesh *mesh, const std::string &file) {
    auto gps = GaussPoints(mat, mesh);
    std::ofstream out(file);
    out << "# vtk DataFile Version 3.0\nintegration points\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << gps.size() << " double\n" << std::setprecision(10);
    for (auto &g : gps) out << g.fX[0] << " " << g.fX[1] << " " << g.fX[2] << "\n";
    out << "CELLS " << gps.size() << " " << 2 * gps.size() << "\n";
    for (size_t i = 0; i < gps.size(); ++i) out << "1 " << i << "\n";
    out << "CELL_TYPES " << gps.size() << "\n";
    for (size_t i = 0; i < gps.size(); ++i) out << "1\n";
    out << "POINT_DATA " << gps.size() << "\n";
    auto scalar = [&](const char *name, const std::function<REAL(const TGaussPoint &)> &f) {
        out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        for (auto &g : gps) out << f(g) << "\n";
    };
    scalar("MeanEffectiveStress", [](const TGaussPoint &g) { return MeanEffectiveStress(g.fSigma); });
    scalar("DeviatoricStress", [](const TGaussPoint &g) { return DeviatoricStress(g.fSigma); });
    scalar("PreconsolidationPressure", [](const TGaussPoint &g) { return g.fPc; });
    scalar("SpecificVolume", [](const TGaussPoint &g) { return g.fV0; });
    scalar("PlasticType", [](const TGaussPoint &g) { return REAL(g.fType); });
    out << "TENSORS EffectiveStress double\n";
    for (auto &g : gps) {
        for (int r = 0; r < 3; ++r) out << g.fSigma(r, 0) << " " << g.fSigma(r, 1) << " " << g.fSigma(r, 2) << "\n";
    }
}

/**
 * @brief VTK output of the nodal fields (displacement and pore pressure) with the native graphical mesh
 * of NeoPZ (TPZAnalysis::DefineGraphMesh/PostProcess); step is the index of the file
 */
inline void WriteNodalVTK(TPZLinearAnalysis &an, int dim, const std::string &file, int step, int resolution = 0) {
    TPZManVector<std::string, 2> scal(1), vec(1);
    scal[0] = "PorePressure";
    vec[0] = "Displacement";
    an.DefineGraphMesh(dim, scal, vec, file);
    an.SetStep(step);
    an.PostProcess(resolution, dim);
}

/**
 * @brief Volume element that contains the point x and its parametric coordinates (LocatePoint): the
 * first element in the order of the geometric mesh (as in the Python code), checked by bounding box and
 * inversion of the geometric map
 */
inline TPZGeoEl *LocatePoint(TPZGeoMesh *gmesh, int volid, const TPZVec<REAL> &xpoint, TPZVec<REAL> &qsi) {
    const int dim = gmesh->Dimension();
    TPZManVector<REAL, 3> x(3, 0.);
    for (int k = 0; k < xpoint.size() && k < 3; ++k) x[k] = xpoint[k];
    for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel || gel->MaterialId() != volid || gel->Dimension() != dim || gel->HasSubElement()) continue;
        // bounding box of the nodes
        bool inside = true;
        for (int c = 0; c < dim && inside; ++c) {
            REAL lo = 1e300, hi = -1e300;
            for (int i = 0; i < gel->NNodes(); ++i) {
                const REAL v = gel->NodePtr(i)->Coord(c);
                lo = std::min(lo, v);
                hi = std::max(hi, v);
            }
            if (x[c] < lo - 1e-12 || x[c] > hi + 1e-12) inside = false;
        }
        if (!inside) continue;
        qsi.Resize(dim);
        qsi.Fill(0.);
        gel->ComputeXInverse(x, qsi, 1.e-12);
        bool ok = true;
        for (int c = 0; c < dim; ++c)
            if (std::fabs(qsi[c]) > 1. + 1e-7) ok = false;
        if (ok) return gel;
    }
    return nullptr;
}

/**
 * @brief Effective stress at a point interpolated from the integration points of the element that
 * contains it with Lagrange polynomials on the Gauss abscissas (LocatePoint and GaussPointWeights of
 * fe_user.py); valid for the tensor-product Gauss rules of quadrilaterals and hexahedra
 */
inline TPZTensor<REAL> StressAtPoint(TPoroMaterial *mat, TPZCompMesh *mesh, const TPZVec<REAL> &xpoint) {
    TPZGeoMesh *gmesh = mesh->Reference();
    const int dim = gmesh->Dimension();
    TPZManVector<REAL, 3> qsi(dim, 0.);
    TPZGeoEl *gel = LocatePoint(gmesh, mat->Id(), xpoint, qsi);
    if (!gel) DebugStop();
    std::vector<TGaussPoint> pts;
    for (auto &g : GaussPoints(mat, mesh)) {
        if (mesh->Element(g.fElement)->Reference() == gel) pts.push_back(g);
    }
    if (pts.empty()) DebugStop();
    std::vector<REAL> absc;
    for (auto &g : pts) absc.push_back(g.fQsi[0]);
    std::sort(absc.begin(), absc.end());
    absc.erase(std::unique(absc.begin(), absc.end(), [](REAL a, REAL b) { return std::fabs(a - b) < 1e-10; }),
               absc.end());
    auto lag = [&](int l, REAL t) {
        REAL v = 1.;
        for (size_t m = 0; m < absc.size(); ++m)
            if ((int)m != l) v *= (t - absc[m]) / (absc[l] - absc[m]);
        return v;
    };
    auto index = [&](REAL t) {
        int best = 0;
        for (size_t m = 0; m < absc.size(); ++m)
            if (std::fabs(absc[m] - t) < std::fabs(absc[best] - t)) best = m;
        return best;
    };
    TPZTensor<REAL> sig;
    for (auto &g : pts) {
        REAL w = 1.;
        for (int c = 0; c < dim; ++c) w *= lag(index(g.fQsi[c]), qsi[c]);
        sig.Add(g.fSigma, w);
    }
    return sig;
}

/** @brief Writes a table (header and rows) as comma separated values */
inline void WriteCSV(const std::string &file, const std::vector<std::string> &header,
                     const std::vector<std::vector<REAL>> &rows) {
    std::ofstream out(file);
    for (size_t i = 0; i < header.size(); ++i) out << header[i] << (i + 1 < header.size() ? "," : "\n");
    out << std::setprecision(12);
    for (auto &r : rows)
        for (size_t i = 0; i < r.size(); ++i) out << r[i] << (i + 1 < r.size() ? "," : "\n");
}

/**
 * @brief Restores the references of a geometric mesh when it goes out of scope
 *
 * TPZPostProcAnalysis (SetCompMesh, SetPostProcessVariables, TransferSolution and its destructor) resets and
 * loads the references of the geometric mesh (TPZGeoEl::Reference and TPZGeoMesh::Reference). The guard saves
 * them at construction and restores them at destruction, so that the post-processing leaves the geometric
 * mesh exactly as it found it.
 */
class TReferenceGuard {
public:
    /** @brief Saves the references of the elements of gmesh and of gmesh itself */
    explicit TReferenceGuard(TPZGeoMesh *gmesh) : fGMesh(gmesh), fMesh(gmesh->Reference()) {
        fRefs.resize(gmesh->NElements(), nullptr);
        for (int64_t i = 0; i < gmesh->NElements(); ++i)
            if (TPZGeoEl *gel = gmesh->Element(i)) fRefs[i] = gel->Reference();
    }
    /** @brief Restores the saved references */
    ~TReferenceGuard() {
        fGMesh->SetReference(fMesh);
        for (int64_t i = 0; i < fGMesh->NElements() && i < int64_t(fRefs.size()); ++i)
            if (TPZGeoEl *gel = fGMesh->Element(i)) gel->SetReference(fRefs[i]);
    }
    TReferenceGuard(const TReferenceGuard &) = delete;
    TReferenceGuard &operator=(const TReferenceGuard &) = delete;

private:
    TPZGeoMesh *fGMesh;
    TPZCompMesh *fMesh;
    std::vector<TPZCompEl *> fRefs;
};

/**
 * @brief Writes a ParaView file series (JSON, "\<file\>.vtk.series")
 * @param file name of the series file
 * @param names names of the data files, relative to the directory of the series file
 * @param times time of each data file (shown by ParaView in the time controls)
 */
inline void WriteFileSeries(const std::string &file, const std::vector<std::string> &names,
                            const std::vector<REAL> &times) {
    std::ofstream out(file);
    out << "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n" << std::setprecision(15);
    for (size_t i = 0; i < names.size() && i < times.size(); ++i)
        out << "    { \"name\" : \"" << names[i] << "\", \"time\" : " << times[i] << " }"
            << (i + 1 < names.size() ? "," : "") << "\n";
    out << "  ]\n}\n";
}

/**
 * @brief VTK output of every converged state of a u-p run, as ParaView file series
 *
 * Built once per run, after the multiphysics mesh and the memory of the integration points are initialized,
 * and destroyed before the meshes. Each call of Write (one per converged state k = 0, 1, ...) writes in the
 * directory of the run (created with std::filesystem):
 *  - \<name\>_nodal.scal_vec.\<k\>.vtk: Displacement (vector) and PorePressure at the nodes, written by the native
 *    graph mesh of the multiphysics mesh (TPZAnalysis::DefineGraphMesh once, SetStep and PostProcess for each
 *    state); the graph mesh belongs to an analysis object of its own, so that it is independent of other
 *    post-processing calls on the analysis that solves the problem (e.g. WriteNodalVTK);
 *  - \<name\>_intpoints.scal_vec.\<k\>.vtk: the variables of the integration points of TPZMatPoroElastoPlasticUP
 *    (IntegrationPointVariables) projected element by element on a discontinuous mesh by TPZPostProcAnalysis,
 *    as in the footing example of NeoPZ (SetPostProcessVariables once, TransferSolution and PostProcess for
 *    each state). The order of the projection is n-1 for n x n (x n) Gauss points, so the nodal values are
 *    the Lagrange extrapolation of the values at the points (the interpolation of StressAtPoint). Being
 *    extrapolated, the vertex values can overshoot: PlasticType is not an integer there and the principal
 *    stresses, sorted at the points, can be out of order where two of them are close;
 *  - \<name\>_gausspoints.\<k\>.vtk: the integration points as a point cloud (WriteGaussPointsVTK), optional;
 *  - \<name\>_nodal.vtk.series, \<name\>_intpoints.vtk.series and \<name\>_gausspoints.vtk.series: ParaView file
 *    series with the time of each file (rewritten after each state, so that an interrupted run is readable);
 *  - \<name\>_states.csv: index k of the state, series time and the values given to Write.
 *
 * The output does not change the analysis: the solution and the memory are only read, and the references
 * of the geometric mesh are restored after each operation (TReferenceGuard).
 */
class TVTKSeries {
public:
    /**
     * @brief Creates the directory, the graph mesh of the nodal fields and the post-processing mesh
     * @param mphys multiphysics mesh of the run
     * @param mat u-p material (its id selects the elements of the post-processing mesh)
     * @param dir directory of the files (created if it does not exist)
     * @param name prefix of the file names
     * @param timename name of the series time in the CSV file of the states
     * @param columns names of the values given to Write (columns of the CSV file of the states)
     * @param gausspoints also write the point cloud of the integration points
     * @param resolution resolution of the graph meshes (0: vertices of the elements; r: 2^r subdivisions per side)
     */
    TVTKSeries(TPZMultiphysicsCompMesh *mphys, TPoroMaterial *mat, const std::string &dir, const std::string &name,
               const std::string &timename, const std::vector<std::string> &columns = {}, bool gausspoints = true,
               int resolution = 0)
        : fMPhys(mphys), fMat(mat), fDir(dir), fName(name), fTimeName(timename), fColumns(columns),
          fGaussPoints(gausspoints), fResolution(resolution), fDim(mat->Dimension()) {
        std::filesystem::create_directories(fDir);
        TReferenceGuard guard(mphys->Reference());
        // nodal fields: native graph mesh of the multiphysics mesh (no renumbering of the equations)
        fNodal = std::make_unique<TPZLinearAnalysis>(mphys, false);
        TPZManVector<std::string, 1> scal(1, "PorePressure"), vec(1, "Displacement");
        fNodal->DefineGraphMesh(fDim, scal, vec, Path("_nodal.vtk"));
        // integration point fields: post-processing mesh of TPZPostProcAnalysis (footing example)
        TPZManVector<std::string, 32> scalip, vecip, tensip, all;
        IntegrationPointVariables(fDim, scalip, vecip, tensip);
        for (auto *names : {&scalip, &vecip, &tensip})
            for (int i = 0; i < names->size(); ++i) {
                all.Resize(all.size() + 1);
                all[all.size() - 1] = (*names)[i];
            }
        TPZManVector<int, 1> matids(1, mat->Id());
        fPost = std::make_unique<TPZPostProcAnalysis>();
        fPost->SetCompMesh(mphys);
        fPost->SetPostProcessVariables(matids, all);
        TPZFStructMatrix<STATE> sm(fPost->Mesh());
        sm.SetNumThreads(0);
        fPost->SetStructuralMatrix(sm);
        fPost->DefineGraphMesh(fDim, scalip, vecip, tensip, Path("_intpoints.vtk"));
    }

    /** @brief Deletes the post-processing analyses (before the meshes) and restores the references */
    ~TVTKSeries() {
        TReferenceGuard guard(fMPhys->Reference());
        fPost.reset();
        fNodal.reset();
    }

    TVTKSeries(const TVTKSeries &) = delete;
    TVTKSeries &operator=(const TVTKSeries &) = delete;

    /**
     * @brief Writes the files of the current state (solution loaded in the mesh, memory of the converged step)
     * @param time time of the state in the file series
     * @param values values of the columns of the CSV file of the states
     */
    void Write(REAL time, const std::vector<REAL> &values = {}) {
        const int k = static_cast<int>(fTimes.size());
        {
            TReferenceGuard guard(fMPhys->Reference());
            fNodal->SetStep(k);
            fNodal->PostProcess(fResolution, fDim);
            fPost->TransferSolution();
            fPost->SetStep(k);
            fPost->PostProcess(fResolution, fDim);
        }
        if (fGaussPoints) WriteGaussPointsVTK(fMat, fMPhys, fDir + "/" + GaussPointsFile(k));
        fTimes.push_back(time);
        fValues.push_back(values);
        // file series and table of the states, rewritten after each state
        std::vector<std::string> nodal, intpoints, gauss;
        for (int i = 0; i <= k; ++i) {
            nodal.push_back(fName + "_nodal.scal_vec." + std::to_string(i) + ".vtk");
            intpoints.push_back(fName + "_intpoints.scal_vec." + std::to_string(i) + ".vtk");
            gauss.push_back(GaussPointsFile(i));
        }
        WriteFileSeries(Path("_nodal.vtk.series"), nodal, fTimes);
        WriteFileSeries(Path("_intpoints.vtk.series"), intpoints, fTimes);
        if (fGaussPoints) WriteFileSeries(Path("_gausspoints.vtk.series"), gauss, fTimes);
        std::vector<std::string> header = {"index", fTimeName};
        header.insert(header.end(), fColumns.begin(), fColumns.end());
        std::vector<std::vector<REAL>> rows;
        for (int i = 0; i <= k; ++i) {
            rows.push_back({REAL(i), fTimes[i]});
            rows.back().insert(rows.back().end(), fValues[i].begin(), fValues[i].end());
        }
        WriteCSV(Path("_states.csv"), header, rows);
    }

    /** @brief Number of states written */
    int NStates() const { return static_cast<int>(fTimes.size()); }

    /** @brief Directory of the files */
    const std::string &Directory() const { return fDir; }

    /**
     * @brief Variables of TPZMatPoroElastoPlasticUP written at the integration points: scalars (p', q, pc, type of
     * response, volumetric strain, specific volume and the components of the effective and total stresses;
     * XZ and YZ only in 3D), vectors (principal effective stresses) and tensors (effective and total stress)
     */
    static void IntegrationPointVariables(int dim, TPZVec<std::string> &scal, TPZVec<std::string> &vec,
                                          TPZVec<std::string> &tens) {
        std::vector<std::string> s = {"MeanEffectiveStress", "DeviatoricStress", "PreconsolidationPressure",
                                      "PlasticType", "VolumetricStrain", "SpecificVolume"};
        std::vector<std::string> comps = {"XX", "YY", "ZZ", "XY"};
        if (dim == 3) comps.insert(comps.end(), {"XZ", "YZ"});
        for (const char *stress : {"EffectiveStress", "TotalStress"})
            for (auto &c : comps) s.push_back(stress + c);
        scal.Resize(s.size());
        for (size_t i = 0; i < s.size(); ++i) scal[i] = s[i];
        vec.Resize(1);
        vec[0] = "PrincipalEffectiveStress";
        tens.Resize(2);
        tens[0] = "EffectiveStress";
        tens[1] = "TotalStress";
    }

private:
    /** @brief Path of a file of the run: directory, prefix and suffix */
    std::string Path(const std::string &suffix) const { return fDir + "/" + fName + suffix; }
    /** @brief Name of the point cloud of the state k */
    std::string GaussPointsFile(int k) const { return fName + "_gausspoints." + std::to_string(k) + ".vtk"; }

    TPZMultiphysicsCompMesh *fMPhys;
    TPoroMaterial *fMat;
    std::string fDir, fName, fTimeName;
    std::vector<std::string> fColumns;
    bool fGaussPoints;
    int fResolution;
    int fDim;
    std::unique_ptr<TPZLinearAnalysis> fNodal;  ///< analysis of the graph mesh of the nodal fields
    std::unique_ptr<TPZPostProcAnalysis> fPost; ///< projection of the integration point variables
    std::vector<REAL> fTimes;                   ///< series time of each state
    std::vector<std::vector<REAL>> fValues;     ///< values of each state for the CSV file
};

/** @brief Deletes a multiphysics mesh, its atomic meshes and the geometric mesh */
inline void DeleteMeshes(TPZMultiphysicsCompMesh *mphys) {
    TPZGeoMesh *gmesh = mphys->Reference();
    TPZManVector<TPZCompMesh *, 2> meshvec(mphys->MeshVector());
    delete mphys;
    for (auto *m : meshvec) delete m;
    delete gmesh;
}

/** @brief Mean number of residual evaluations per converged increment of a step log */
inline REAL MeanEvaluations(const std::vector<TAnalysis::TStepLog> &log, size_t first = 0) {
    if (log.size() <= first) return 0.;
    REAL sum = 0.;
    for (size_t i = first; i < log.size(); ++i) sum += log[i].fResiduals.size();
    return sum / (log.size() - first);
}
/** @} */

/** @name Material point drivers (TriaxialPointCC) */
/** @{ */

/** @brief Statistics of the local Newton iterations of the projections */
struct TLocalStats {
    int fCalls = 0; ///< number of plastic projections
    int fIts = 0;   ///< total number of local iterations
    int fMax = 0;   ///< maximum number of local iterations
    REAL Mean() const { return fCalls ? REAL(fIts) / fCalls : 0.; }
};

/** @brief State of a material point: stress, preconsolidation pressure and specific volume */
struct TPointState {
    TPZTensor<REAL> fSigma; ///< effective stress
    REAL fPc;               ///< preconsolidation pressure
    REAL fV0;               ///< specific volume
    TPointState(const TPZTensor<REAL> &sigma, REAL pc, REAL v0) : fSigma(sigma), fPc(pc), fV0(v0) {}
};

/**
 * @brief Applies a total strain eps from the converged state (eps_n, state) and returns the stress and
 * the tangent; the plastic state of the model is not changed (ApplyStrainComputeSigmaDepCC)
 * @return false if the local projection failed
 */
inline bool ApplyStrain(const TPlastic &model, const TPZTensor<REAL> &epsn, const TPointState &sn,
                        const TPZTensor<REAL> &eps, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep, REAL &pc,
                        int &type, TLocalStats *stats = nullptr) {
    TPlastic m(model);
    TPZPlasticState<REAL> st = m.GetState();
    st.m_eps_t = epsn;
    st.m_hardening = sn.fPc;
    m.SetState(st);
    m.SetSpecificVolume(sn.fV0);
    sigma = sn.fSigma;
    m.ApplyStrainComputeSigma(eps, sigma, &Dep);
    if (m.LastProjectionFailed()) return false;
    pc = m.GetState().m_hardening;
    type = m.GetState().m_m_type;
    if (stats && type != 0) {
        stats->fCalls++;
        stats->fIts += m.LastNewtonIterations();
        stats->fMax = std::max(stats->fMax, m.LastNewtonIterations());
    }
    return true;
}

/**
 * @brief Drained triaxial test at a material point (TriaxialPointCC): sigma_xx = sigma_yy kept equal to
 * -p0, axial strain eps_zz prescribed in nsteps equal increments up to eamax (compression); the lateral
 * strains follow from a Newton iteration with the xx-yy block of the consistent tangent.
 * @return rows (eps_a, p', q, eps_v, eps_q, sigma_a) with compression positive; empty on failure
 */
inline std::vector<std::array<REAL, 6>> TriaxialDrained(const TPlastic &model, REAL p0, REAL pc0, REAL v0, REAL eamax,
                                                        int nsteps, TLocalStats *stats = nullptr) {
    std::vector<std::array<REAL, 6>> out;
    TPointState st(IsotropicTensor(-p0), pc0, v0);
    TPZTensor<REAL> eps;
    const REAL sig3 = st.fSigma.XX();
    out.push_back({0., p0, 0., 0., 0., p0});
    const REAL dea = eamax / nsteps;
    REAL er = 0.3 * dea;
    TPZFNMatrix<36, REAL> D(6, 6, 0.);
    for (int k = 0; k < nsteps; ++k) {
        TPZTensor<REAL> epst(eps);
        epst.ZZ() -= dea;
        epst.XX() += er;
        epst.YY() += er;
        TPZTensor<REAL> s;
        REAL pc = pc0;
        int type = 0;
        for (int it = 0; it < 30; ++it) {
            if (!ApplyStrain(model, eps, st, epst, s, D, pc, type, stats)) return {};
            const REAL r0 = s.XX() - sig3, r1 = s.YY() - sig3;
            if (std::max(std::fabs(r0), std::fabs(r1)) < 1e-10 * std::fabs(sig3)) break;
            const REAL a = D(_XX_, _XX_), b = D(_XX_, _YY_), c = D(_YY_, _XX_), d = D(_YY_, _YY_);
            const REAL det = a * d - b * c;
            epst.XX() += (-r0 * d + r1 * b) / det;
            epst.YY() += (-r1 * a + r0 * c) / det;
        }
        er = epst.XX() - eps.XX();
        eps = epst;
        st.fSigma = s;
        st.fPc = pc;
        out.push_back({-eps.ZZ(), MeanEffectiveStress(s), DeviatoricStress(s), -eps.I1(),
                       2. / 3. * std::fabs(eps.ZZ() - eps.XX()), -s.ZZ()});
    }
    return out;
}

/**
 * @brief Undrained triaxial test at a material point: strain increments with
 * \f$\Delta\varepsilon_{xx}=\Delta\varepsilon_{yy}=-\Delta\varepsilon_{zz}/2\f$ (eps_v = 0)
 * @return rows (eps_a, p', q); empty on failure
 */
inline std::vector<std::array<REAL, 3>> TriaxialUndrained(const TPlastic &model, REAL p0, REAL pc0, REAL v0,
                                                          REAL eamax, int nsteps, TLocalStats *stats = nullptr) {
    std::vector<std::array<REAL, 3>> out;
    TPointState st(IsotropicTensor(-p0), pc0, v0);
    TPZTensor<REAL> eps;
    out.push_back({0., p0, 0.});
    TPZFNMatrix<36, REAL> D(6, 6, 0.);
    for (int k = 0; k < nsteps; ++k) {
        TPZTensor<REAL> epst(eps);
        epst.XX() += 0.5 * eamax / nsteps;
        epst.YY() += 0.5 * eamax / nsteps;
        epst.ZZ() -= eamax / nsteps;
        TPZTensor<REAL> s;
        REAL pc;
        int type;
        if (!ApplyStrain(model, eps, st, epst, s, D, pc, type, stats)) return {};
        eps = epst;
        st.fSigma = s;
        st.fPc = pc;
        out.push_back({-eps.ZZ(), MeanEffectiveStress(s), DeviatoricStress(s)});
    }
    return out;
}
/** @} */

/** @name Closed-form solutions (Appendix B) */
/** @{ */

/**
 * @brief Drained triaxial test with constant cell pressure, eqs. (B.1)-(B.6) (TriaxialDrainedClosedCC);
 * pt = 0, omega = 1, porous elasticity, G constant (G > 0) or from the Poisson ratio nu (G <= 0)
 * @return rows (eps_a, p', q, eps_v, eps_q, sigma_a), compression positive
 */
inline std::vector<std::array<REAL, 6>> TriaxialDrainedClosed(REAL p0, REAL pc0, REAL v0, REAL M, REAL lambda,
                                                              REAL kappa, REAL G, REAL nu, int npts = 300) {
    const REAL r = 3. * (1. - 2. * nu) / (2. * (1. + nu));
    auto F = [M](REAL x) {
        return (1. / M) * std::log(std::fabs((M + x) / (M - x))) - (2. / M) * std::atan(x / M) -
               std::log(std::fabs(M - x)) / (3. - M) - std::log(M + x) / (3. + M) + 6. * std::log(3. - x) / (9. - M * M);
    };
    auto el = [&](REAL pp, REAL q, REAL &ev, REAL &eq) {
        ev = kappa / v0 * std::log(pp / p0);
        eq = (G > 0.) ? q / (3. * G) : kappa / (r * v0) * std::log(pp / p0);
    };
    const REAL A = 9. + M * M, B = -(18. * p0 + M * M * pc0), C = 9. * p0 * p0;
    const REAL disc = std::sqrt(B * B - 4. * A * C);
    REAL py = 1e300;
    for (REAL root : {(-B + disc) / (2. * A), (-B - disc) / (2. * A)})
        if (root >= p0 - 1e-9) py = std::min(py, root);
    const REAL qy = 3. * (py - p0), etay = qy / py;
    std::vector<std::array<REAL, 4>> rows;
    if (py - p0 > 1e-9 * p0) {
        for (int i = 0; i <= 40; ++i) {
            const REAL pp = p0 + (py - p0) * i / 40.;
            REAL ev, eq;
            el(pp, 3. * (pp - p0), ev, eq);
            rows.push_back({pp, 3. * (pp - p0), ev, eq});
        }
    } else {
        rows.push_back({p0, 0., 0., 0.});
    }
    std::vector<REAL> etas;
    for (int i = 1; i < npts; ++i) etas.push_back(etay + (M - etay) * REAL(i) / npts);
    for (int i = 1; i <= npts; ++i) etas.push_back(M - (M - etay) * std::exp(-12. * REAL(i) / npts));
    std::sort(etas.begin(), etas.end());
    etas.erase(std::unique(etas.begin(), etas.end()), etas.end());
    if (etay > M) std::reverse(etas.begin(), etas.end());
    for (REAL eta : etas) {
        const REAL pp = 3. * p0 / (3. - eta), q = eta * pp, pc = pp * (1. + eta * eta / (M * M));
        REAL ev, eq;
        el(pp, q, ev, eq);
        ev += (lambda - kappa) / v0 * std::log(pc / pc0);
        eq += (lambda - kappa) / v0 * (F(eta) - F(etay));
        rows.push_back({pp, q, ev, eq});
    }
    std::vector<std::array<REAL, 6>> out;
    for (auto &w : rows) out.push_back({w[3] + w[2] / 3., w[0], w[1], w[2], w[3], w[0] + 2. * w[1] / 3.});
    return out;
}

/** @brief Linear interpolation of column c of a table at the abscissa x of column 0 (numpy.interp) */
template <size_t N>
inline REAL Interpolate(const std::vector<std::array<REAL, N>> &tab, REAL x, int c) {
    if (tab.empty()) return 0.;
    if (x <= tab.front()[0]) return tab.front()[c];
    if (x >= tab.back()[0]) return tab.back()[c];
    // binary search of the interval and the formula of numpy.interp
    const auto it = std::upper_bound(tab.begin(), tab.end(), x,
                                     [](REAL v, const std::array<REAL, N> &row) { return v < row[0]; });
    const auto &a = *(it - 1), &b = *it;
    const REAL slope = (b[c] - a[c]) / (b[0] - a[0]);
    return slope * (x - a[0]) + a[c];
}

/**
 * @brief Undrained stress path (B.7): \f$p'_0/p' = ((M^2+\eta^2)/(M^2R))^\Lambda\f$, \f$\Lambda=(\lambda-\kappa)/\lambda\f$
 * @return p' for the stress ratio eta
 */
inline REAL UndrainedClosedP(REAL p0, REAL R, REAL M, REAL lambda, REAL kappa, REAL eta) {
    const REAL Lam = (lambda - kappa) / lambda;
    return p0 / std::pow((M * M + eta * eta) / (M * M * R), Lam);
}

/**
 * @brief Terzaghi's series (B.9), 200 terms
 * @param z distance from the drained face divided by the thickness H
 * @param T time factor \f$c_vt/H^2\f$
 * @param[out] U degree of consolidation
 * @return normalized excess pore pressure \f$p_w/q\f$
 */
inline REAL TerzaghiPressure(REAL z, REAL T, REAL &U) {
    REAL p = 0., s = 0.;
    for (int m = 0; m < 200; ++m) {
        const REAL mu = (2. * m + 1.) * M_PI / 2.;
        p += 2. / mu * std::sin(mu * z) * std::exp(-mu * mu * T);
        s += 2. / (mu * mu) * std::exp(-mu * mu * T);
    }
    U = 1. - s;
    return p;
}
/** @} */

} // namespace mcc

#endif
