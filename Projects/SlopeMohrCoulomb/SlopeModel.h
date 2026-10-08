// Slope of SlopeMohrCoulomb (moved from main.cpp, unchanged, so that other projects reuse it):
// soil data, Mohr-Coulomb models, geometric mesh and elastoplastic computational mesh.
#ifndef SLOPEMODEL_H
#define SLOPEMODEL_H

#include "Plasticity/TPZMatElastoPlastic2D.h"
#include "Plasticity/TPZPlasticStepPV.h"
#include "Plasticity/TPZPlasticStepVoigt.h"
#include "Plasticity/TPZYCMohrCoulombPV.h"
#include "Plasticity/TPZYCMohrCoulombPV2.h"
#include "TPZGeoLinear.h"
#include "pzcmesh.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"

#include <cmath>
#include <vector>

using TMCVoigt = TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>;
using TMCPV = TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;

struct Soil {
    REAL E = 20000., nu = 0.49, c = 10., phi = 30. * M_PI / 180., gamma = 20.;
};

inline TMCVoigt ModelVoigt(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    return TMCVoigt(TPZYCMohrCoulombPV2(s.phi, s.phi, s.c, ER), ER);
}

inline TMCPV ModelPV(const Soil &s) {
    TPZElasticResponse ER;
    ER.SetEngineeringData(s.E, s.nu);
    TMCPV pv;
    pv.fYC.SetUp(s.phi, s.phi, s.c, ER);
    pv.fER = ER;
    return pv;
}

/// Slope 70 x 40 m, height 10 m at 45 degrees. BC ids: -1 base, -2 right, -3 top right,
/// -4 top left, -5 left, -6 face. ref = uniform refinements.
inline TPZGeoMesh *TriGMesh(int ref) {
    auto *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    const std::vector<std::vector<REAL>> co = {
        {0, 0}, {10, 0}, {20, 0}, {30, 0}, {40, 0}, {50, 0}, {60, 0}, {70, 0},
        {0, 10}, {10, 10}, {20, 10}, {30, 10}, {40, 10}, {50, 10}, {60, 10}, {70, 10},
        {0, 20}, {10, 20}, {20, 20}, {30, 20}, {40, 20}, {50, 20}, {60, 20}, {70, 20},
        {0, 30}, {10, 30}, {20, 30}, {30, 30}, {40, 30}, {50, 30}, {60, 30}, {70, 30},
        {0, 40}, {10, 40}, {20, 40}, {30, 40}};
    const std::vector<std::vector<int64_t>> tri = {
        {0, 1, 8}, {1, 9, 8}, {1, 2, 9}, {2, 10, 9}, {2, 3, 10}, {3, 11, 10}, {3, 4, 11},
        {4, 12, 11}, {4, 5, 12}, {5, 13, 12}, {5, 6, 13}, {6, 14, 13}, {6, 7, 14}, {7, 15, 14},
        {8, 9, 16}, {9, 17, 16}, {9, 10, 17}, {10, 18, 17}, {10, 11, 18}, {11, 19, 18}, {11, 12, 19},
        {12, 20, 19}, {12, 13, 20}, {13, 21, 20}, {13, 14, 21}, {14, 22, 21}, {14, 15, 22}, {15, 23, 22},
        {16, 17, 24}, {17, 25, 24}, {17, 18, 25}, {18, 26, 25}, {18, 19, 26}, {19, 27, 26}, {19, 20, 27},
        {20, 28, 27}, {20, 21, 28}, {21, 29, 28}, {21, 22, 29}, {22, 30, 29}, {22, 23, 30}, {23, 31, 30},
        {24, 25, 32}, {25, 33, 32}, {25, 26, 33}, {26, 34, 33}, {26, 27, 34}, {27, 35, 34}, {27, 28, 35}};
    const std::vector<std::pair<std::vector<int64_t>, int>> lines = {
        {{0, 1}, -1}, {{1, 2}, -1}, {{2, 3}, -1}, {{3, 4}, -1}, {{4, 5}, -1}, {{5, 6}, -1}, {{6, 7}, -1},
        {{7, 15}, -2}, {{15, 23}, -2}, {{23, 31}, -2},
        {{31, 30}, -3}, {{30, 29}, -3}, {{29, 28}, -3},
        {{35, 34}, -4}, {{34, 33}, -4}, {{33, 32}, -4},
        {{32, 24}, -5}, {{24, 16}, -5}, {{16, 8}, -5}, {{8, 0}, -5},
        {{28, 35}, -6}};
    gmesh->NodeVec().Resize(co.size());
    for (size_t i = 0; i < co.size(); i++) {
        TPZManVector<REAL, 3> x = {co[i][0], co[i][1], 0.};
        gmesh->NodeVec()[i] = TPZGeoNode(i, x, *gmesh);
    }
    for (auto &t : tri) {
        TPZManVector<int64_t, 3> nodes = {t[0], t[1], t[2]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodes, 1, *gmesh);
    }
    for (auto &l : lines) {
        TPZManVector<int64_t, 2> nodes = {l.first[0], l.first[1]};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, l.second, *gmesh);
    }
    gmesh->BuildConnectivity();
    for (int d = 0; d < ref; d++) {
        const int64_t nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> sub;
        for (int64_t iel = 0; iel < nel; iel++)
            if (!gmesh->Element(iel)->HasSubElement()) gmesh->Element(iel)->Divide(sub);
    }
    return gmesh;
}

/// H1 mesh with memory; base fixed, lateral sides on rollers (BC type 3: val2 = constrained directions)
template <class TPlastic>
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int porder, const TPlastic &model, const Soil &s) {
    auto *mat = new TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem>(1, 1);
    TPlastic m(model);
    mat->SetPlasticityModel(m);
    mat->SetBodyForce({0., -s.gamma, 0.});
    mat->SetBigNumber(1.e12);
    auto *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetDimModel(2);
    cmesh->InsertMaterialObject(mat);
    TPZFMatrix<STATE> val1(2, 2, 0.);
    const int bcid[3] = {-1, -2, -5};
    const REAL fixy[3] = {1., 0., 0.};
    for (int i = 0; i < 3; i++) {
        TPZManVector<STATE, 2> val2 = {1., fixy[i]};
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcid[i], 3, val1, val2));
    }
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

#endif
