// FEM stability with seepage forces: gravity increase of Projects/SlopeMohrCoulomb/SlopeAnalysis.h (unchanged) on
// the stability mesh with the body force
//     b = lambda (gamma_ref g + f(x)),   g = (0, -1),
// and a stress-free effective boundary. With f = -grad u (excess pore pressure of SeepageFE.h) and
// gamma_ref = gamma' = gamma - gamma_w this is the loading of Ceron et al. (2025) and Gamma_FEM = lambda_crit
// (H_crit = lambda_crit H). With f = -grad p+ (p = u - gamma_w y) and gamma_ref = gamma_sat it is the form of
// Projects/SlopeDrawdown (b = lambda (gamma_sat g - grad p+)); both coincide where p >= 0. Dry slope: f = 0,
// gamma_ref = gamma.
// How lambda reaches the forcing function: TPZMatElastoPlastic2D::Contribute starts the local body force from the
// material body force m_force (which SlopeAnalysis sets to lambda (0, -gamma_ref, 0)) and lets the forcing function
// overwrite it, so the forcing function recovers lambda = -m_force[1] / gamma_ref (as in SlopeDrawdown).
#ifndef FEMSTABILITY_H
#define FEMSTABILITY_H

#include "../SlopeMohrCoulomb/SlopeAnalysis.h"
#include "../SlopeMohrCoulomb/SlopeModel.h"
#include "SeepageForceField.h"
#include "SlopeGeometry.h"

#include <chrono>
#include <memory>
#include <string>
#include <vector>

namespace slope {

/// f = -grad p+ of the total pore pressure p = u - gamma_w y (y up, crest y = 0; p+ = max(p, 0)): the seepage term of
/// SlopeDrawdown, to be used with gamma_ref = gamma_sat
inline ForceField TotalPressureForce(std::shared_ptr<const PoreField> pf, REAL gammaw) {
    return [pf, gammaw](const TPZVec<REAL> &x, REAL f[2]) {
        REAL u, g[2];
        if (!pf->Evaluate(x, u, g) || u - gammaw * x[1] <= 0.) {
            f[0] = f[1] = 0.;
            return;
        }
        f[0] = -g[0];
        f[1] = -(g[1] - gammaw);
    };
}

/// Body force b = lambda (gamma_ref g + f(x)) of the soil (material 1); lambda is the gravity factor that
/// SlopeAnalysis sets through the body force of the material (adapted from SetSeepageForce of SlopeDrawdown).
/// Must be called before SlopeAnalysis is constructed (it stores the reference body force).
template <class TPlastic>
void SetSeepageForce(TPZCompMesh *cmesh, const ForceField &f, REAL gammaRef) {
    auto *mat = dynamic_cast<TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem> *>(cmesh->FindMaterial(ESoil));
    if (!mat || gammaRef <= 0.) DebugStop();
    mat->SetBodyForce({0., -gammaRef, 0.});
    mat->SetForcingFunction(
        [mat, f, gammaRef](const TPZVec<REAL> &x, TPZVec<STATE> &F) {
            const REAL lambda = -mat->GetBodyForce()[1] / gammaRef;
            REAL fs[2];
            f(x, fs);
            F[0] = lambda * fs[0];
            F[1] = lambda * (fs[1] - gammaRef);
            F[2] = 0.;
        },
        0);
}

struct FSCycle {
    int cycle = 0;
    int64_t neq = 0;
    REAL gi = 0., srm = -1.; ///< gravity-increase factor; SRM factor (-1: not computed)
    REAL zone[4] = {0., 0., 0., 0.}; ///< xmin, ymin, xmax, ymax of the plastic zone at the GI collapse
    double seconds = 0.;     ///< wall time of the cycle
};

/// Bounding box of the elements with sqrt(J2(eps_p)) >= frac * max (the indicator of SlopeAnalysis::MarkPlasticZone)
/// in the last accepted state: the extent of the failure mechanism, to check that the box contains it
template <class TPlastic>
void PlasticZoneBox(TPZCompMesh *cmesh, REAL frac, REAL box[4]) {
    auto *mat = dynamic_cast<TPZMatElastoPlastic2D<TPlastic, TPZElastoPlasticMem> *>(cmesh->FindMaterial(ESoil));
    std::vector<std::pair<TPZCompEl *, REAL>> ind;
    REAL vmax = 0.;
    for (int64_t el = 0; el < cmesh->NElements(); el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel || cel->Material() != mat) continue;
        TPZManVector<int64_t> mem;
        cel->GetMemoryIndices(mem);
        REAL v = 0.;
        for (int64_t m : mem) {
            if (m < 0) continue;
            TPZTensor<REAL> ep = mat->MemItem(m).m_elastoplastic_state.m_eps_p;
            ep.XY() *= 0.5; ep.XZ() *= 0.5; ep.YZ() *= 0.5;
            v = std::max(v, std::sqrt(std::max<REAL>(ep.J2(), 0.)));
        }
        ind.push_back({cel, v});
        vmax = std::max(vmax, v);
    }
    box[0] = box[1] = 1.e300, box[2] = box[3] = -1.e300;
    for (auto &e : ind) {
        if (vmax <= 0. || e.second < frac * vmax) continue;
        TPZGeoEl *gel = e.first->Reference();
        for (int k = 0; k < gel->NCornerNodes(); k++) {
            box[0] = std::min(box[0], gel->NodePtr(k)->Coord(0)), box[2] = std::max(box[2], gel->NodePtr(k)->Coord(0));
            box[1] = std::min(box[1], gel->NodePtr(k)->Coord(1)), box[3] = std::max(box[3], gel->NodePtr(k)->Coord(1));
        }
    }
}

/// Gravity increase with nref refinement cycles of the plastic zone (adapted from FactorOfSafety of SlopeDrawdown):
/// each cycle marks the elements with sqrt(J2(eps_p)) >= 10% of the max at collapse. With srm = true the strength
/// reduction is also computed and its plastic zone is marked too (the procedure of SlopeDrawdown, whose mesh
/// sequence it reproduces). The geometric mesh is refined in place.
template <class TPlastic>
std::vector<FSCycle> GravityIncreaseFS(TPZGeoMesh *gmesh, const TPlastic &model, const Soil &s, const ForceField &f,
                                       REAL gammaRef, int nref, bool srm, const std::string &vtk = "") {
    TPZCompMesh *cmesh = CreateCMesh(gmesh, 2, model, s);
    SetSeepageForce<TPlastic>(cmesh, f, gammaRef);
    std::vector<FSCycle> out;
    {
        SlopeAnalysis<TPlastic> slope(cmesh);
        for (int k = 0;; k++) {
            const auto t0 = std::chrono::steady_clock::now();
            FSCycle c;
            c.cycle = k;
            c.neq = cmesh->NEquations();
            c.gi = slope.GravityIncrease();
            PlasticZoneBox<TPlastic>(cmesh, 0.1, c.zone);
            if (!vtk.empty()) slope.PostPlasticity(vtk + "_GI_ref" + std::to_string(k) + ".vtk");
            slope.MarkPlasticZone(0.1);
            if (srm) {
                c.srm = slope.StrengthReduction();
                slope.MarkPlasticZone(0.1);
            }
            c.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "[FS] cycle " << k << ": " << c.neq << " equations, lambda GI " << c.gi;
            if (srm) std::cout << ", FS SRM " << c.srm;
            std::cout << " (" << c.seconds << " s); GI plastic zone x in [" << c.zone[0] << ", " << c.zone[2] << "], y in ["
                      << c.zone[1] << ", " << c.zone[3] << "]\n";
            out.push_back(c);
            if (k >= nref) break;
            slope.Refine();
        }
    }
    delete cmesh;
    return out;
}

} // namespace slope

#endif
