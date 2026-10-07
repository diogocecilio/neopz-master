/**
 * @file EmbankmentConsolidation.h
 * @brief Sect. 6.6 of the article: embankment loading on a Modified Cam-Clay foundation (FLAC3D example).
 * A 1 m slice of the foundation is modelled with 20 x 10 x 1 Hex20-Hex8 u-p elements in plane strain, with
 * geostatic initial state, undrained loading and consolidation (Figs. 12 to 14, Table 8 and the embankment
 * column of Table 10).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylnsymmat.h"
#include "pzstepsolver.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Embankment on a Cam-Clay foundation (function aterro() of gen_data.py, aterro_elastic.py and the
 * embankment part of tangentes() of gen_data.py), three-dimensional model.
 *
 * Half of the problem: the slab \f$[0,20]\times[0,10]\times[0,1]\f$ m (x horizontal, y vertical, z the thickness),
 * the 1 m slice of the FLAC3D model, divided into 20 x 10 x 1 Hex20-Hex8 elements of 1 m (serendipity quadratic
 * displacement, trilinear pore pressure) with 3 x 3 x 3 Gauss points (Fig. 12). Plane strain is imposed with
 * \f$u_z=0\f$ on the faces z = 0 and z = 1 m.
 *
 * Material (Table 1): \f$M=0.888\f$, \f$\lambda=0.161\f$, \f$\kappa=0.062\f$, \f$v_\lambda=2.858\f$,
 * porous elasticity with constant Poisson ratio \f$\nu=0.3\f$ and uniform preconsolidation pressure
 * \f$p'_{c0}=160\f$ kPa. Unit weights \f$\gamma_{sat}=23\f$ and \f$\gamma_w=10\f$ kN/m\f$^3\f$, pore fluid with
 * \f$K_f=2\times10^5\f$ kPa and porosity \f$n=0.3\f$ (\f$\alpha_B=1\f$, \f$1/M_B=n/K_f\f$), mobility
 * \f$k=10^{-9}\f$ m\f$^2\f$/(kPa s); body force \f$(0,-\gamma_{sat},0)\f$ and fluid weight \f$(0,-\gamma_w,0)\f$.
 *
 * Initial state at the depth \f$d=10-y\f$: hydrostatic pore pressure \f$p_w=\gamma_w d\f$, effective stresses
 * \f$\sigma'_{yy}=-(\gamma_{sat}-\gamma_w)d\f$ and \f$\sigma'_{xx}=\sigma'_{zz}=0.7(-\gamma_{sat}d)+\gamma_w d\f$ (the
 * total horizontal stress is 0.7 times the total vertical stress), \f$p'_c=160\f$ kPa and the specific volume
 * \f$v_0=v_\lambda-\lambda\ln p'_{c0}+\kappa\ln(p'_{c0}/p'_0)\f$ at each integration point.
 *
 * Boundary conditions: \f$u_x=0\f$ at x = 0 (symmetry) and x = 20 m, \f$u_z=0\f$ at z = 0 and z = 1 m (plane
 * strain), fixed and impermeable base, drained top (\f$p_w=0\f$) and the strip load q = 50 kPa on
 * \f$0\le x\le 4\f$ m of the top, scaled by the load factor.
 *
 * Loading: ten undrained increments of the load factor (\f$\Delta t=0\f$), followed by 25 consolidation steps
 * with four steps per decade from \f$t=10^2\f$ to \f$10^8\f$ s; no predictor. Monitoring: settlements of the
 * top at x = 0, 2, 4 and 6 m (vertices of the face z = 0) and the pore pressures pp1 and pp2, averages of the
 * eight vertices of the elements centred at (0.5, 9.5, 0.5) and (1.5, 7.5, 0.5) m.
 *
 * The plane strain solution of the 20 x 10 Q8-Q4 model of the article v0.6 belongs to the Hex20-Hex8 space of
 * the slab (it does not depend on z and has \f$u_z=0\f$) and satisfies its discrete equations: with a
 * z-independent stress and \f$\sigma'_{xz}=\sigma'_{yz}=0\f$, the integral over the thickness of each Hex20 test
 * function is a combination of Q8 test functions of the same boundary conditions, and that of each Hex8 test
 * function is half of a Q4 one. The 3D model therefore reproduces the 2D solution; only the normalized residual
 * \f$\|R\|/\|f_{ext}\|\f$, which is computed with the nodal components of a different basis, changes, and with it
 * possibly the number of iterations of an increment and the converged solution at the level of the tolerance.
 *
 * Models (EVariant): the Cam-Clay foundation and the elastic one of aterro_elastic.py, i.e. the same model with
 * \f$p_c=10^7\f$ kPa in the state (no yielding) and \f$v_0\f$ computed with \f$p'_{c0}=160\f$ kPa, which isolates
 * the plastic share of the settlement. The comparison of the tangent operators (Table 10) repeats the Cam-Clay
 * analysis with each operator returned by TPZPlasticStepModifiedCamClay::SetTangentMode; its run with the
 * transpose \f$D^T\f$ also gives the undrained loading with the transposed tangent of the end of aterro().
 *
 * The class follows the structure of the NeoPZ examples: geometric mesh, computational meshes
 * (displacement, pore pressure and multiphysics), analysis with structural matrix and solver,
 * incremental solution and post-processing.
 */
class EmbankmentConsolidation {
public:
    /** @brief Material and boundary ids (the Python markers are 1 bottom, 2 right, 3 top, 4 left, 5 loaded top) */
    enum {
        EMatId = 1,   ///< clay layer (u-p material)
        EBottom = -1, ///< base y = 0: \f$u=0\f$, impermeable
        ERight = -2,  ///< x = 20 m: \f$u_x=0\f$
        ELeft = -4,   ///< x = 0 (symmetry): \f$u_x=0\f$
        ELoad = -5,   ///< loaded strip of the top, \f$0\le x\le 4\f$ m: traction (0, -q, 0) times the load factor
        EZ0 = -6,     ///< face z = 0: \f$u_z=0\f$ (plane strain)
        EZ1 = -7,     ///< face z = 1 m: \f$u_z=0\f$ (plane strain)
        EPTop = -13   ///< whole top y = 10 m (coincident faces): drained, \f$p_w=0\f$
    };

    /** @brief Models solved by Run */
    enum EVariant {
        ECamClay, ///< Modified Cam-Clay foundation, undrained loading and consolidation (aterro())
        EElastic  ///< no yielding (\f$p_c=10^7\f$ kPa), undrained loading and consolidation (aterro_elastic.py)
    };

    /** @brief Tangent operator returned by the stress update (Table 10) */
    typedef TPZPlasticStepModifiedCamClay::ETangentMode ETangentMode;

    /** @brief Parts of RunAll (bit mask) */
    enum EPart {
        EPartModel = 1,    ///< Cam-Clay model with the consistent tangent: Table 8, Figs. 13 and 14, CSV and VTK files
        EPartElastic = 2,  ///< elastic foundation (plastic share of the settlement)
        EPartTangents = 4, ///< comparison of the tangent operators (Table 10, embankment column)
        EPartAll = 7       ///< all of them
    };

    /** @name Loading (Sect. 6.6) */
    /** @{ */
    static constexpr int kNUndrained = 10;    ///< undrained increments of the load factor
    static constexpr int kNConsolidation = 25; ///< consolidation steps, \f$t=10^{2+j/4}\f$ s, j = 0..24
    static constexpr int kNStates = kNUndrained + kNConsolidation + 1; ///< monitored states (initial state included)
    /** @} */

    /** @brief Results of a run */
    struct TResult {
        ETangentMode fMode = TPZPlasticStepModifiedCamClay::EConsistentTangent; ///< tangent operator of the run
        /** @brief Monitored states: rows (t, lambda, s(x=0), s(x=2), s(x=4), s(x=6), pp1, pp2), the initial state first */
        std::vector<std::vector<REAL>> fHistory;
        size_t fNUndrained = 0;                                  ///< rows of the undrained stage (initial state included)
        std::vector<REAL> fLayerZ;                               ///< z of the layers of integration points
        std::vector<std::array<int, 3>> fTypesUndrained;         ///< elastic, subcritical, supercritical points per layer, end of the loading
        std::vector<std::array<int, 3>> fTypesFinal;             ///< elastic, subcritical, supercritical points per layer, t = 1e8 s
        std::array<int, 3> fTotalUndrained = {0, 0, 0};          ///< totals of fTypesUndrained
        std::array<int, 3> fTotalFinal = {0, 0, 0};              ///< totals of fTypesFinal
        std::vector<mcc::TAnalysis::TStepLog> fLogUndrained;     ///< convergence records of the undrained increments
        std::vector<mcc::TAnalysis::TStepLog> fLogConsolidation; ///< convergence records of the consolidation steps
        int64_t fNGlobalIterations = 0; ///< global iterations, failed attempts included (TPZPoroElastoPlasticUPAnalysis::NGlobalIterations)
        int64_t fNBisections = 0;       ///< failed (sub)increments (TPZPoroElastoPlasticUPAnalysis::NBisections)
        REAL fR0 = 0.;                ///< largest initial residual at the free displacement equations (nodal basis)
        REAL fBaseReaction0 = 0.;     ///< initial vertical reaction of the base (weight of the soil)
        REAL fReactionX0 = 0.;        ///< horizontal reaction at x = 0 (base nodes excluded), t = 1e8 s
        REAL fReactionX20 = 0.;       ///< horizontal reaction at x = 20 m (base nodes excluded), t = 1e8 s
        REAL fReactionBaseX = 0.;     ///< horizontal reaction of the base, t = 1e8 s
        REAL fReactionBaseY = 0.;     ///< vertical reaction of the base, t = 1e8 s
        /** @brief Out-of-plane reaction of the face z = 0 at t = 1e8 s: all its Lagrange nodes, those of the base line
         * included, so that it equals \f$-\int\sigma_{zz}\,dA\f$ over the face (total stress), the force that keeps the
         * plane strain */
        REAL fReactionZ0 = 0.;
        REAL fReactionZ1 = 0.;        ///< the same for the face z = 1 m (\f$+\int\sigma_{zz}\,dA\f$)
        /** @brief Check of fReactionZ0 at t = 1e8 s: \f$-\frac{1}{T}\int_V(\sigma'_{zz}-\alpha_B p_w)\,dV\f$ summed over the
         * integration points (T = 1 m, the fields do not depend on z) */
        REAL fSigmaZZForce = 0.;
        REAL fExcessUndrained = 0.;   ///< largest excess pore pressure at the end of the loading (Fig. 14a)
        std::array<REAL, 3> fExcessUndrainedX = {0., 0., 0.}; ///< vertex where it occurs
        REAL fExcess1e6 = 0.;         ///< largest excess pore pressure at t = 1e6 s (Fig. 14b)
        std::array<REAL, 3> fExcess1e6X = {0., 0., 0.};       ///< vertex where it occurs
        REAL fHeave = 0.;             ///< largest heave of the top at t = 1e8 s (negative settlement, Fig. 14d)
        REAL fHeaveStart = -1.;       ///< first x of the top with heave at t = 1e8 s (-1 if none)
        REAL fCvUndrained = 0.;       ///< consolidation coefficient of the skeleton in zone pp2, end of the loading
        REAL fCvFinal = 0.;           ///< consolidation coefficient of the skeleton in zone pp2, t = 1e8 s
        /** @brief Plane strain check at t = 1e8 s: largest \f$|u(x,y,0)-u(x,y,1)|\f$ over the nodes of the two faces */
        REAL fPlaneStrainU = 0.;
        /** @brief Plane strain check at t = 1e8 s: largest \f$|p(x,y,0)-p(x,y,1)|\f$ over the pore pressure nodes */
        REAL fPlaneStrainP = 0.;
        /** @brief Plane strain check at t = 1e8 s: largest \f$|u_z|\f$ at the mid-edge nodes z = 0.5 m (free u_z) */
        REAL fPlaneStrainUz = 0.;
        /** @brief Plane strain check at t = 1e8 s: largest difference of \f$u_x,u_y\f$ between a mid-edge node z = 0.5 m and the vertex below */
        REAL fPlaneStrainMid = 0.;
        int64_t fNEquations = 0;         ///< number of equations of the mesh
        int64_t fNFreeEquations = 0;     ///< equations of the linear systems (Dirichlet equations eliminated)
        int64_t fNDisplacementNodes = 0; ///< number of displacement nodes (vertices and mid-edge nodes)
        int64_t fNPressureNodes = 0;     ///< number of pore pressure nodes
        int64_t fProfile = 0;            ///< entries of the upper (= lower) skyline profile of the linear system
        double fLUFlops = 0.;            ///< multiply-adds of one LU decomposition of the skyline matrix (estimate)
        double fSeconds = 0.;            ///< run time
        double fOutputSeconds = 0.;      ///< part of the run time spent writing CSV and VTK files
        bool fConverged = true;          ///< false if an increment failed
    };

    /** @name Data of the example (Table 1 and Sect. 6.6) */
    /** @{ */
    REAL fM = 0.888;        ///< slope of the critical state line
    REAL fLambda = 0.161;   ///< slope of the normal compression line
    REAL fKappa = 0.062;    ///< slope of the unloading-reloading line
    REAL fVLambda = 2.858;  ///< specific volume of the normal compression line at p' = 1 kPa
    REAL fNu = 0.3;         ///< Poisson ratio of the porous elasticity, \f$G=3K(1-2\nu)/(2(1+\nu))\f$
    REAL fPc0 = 160.;       ///< uniform preconsolidation pressure (kPa)
    REAL fPcElastic = 1.e7; ///< preconsolidation pressure of the elastic variant (kPa)
    REAL fGammaSat = 23.;   ///< saturated unit weight (kN/m3)
    REAL fGammaW = 10.;     ///< unit weight of the water (kN/m3)
    REAL fKf = 2.e5;        ///< bulk modulus of the pore fluid (kPa)
    REAL fPorosity = 0.3;   ///< porosity, \f$1/M_B=n/K_f\f$
    REAL fMobility = 1.e-9; ///< mobility k (m2/(kPa s))
    REAL fQ = 50.;          ///< embankment load (kPa)
    REAL fWidth = 20.;      ///< width of the model (m)
    REAL fHeight = 10.;     ///< thickness of the clay layer (m)
    REAL fThickness = 1.;   ///< thickness of the slice (z direction, m)
    REAL fLoadWidth = 4.;   ///< loaded strip \f$0\le x\le 4\f$ m
    REAL fRatioTotal = 0.7; ///< ratio of the total horizontal to the total vertical geostatic stress
    int fNx = 20;           ///< elements along x
    int fNy = 10;           ///< elements along y
    /** @} */

    /** @name Options of RunAll (command line of main.cpp) */
    /** @{ */
    /**
     * @brief Write the VTK file series of every converged state (mcc::TVTKSeries) of the Cam-Clay and elastic runs,
     * in the directories vtk/embankment and vtk/embankment_elastic (the command line argument "novtk" disables it)
     */
    bool fWriteVTK = true;
    /** @brief Parts of RunAll (EPart bit mask; command line arguments "model", "elastic", "tangents") */
    int fParts = EPartAll;
    /** @brief Tangent operators compared in Table 10 (command line argument "modes=D,sym,cont,DT,fd") */
    std::vector<ETangentMode> fTangentModes = {
        TPZPlasticStepModifiedCamClay::EConsistentTangent, TPZPlasticStepModifiedCamClay::ESymmetricTangent,
        TPZPlasticStepModifiedCamClay::EContinuumTangent, TPZPlasticStepModifiedCamClay::ETransposedTangent,
        TPZPlasticStepModifiedCamClay::EFiniteDifferenceTangent};
    /** @brief Verbosity of the analysis (TPZPoroElastoPlasticUPAnalysis::SetVerbose; argument "verbose" sets 1) */
    int fVerbose = 0;
    /** @} */

    /**
     * @brief Effective geostatic stresses (tension positive) at the height y
     * @param y height (m); the depth is \f$d=10-y\f$
     * @param[out] sv vertical effective stress \f$-(\gamma_{sat}-\gamma_w)d\f$
     * @param[out] sh horizontal effective stress \f$0.7(-\gamma_{sat}d)+\gamma_w d\f$ (x and z)
     */
    void GeostaticStress(REAL y, REAL &sv, REAL &sh) const {
        const REAL d = fHeight - y;
        sv = -(fGammaSat - fGammaW) * d;
        sh = fRatioTotal * (-fGammaSat * d) + fGammaW * d;
    }

    /** @brief Hydrostatic pore pressure \f$\gamma_w(10-y)\f$ (water table at the surface) at the height y */
    REAL Hydrostatic(REAL y) const { return fGammaW * (fHeight - y); }

    /**
     * @brief Geometric mesh: 20 x 10 x 1 trilinear hexahedra (mcc::CreateSlabMesh) with the boundary
     * quadrilaterals EBottom, ERight, ELeft, EZ0, EZ1, EPTop (whole top) and ELoad (coincident with EPTop on the
     * faces of the loaded strip)
     * @return the geometric mesh
     */
    TPZGeoMesh *CreateGeoMesh();

    /**
     * @brief Displacement, pore pressure and multiphysics meshes, material, boundary conditions and initial
     * state (memory of the integration points and pore pressure field)
     * @param gmesh geometric mesh
     * @param variant model (elastic variant: pc = 1e7 kPa in the state)
     * @param mode tangent operator returned by the stress update
     * @param[out] mat the u-p material of the clay layer
     * @return the multiphysics mesh
     */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, EVariant variant, ETangentMode mode,
                                            mcc::TPoroMaterial *&mat);

    /**
     * @brief Largest component, in the nodal basis, of the initial residual \f$F_{int}-QP-f_b\f$ at the free
     * displacement equations (equilibrium of the geostatic state)
     * @param an analysis of the model, with the initial solution loaded
     * @param mat the u-p material
     * @param[out] basereaction initial vertical reaction of the base
     * @return the largest absolute value of the residual at the free displacement equations
     */
    REAL InitialResidual(mcc::TAnalysis &an, mcc::TPoroMaterial *mat, REAL &basereaction);

    /**
     * @brief Mean over the integration points of an element of the consolidation coefficient of the skeleton,
     * \f$c_v=k(K+4G/3)\f$ with \f$K=v_0p'/\kappa\f$ and \f$G=3K(1-2\nu)/(2(1+\nu))\f$ (Sect. 6.6)
     * @param mat the u-p material
     * @param mphys multiphysics mesh
     * @param gel volume element
     * @return mean of \f$c_v\f$ over the integration points of gel (m2/s)
     */
    REAL ConsolidationCoefficient(mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys, TPZGeoEl *gel) const;

    /**
     * @brief Displacement at the nodes of the serendipity (Lagrange) basis: the vertices and the mid-edge nodes
     * of the Hex20 elements; the value at a mid-edge node is the coefficient of the hierarchical edge function plus
     * the mean of the vertex values (TPZPoroElastoPlasticUPAnalysis::SetNodalResidualNorm)
     * @param an analysis (solution loaded in the mesh)
     * @param mphys multiphysics mesh
     * @return rows (x, y, z, u_x, u_y, u_z), vertices first
     */
    std::vector<std::array<REAL, 6>> LagrangeDisplacements(mcc::TAnalysis &an, TPZMultiphysicsCompMesh *mphys) const;

    /**
     * @brief Plane strain check: differences of the solution between the faces z = 0 and z = 1 m and the
     * out-of-plane displacement of the mid-edge nodes z = 0.5 m (TResult::fPlaneStrainU, fPlaneStrainP,
     * fPlaneStrainUz, fPlaneStrainMid)
     */
    void PlaneStrainCheck(mcc::TAnalysis &an, TPZMultiphysicsCompMesh *mphys, TResult &res) const;

    /**
     * @brief Types of response of the integration points per layer z = const of points
     * @param mat the u-p material
     * @param mphys multiphysics mesh
     * @param[out] z coordinates of the layers (ascending)
     * @param[out] types elastic, subcritical and supercritical points of each layer
     * @param[out] total the same over the whole mesh
     */
    void CountTypesByLayer(mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys, std::vector<REAL> &z,
                           std::vector<std::array<int, 3>> &types, std::array<int, 3> &total) const;

    /**
     * @brief Profile of the skyline matrix of the last linear system (the Dirichlet equations eliminated):
     * writes embankment_matrix_profile.csv (equation, height of its column, 1 for a pore pressure equation)
     * and returns the size of the profile and the estimated multiply-adds of one LU decomposition
     */
    void MatrixProfile(mcc::TAnalysis &an, TResult &res, bool write) const;

    /**
     * @brief Writes the fields of a state (Fig. 14): nodal VTK (displacement and pore pressure), VTK of the
     * integration points, CSV of the vertex values (with the excess pore pressure) and of the integration points
     * @param an analysis (solution of the state loaded in the mesh)
     * @param mat the u-p material
     * @param mphys multiphysics mesh
     * @param tag name of the state (undrained, t1e6, t1e8)
     * @param step index of the nodal VTK file
     * @param[out] excess largest excess pore pressure \f$p-\gamma_w(10-y)\f$ at the vertices
     * @param[out] xexcess vertex where it occurs
     */
    void WriteState(mcc::TAnalysis &an, mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys, const std::string &tag,
                    int step, REAL &excess, std::array<REAL, 3> &xexcess);

    /**
     * @brief Solves a model (analysis, structural matrix, solver, increments) and post-processes it
     * @param variant ECamClay or EElastic: undrained loading and consolidation
     * @param mode tangent operator returned by the stress update (SetTangentMode)
     * @param output write the files of the run: CSV of the mesh, histories, convergence records and matrix
     * profile, fields of Fig. 14 (ECamClay) and the VTK series (with fWriteVTK); without output only the
     * results are returned (comparison of the tangent operators)
     * @return monitored history, convergence records and post-processed quantities of the run
     *
     * With output and fWriteVTK, the run writes the VTK file series of the 36 converged states
     * (mcc::TVTKSeries) in vtk/embankment or vtk/embankment_elastic: the series time is the index of the
     * state (0 initial state, 1 to 10 undrained increments, 11 to 35 consolidation steps) and the file
     * \<name\>_states.csv gives the time t and the load factor of each index.
     */
    TResult Run(EVariant variant, ETangentMode mode, bool output);

    /** @brief Runs the parts selected by fParts and prints the comparison with the article v0.6 (2D model) */
    void RunAll();

    /** @name Reference values: the 20 x 10 Q8-Q4 plane strain model of the article v0.6 (Python data_aterro.pkl,
     * aterro_elastic.pkl and data_tangentes.pkl, reproduced to round-off by the 2D NeoPZ model) */
    /** @{ */
    /** @brief History of the Cam-Clay model: rows (t, s0, s2, s4, s6, pp1, pp2); the first 11 rows are the undrained loading */
    static constexpr REAL kRefHistory[kNStates][7] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 24.999999999999993},
        {0.0, 0.02269168900824459, 0.022452700209367543, 0.009543648886918481, -0.004050064901332618, 8.41193222759174, 27.751256612641363},
        {0.0, 0.04062291800227912, 0.040227881037782263, 0.018107675398208657, -0.007955844264457473, 11.310309419058186, 30.722881909118733},
        {0.0, 0.056580232621990544, 0.05606140155254447, 0.02571668756381279, -0.011828997808386624, 14.087809910177807, 33.76606133634975},
        {0.0, 0.07152791128242972, 0.0709166328502578, 0.03264574880428218, -0.01568726340085298, 16.825610713869384, 36.84544094554266},
        {0.0, 0.0858584737509302, 0.0851875025653235, 0.03909184206965954, -0.019532765441910224, 19.54612445474649, 39.9495680457187},
        {0.0, 0.09976517863923912, 0.09906760546764498, 0.04518243939605793, -0.023365062272500427, 22.257285952143555, 43.073476084509025},
        {0.0, 0.11335705311799237, 0.11266545911992026, 0.051001864831312076, -0.027183456285763494, 24.96247992747473, 46.21447879751899},
        {0.0, 0.12670189098560644, 0.1260478668123638, 0.05660820619508017, -0.030987384842819213, 27.66340246474848, 49.370855099357485},
        {0.0, 0.1398449969834096, 0.13925889486532772, 0.06204313035989801, -0.034776436368980035, 30.36104326862577, 52.54135174268359},
        {0.0, 0.1528183661090258, 0.15232918484645108, 0.06733763940150474, -0.038550304250082205, 33.05606629229955, 55.72496757943132},
        {100.0, 0.1528258631511103, 0.15233730370983228, 0.06734250649593855, -0.03854946707846289, 33.05301852825992, 55.72657879198736},
        {177.82794100389228, 0.15283169659784718, 0.15234362088668668, 0.06734629359444932, -0.03854881567121559, 33.05064646809869, 55.72783263987889},
        {316.22776601683796, 0.1528420669613567, 0.15235485099183085, 0.06735302583713025, -0.038547658012545394, 33.04643015702929, 55.73006277399918},
        {562.341325190349, 0.152860497981991, 0.15237480954283217, 0.06736498998748398, -0.03854560187488167, 33.03893574560506, 55.73402872627178},
        {1000.0, 0.15289324152568923, 0.15241026502431226, 0.06738624218691215, -0.03854195303101513, 33.02562322575098, 55.74108364680009},
        {1778.2794100389228, 0.15295136768764214, 0.15247319985773777, 0.06742396043069379, -0.038535488263957154, 33.001993216295, 55.7536349889873},
        {3162.2776601683795, 0.153054415903792, 0.15258475512102124, 0.0674908017656299, -0.03852406670390687, 32.960107798099074, 55.77597037144001},
        {5623.413251903491, 0.1532366808156036, 0.15278200780263054, 0.06760894142583286, -0.0385039866358712, 32.88604194718404, 55.81572464626204},
        {10000.0, 0.15355777751001493, 0.1531293106752484, 0.06781680255919184, -0.03846898015498537, 32.75560187388602, 55.88646042492434},
        {17782.794100389227, 0.15411968039922855, 0.15373639075382312, 0.06817972636211625, -0.038408804466338244, 32.52738922191555, 56.012022819026434},
        {31622.776601683792, 0.15509233930610447, 0.1547848379721564, 0.0688054766725757, -0.03830767051109124, 32.132121542744045, 56.23295706221749},
        {56234.13251903491, 0.15674774443902442, 0.15656060502022565, 0.06986347208536826, -0.0381433403090453, 31.457068656485205, 56.611908185910906},
        {100000.0, 0.159494295525961, 0.15947689572483828, 0.07160146088058034, -0.03788807124579101, 30.324562976877065, 57.22154719747242},
        {177827.94100389228, 0.16387912313888672, 0.16403779946209635, 0.07434279988308284, -0.03751036779338725, 28.46989262230556, 58.06495880201582},
        {316227.7660168379, 0.17047618866509437, 0.17065075194826873, 0.07843051026334981, -0.0369694206305187, 25.5770507609839, 58.84037885739023},
        {562341.3251903491, 0.179593976131448, 0.1793238223437227, 0.08408840384570536, -0.036181358006602216, 21.56655834408381, 58.578672270512484},
        {1000000.0, 0.19103857713061387, 0.18969463555775956, 0.09134379145407345, -0.03486512925941795, 17.046735853947, 55.81543829962162},
        {1778279.410038923, 0.20415195109638407, 0.20138450293189275, 0.10022768256806915, -0.03221960649696179, 12.96662898847607, 50.08948865392948},
        {3162277.6601683795, 0.2180453624221057, 0.2140317130442482, 0.11081676223648514, -0.02697945575085234, 9.857222302392085, 42.95014406694117},
        {5623413.251903491, 0.23204039140085705, 0.2272295972076745, 0.12289145166891476, -0.01831570060437987, 7.781801907395407, 36.50494962271083},
        {10000000.0, 0.2456883990032203, 0.24044304749715537, 0.1356657144609411, -0.006882144580835944, 6.512422799628082, 31.71724870470962},
        {17782794.100389227, 0.25806699274622275, 0.25256459180952867, 0.1476356459563155, 0.005069201750944501, 5.7540342044903605, 28.48354584571232},
        {31622776.60168379, 0.26751332048374915, 0.2618400002259704, 0.15683077281665525, 0.014702579410799855, 5.3162757163884145, 26.4907386279993},
        {56234132.51903491, 0.2729495343591471, 0.26717956064769127, 0.16212317606615978, 0.020375882183194312, 5.100844120855548, 25.480458732506797},
        {100000000.0, 0.27511301005450095, 0.2693052918781524, 0.1642308621899741, 0.02266768732429943, 5.022406847674765, 25.10745808288039}
    };
    /** @brief History of the elastic model: rows (t, s0, s2, s4, s6) */
    static constexpr REAL kRefElastic[kNStates][5] = {
        {0.0, -0.0, -0.0, -0.0, -0.0},
        {0.0, 0.02269168900824459, 0.022452700209367543, 0.009543648886918481, -0.004050064901332618},
        {0.0, 0.04062291800227912, 0.040227881037782263, 0.018107675398208657, -0.007955844264457473},
        {0.0, 0.056580232621990544, 0.05606140155254447, 0.02571668756381279, -0.011828997808386624},
        {0.0, 0.07152791128242972, 0.0709166328502578, 0.03264574880428218, -0.01568726340085298},
        {0.0, 0.0858584737509302, 0.0851875025653235, 0.03909184206965954, -0.019532765441910224},
        {0.0, 0.09976517863923912, 0.09906760546764498, 0.04518243939605793, -0.023365062272500427},
        {0.0, 0.11335705311799237, 0.11266545911992026, 0.051001864831312076, -0.027183456285763494},
        {0.0, 0.12670189098560644, 0.1260478668123638, 0.05660820619508017, -0.030987384842819213},
        {0.0, 0.1398449969834096, 0.13925889486532772, 0.06204313035989801, -0.034776436368980035},
        {0.0, 0.1528183661090258, 0.15232918484645108, 0.06733763940150474, -0.038550304250082205},
        {100.0, 0.1528258631511103, 0.15233730370983228, 0.06734250649593855, -0.03854946707846289},
        {177.82794100389228, 0.15283169659784718, 0.15234362088668668, 0.06734629359444932, -0.03854881567121559},
        {316.22776601683796, 0.1528420669613567, 0.15235485099183085, 0.06735302583713025, -0.038547658012545394},
        {562.341325190349, 0.152860497981991, 0.15237480954283217, 0.06736498998748398, -0.03854560187488167},
        {1000.0, 0.15289324152568923, 0.15241026502431226, 0.06738624218691215, -0.03854195303101513},
        {1778.2794100389228, 0.15295136768764214, 0.15247319985773777, 0.06742396043069379, -0.038535488263957154},
        {3162.2776601683795, 0.153054415903792, 0.15258475512102124, 0.0674908017656299, -0.03852406670390687},
        {5623.413251903491, 0.1532366808156036, 0.15278200780263054, 0.06760894142583286, -0.0385039866358712},
        {10000.0, 0.15355777751001493, 0.1531293106752484, 0.06781680255919184, -0.03846898015498537},
        {17782.794100389227, 0.15411968039922855, 0.15373639075382312, 0.06817972636211625, -0.038408804466338244},
        {31622.776601683792, 0.15509233930610447, 0.1547848379721564, 0.0688054766725757, -0.03830767051109124},
        {56234.13251903491, 0.15674769618461798, 0.15656056161817128, 0.06986344285783375, -0.03814334895095315},
        {100000.0, 0.15949388906424558, 0.15947653311828267, 0.07160122248232388, -0.03788813553832429},
        {177827.94100389228, 0.16387784872836553, 0.16403667466285038, 0.07434208197643945, -0.037510537627874034},
        {316227.7660168379, 0.1704728026598099, 0.1706478293805849, 0.07842874978553538, -0.03696974382847192},
        {562341.3251903491, 0.17958307171750204, 0.17931493182415428, 0.08408390144201597, -0.036181506620890175},
        {1000000.0, 0.19100919889844847, 0.18967120505857876, 0.09133291149737086, -0.034864312041086},
        {1778279.410038923, 0.20408973327040825, 0.20133490818834898, 0.10020498456961283, -0.032216593929124536},
        {3162277.6601683795, 0.21791339463850273, 0.2139251453899748, 0.11076616204243846, -0.026972116091875137},
        {5623413.251903491, 0.2316628496554387, 0.22691839802017166, 0.12273159944310374, -0.01830421498481221},
        {10000000.0, 0.24453322653526438, 0.23946937428198742, 0.13511433884996854, -0.006920341046461993},
        {17782794.100389227, 0.25528163443350693, 0.2501653843027491, 0.14613990047726766, 0.004717643743686508},
        {31622776.60168379, 0.2625776297875839, 0.2575141022805442, 0.15392356970944976, 0.013606520775622105},
        {56234132.51903491, 0.2662494490366335, 0.26124017529903454, 0.1579381879466301, 0.018410571640896795},
        {100000000.0, 0.2675236606783824, 0.26253965088348463, 0.15935397779352253, 0.02015696404117732}
    };
    /** @brief Evaluations of the residual per increment of the undrained loading (2D, all the tangent operators) */
    static constexpr int kRefItsUndrained[kNUndrained] = {5, 5, 5, 4, 4, 4, 4, 4, 4, 4};
    /**
     * @brief Evaluations of the residual per consolidation step in 2D with the operators D, sym, cont and DT
     * (rows in this order; data_tangentes.pkl; fd was not run on the embankment in v0.6)
     */
    static constexpr int kRefItsConsolidation[4][kNConsolidation] = {
        {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 4},
        {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 5, 5},
        {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8, 7, 6},
        {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 5}};
    /** @brief Total global iterations (itcount) of the 2D runs with D, sym, cont and DT (data_tangentes.pkl) */
    static constexpr int kRefItCount[4] = {132, 139, 153, 140};
    /** @brief Computing time (s) of the 2D Python runs with D, sym, cont and DT (data_tangentes.pkl) */
    static constexpr REAL kRefTime[4] = {20.5, 23.5, 26.4, 23.6};
    /** @} */

private:
    /** @brief Row of the reference tables kRefItsConsolidation, kRefItCount and kRefTime of a tangent mode (-1: none) */
    static int RefRow(ETangentMode mode) {
        switch (mode) {
        case TPZPlasticStepModifiedCamClay::EConsistentTangent:
            return 0;
        case TPZPlasticStepModifiedCamClay::ESymmetricTangent:
            return 1;
        case TPZPlasticStepModifiedCamClay::EContinuumTangent:
            return 2;
        case TPZPlasticStepModifiedCamClay::ETransposedTangent:
            return 3;
        default:
            return -1;
        }
    }

    /** @brief A line of embankment_summary.csv: quantity, value of this run, value of the 2D model of v0.6, unit */
    struct TSummaryItem {
        std::string fKey;  ///< name of the quantity
        REAL fValue;       ///< value of this run
        REAL fRef2D;       ///< value of the 2D model of v0.6 (NaN if none)
        std::string fUnit; ///< unit
    };
    /** @brief Lines of embankment_summary.csv collected by RunAll */
    std::vector<TSummaryItem> fSummary;
    /** @brief Adds a line to the summary (ref2D = NaN when the 2D model has no such value) */
    void Summary(const std::string &key, REAL value, REAL ref2d, const std::string &unit) {
        fSummary.push_back({key, value, ref2d, unit});
    }
    /** @brief Writes embankment_summary.csv */
    void WriteSummary(const std::string &file) const;
};

inline TPZGeoMesh *EmbankmentConsolidation::CreateGeoMesh() {
    const REAL xload = fLoadWidth, W = fWidth, H = fHeight, T = fThickness;
    return mcc::CreateSlabMesh(0., 0., fWidth, fHeight, fThickness, fNx, fNy, EMatId,
                               [=](const std::array<std::array<REAL, 3>, 4> &X) {
                                   if (mcc::FaceOnPlane(X, 1, 0.)) return std::vector<int>{EBottom};
                                   if (mcc::FaceOnPlane(X, 0, W)) return std::vector<int>{ERight};
                                   if (mcc::FaceOnPlane(X, 0, 0.)) return std::vector<int>{ELeft};
                                   if (mcc::FaceOnPlane(X, 2, 0.)) return std::vector<int>{EZ0};
                                   if (mcc::FaceOnPlane(X, 2, T)) return std::vector<int>{EZ1};
                                   if (mcc::FaceOnPlane(X, 1, H)) {
                                       // drained top; strip load on the faces whose centre has x <= 4 m
                                       const REAL xc = 0.25 * (X[0][0] + X[1][0] + X[2][0] + X[3][0]);
                                       if (xc <= xload + 1e-9) return std::vector<int>{EPTop, ELoad};
                                       return std::vector<int>{EPTop};
                                   }
                                   return std::vector<int>{};
                               });
}

inline TPZMultiphysicsCompMesh *EmbankmentConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, EVariant variant,
                                                                       ETangentMode mode, mcc::TPoroMaterial *&mat) {
    const int dim = 3;
    const std::set<int> bcids = {EBottom, ERight, ELeft, ELoad, EZ0, EZ1, EPTop};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, dim, EMatId, bcids); // serendipity Hex20
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, dim, EMatId, bcids);     // trilinear Hex8

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(dim);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EThreeDimensional);
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPorousElasticity();
    model.SetPoissonRatio(fNu);
    model.SetTangentMode(mode);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., fPorosity / fKf);
    mat->SetPermeability(fMobility);
    TPZManVector<REAL, 3> b(dim, 0.), rhowg(dim, 0.);
    b[1] = -fGammaSat;
    rhowg[1] = -fGammaW;
    mat->SetBodyForce(b);
    mat->SetFluidWeight(rhowg);
    mat->SetIntegrationOrder(4); // 3 x 3 x 3 Gauss points
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<9, STATE> val1(dim, dim, 0.);
    TPZManVector<STATE, 3> val2(dim, 0.), zero(1, 0.);
    // fixed base
    mphys->InsertMaterialObject(mat->CreateBC(mat, EBottom, B::EDirichletU, val1, val2));
    // u_x = 0 at x = 0 (symmetry) and x = 20 m; u_z = 0 at z = 0 and z = 1 m (plane strain)
    auto directional = [&](int id, int component) {
        val1.Zero();
        val1(component, component) = 1.;
        mphys->InsertMaterialObject(mat->CreateBC(mat, id, B::EDirichletUDirectional, val1, val2));
    };
    directional(ELeft, 0);
    directional(ERight, 0);
    directional(EZ0, 2);
    directional(EZ1, 2);
    // strip load, scaled by the load factor
    val1.Zero();
    TPZManVector<STATE, 3> load(dim, 0.);
    load[1] = -fQ;
    mphys->InsertMaterialObject(mat->CreateBC(mat, ELoad, B::ENeumannU, val1, load));
    // drained top
    mphys->InsertMaterialObject(mat->CreateBC(mat, EPTop, B::EDirichletP, val1, zero));
    mcc::BuildMultiphysics(mphys, cmeshU, cmeshP, EMatId, bcids);

    // geostatic state of each integration point (v0 always from pc0 = 160 kPa, as in aterro_elastic.py)
    const REAL pcstate = (variant == EElastic) ? fPcElastic : fPc0;
    mat->InitializeMemory(mphys, [&](const TPZVec<REAL> &x, TPZElastoPlasticMem &mem) {
        REAL sv, sh;
        GeostaticStress(x[1], sv, sh);
        const REAL p0 = -(2. * sh + sv) / 3.;
        const REAL v0 = mcc::SpecificVolumeNCL(fVLambda, fLambda, fKappa, fPc0, p0);
        mem.m_sigma.Zero();
        mem.m_sigma.XX() = sh;
        mem.m_sigma.YY() = sv;
        mem.m_sigma.ZZ() = sh;
        mem.m_elastoplastic_state.m_hardening = pcstate;
        mem.m_elastoplastic_state.fmatprop.Resize(1, v0);
        mem.m_elastoplastic_state.fmatprop[0] = v0;
        mem.m_elastoplastic_state.fpressure = Hydrostatic(x[1]);
    });
    // hydrostatic pore pressure at the pressure nodes
    mcc::SetInitialPressure(mphys, [this](const TPZVec<REAL> &x) { return Hydrostatic(x[1]); });
    return mphys;
}

inline REAL EmbankmentConsolidation::InitialResidual(mcc::TAnalysis &an, mcc::TPoroMaterial *mat,
                                                     REAL &basereaction) {
    // geostatic state without the embankment load
    mat->SetLoadFactor(0.);
    mat->SetTimeStep(0.);
    // flags of the constrained and pressure equations and edges of the nodal basis (built on demand by the analysis)
    an.IdentifyEquations();
    basereaction = an.Reaction({EBottom}, 1);
    TPZFMatrix<STATE> Ru;
    an.ComputeReactionVector(Ru);
    an.ToNodalBasis(Ru);
    const std::vector<bool> &constrained = an.ConstrainedEquations();
    const std::vector<bool> &pressure = an.PressureEquations();
    REAL rmax = 0.;
    for (int64_t eq = 0; eq < Ru.Rows(); ++eq)
        if (!constrained[eq] && !pressure[eq]) rmax = std::max(rmax, std::fabs(Ru.GetVal(eq, 0)));
    return rmax;
}

inline REAL EmbankmentConsolidation::ConsolidationCoefficient(mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys,
                                                              TPZGeoEl *gel) const {
    const REAL gfac = 3. * (1. - 2. * fNu) / (2. * (1. + fNu)); // G / K
    REAL sum = 0.;
    int n = 0;
    for (auto &g : mcc::GaussPoints(mat, mphys)) {
        if (mphys->Element(g.fElement)->Reference() != gel) continue;
        const REAL K = g.fV0 * mcc::MeanEffectiveStress(g.fSigma) / fKappa;
        sum += fMobility * K * (1. + 4. / 3. * gfac);
        n++;
    }
    return n ? sum / n : 0.;
}

inline std::vector<std::array<REAL, 6>> EmbankmentConsolidation::LagrangeDisplacements(mcc::TAnalysis &an,
                                                                                    TPZMultiphysicsCompMesh *mphys) const {
    std::vector<std::array<REAL, 6>> out;
    TPZGeoMesh *gmesh = mphys->Reference();
    TPZCompMesh *umesh = mphys->MeshVector()[0];
    const TPZFMatrix<STATE> &sol = mphys->Solution();
    std::set<int64_t> vertices;
    std::set<std::pair<int64_t, int64_t>> edges;
    auto vertex = [&](int64_t n) {
        std::array<REAL, 6> row;
        for (int k = 0; k < 3; ++k) row[k] = gmesh->NodeVec()[n].Coord(k);
        for (int d = 0; d < 3; ++d) row[3 + d] = an.NodalValue(n, 0, d);
        return row;
    };
    std::vector<std::array<REAL, 6>> mids;
    for (int64_t iel = 0; iel < umesh->NElements(); ++iel) {
        auto *intel = dynamic_cast<TPZInterpolatedElement *>(umesh->Element(iel));
        if (!intel || !intel->Reference() || intel->Reference()->Dimension() != 3) continue;
        TPZGeoEl *gel = intel->Reference();
        for (int i = 0; i < gel->NCornerNodes(); ++i)
            if (vertices.insert(gel->NodeIndex(i)).second) out.push_back(vertex(gel->NodeIndex(i)));
        for (int side = gel->NCornerNodes(); side < gel->NSides(); ++side) {
            if (gel->SideDimension(side) != 1) continue;
            const int64_t a = gel->SideNodeIndex(side, 0), b = gel->SideNodeIndex(side, 1);
            if (!edges.insert({std::min(a, b), std::max(a, b)}).second) continue;
            // the connect indices of the displacement mesh are those of the multiphysics mesh (first space)
            const int64_t ic = intel->ConnectIndex(intel->MidSideConnectLocId(side));
            const int64_t seq = mphys->ConnectVec()[ic].SequenceNumber();
            const int64_t pos = mphys->Block().Position(seq);
            const std::array<REAL, 6> va = vertex(a), vb = vertex(b);
            std::array<REAL, 6> row;
            for (int k = 0; k < 3; ++k) row[k] = 0.5 * (va[k] + vb[k]);
            for (int d = 0; d < 3; ++d) row[3 + d] = sol.GetVal(pos + d, 0) + 0.5 * (va[3 + d] + vb[3 + d]);
            mids.push_back(row);
        }
    }
    out.insert(out.end(), mids.begin(), mids.end());
    return out;
}

inline void EmbankmentConsolidation::PlaneStrainCheck(mcc::TAnalysis &an, TPZMultiphysicsCompMesh *mphys,
                                                      TResult &res) const {
    auto nodes = LagrangeDisplacements(an, mphys);
    auto key = [](REAL x, REAL y) { return std::make_pair(std::llround(x * 1e6), std::llround(y * 1e6)); };
    std::map<std::pair<int64_t, int64_t>, std::array<REAL, 6>> front, back;
    std::vector<std::array<REAL, 6>> middle;
    for (auto &n : nodes) {
        if (std::fabs(n[2]) < 1e-9) front[key(n[0], n[1])] = n;
        else if (std::fabs(n[2] - fThickness) < 1e-9) back[key(n[0], n[1])] = n;
        else middle.push_back(n);
    }
    res.fPlaneStrainU = res.fPlaneStrainUz = res.fPlaneStrainMid = res.fPlaneStrainP = 0.;
    for (auto &f : front) {
        auto it = back.find(f.first);
        if (it == back.end()) DebugStop();
        for (int d = 0; d < 3; ++d) res.fPlaneStrainU = std::max(res.fPlaneStrainU, std::fabs(f.second[3 + d] - it->second[3 + d]));
    }
    for (auto &m : middle) {
        res.fPlaneStrainUz = std::max(res.fPlaneStrainUz, std::fabs(m[5]));
        auto it = front.find(key(m[0], m[1]));
        if (it == front.end()) DebugStop();
        for (int d = 0; d < 2; ++d) res.fPlaneStrainMid = std::max(res.fPlaneStrainMid, std::fabs(m[3 + d] - it->second[3 + d]));
    }
    // pore pressure at the vertices of the two faces
    TPZGeoMesh *gmesh = mphys->Reference();
    std::map<std::pair<int64_t, int64_t>, REAL> pfront;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n)
        if (std::fabs(gmesh->NodeVec()[n].Coord(2)) < 1e-9)
            pfront[key(gmesh->NodeVec()[n].Coord(0), gmesh->NodeVec()[n].Coord(1))] = an.NodalValue(n, 1, 0);
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        if (std::fabs(gmesh->NodeVec()[n].Coord(2) - fThickness) > 1e-9) continue;
        auto it = pfront.find(key(gmesh->NodeVec()[n].Coord(0), gmesh->NodeVec()[n].Coord(1)));
        if (it == pfront.end()) DebugStop();
        res.fPlaneStrainP = std::max(res.fPlaneStrainP, std::fabs(an.NodalValue(n, 1, 0) - it->second));
    }
}

inline void EmbankmentConsolidation::CountTypesByLayer(mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys,
                                                       std::vector<REAL> &z, std::vector<std::array<int, 3>> &types,
                                                       std::array<int, 3> &total) const {
    std::map<int64_t, std::array<int, 3>> layers;
    total = {0, 0, 0};
    for (auto &g : mcc::GaussPoints(mat, mphys)) {
        auto &l = layers[std::llround(g.fX[2] * 1e9)];
        if (g.fType >= 0 && g.fType < 3) {
            l[g.fType]++;
            total[g.fType]++;
        }
    }
    z.clear();
    types.clear();
    for (auto &l : layers) {
        z.push_back(l.first * 1e-9);
        types.push_back(l.second);
    }
}

inline void EmbankmentConsolidation::MatrixProfile(mcc::TAnalysis &an, TResult &res, bool write) const {
    auto *sky = dynamic_cast<TPZSkylNSymMatrix<STATE> *>(an.MatrixSolver<STATE>().Matrix().operator->());
    if (!sky) return;
    const int64_t n = sky->Rows();
    res.fNFreeEquations = n;
    // full equations of the rows of the filtered system (the active equations in ascending order)
    const std::vector<bool> &constrained = an.ConstrainedEquations();
    const std::vector<bool> &pressure = an.PressureEquations();
    std::vector<int> ispressure;
    for (size_t eq = 0; eq < constrained.size(); ++eq)
        if (!constrained[eq]) ispressure.push_back(pressure[eq] ? 1 : 0);
    if ((int64_t)ispressure.size() != n) ispressure.assign(n, -1);
    std::vector<int64_t> top(n);
    res.fProfile = 0;
    for (int64_t j = 0; j < n; ++j) {
        const int64_t h = sky->SkyHeight(j);
        top[j] = j - h;
        res.fProfile += h;
    }
    // multiply-adds of the skyline LU: for each entry (i, j) of the profile, a dot product of length
    // i - max(top_i, top_j); twice for the non-symmetric factorization (L and U)
    double flops = 0.;
    for (int64_t j = 0; j < n; ++j)
        for (int64_t i = top[j]; i < j; ++i) flops += double(i - std::max(top[i], top[j]));
    res.fLUFlops = 2. * flops;
    if (!write) return;
    std::ofstream out("embankment_matrix_profile.csv");
    out << "equation,height,pressure\n";
    for (int64_t j = 0; j < n; ++j) out << j << "," << j - top[j] << "," << ispressure[j] << "\n";
}

inline void EmbankmentConsolidation::WriteState(mcc::TAnalysis &an, mcc::TPoroMaterial *mat,
                                                TPZMultiphysicsCompMesh *mphys, const std::string &tag, int step,
                                                REAL &excess, std::array<REAL, 3> &xexcess) {
    // native NeoPZ post-processing of the nodal fields and the integration points
    mcc::WriteNodalVTK(an, 3, "embankment_nodal.vtk", step);
    mcc::WriteGaussPointsVTK(mat, mphys, "embankment_gauss_" + tag + ".vtk");
    // vertex values: displacement, pore pressure and excess pore pressure (Fig. 14a, b and d use the face z = 0)
    TPZGeoMesh *gmesh = mphys->Reference();
    std::vector<std::vector<REAL>> rows;
    excess = -1e300;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        const REAL p = an.NodalValue(n, 1, 0);
        const REAL pex = p - Hydrostatic(x[1]);
        rows.push_back({x[0], x[1], x[2], an.NodalValue(n, 0, 0), an.NodalValue(n, 0, 1), an.NodalValue(n, 0, 2), p, pex});
        if (pex > excess + 1e-12) { // first vertex in the node order (face z = 0 first)
            excess = pex;
            xexcess = {x[0], x[1], x[2]};
        }
    }
    mcc::WriteCSV("embankment_nodal_" + tag + ".csv", {"x", "y", "z", "ux", "uy", "uz", "p", "p_excess"}, rows);
    // integration points (Fig. 14c)
    rows.clear();
    for (auto &g : mcc::GaussPoints(mat, mphys))
        rows.push_back({g.fX[0], g.fX[1], g.fX[2], mcc::MeanEffectiveStress(g.fSigma), mcc::DeviatoricStress(g.fSigma),
                        g.fPc, g.fV0, REAL(g.fType)});
    mcc::WriteCSV("embankment_gauss_" + tag + ".csv", {"x", "y", "z", "p_eff", "q", "pc", "v0", "type"}, rows);
}

inline EmbankmentConsolidation::TResult EmbankmentConsolidation::Run(EVariant variant, ETangentMode mode, bool output) {
    const auto clock0 = std::chrono::steady_clock::now();
    double outsec = 0.;
    auto timed = [&outsec](const std::function<void()> &f) { // time spent writing files
        const auto c0 = std::chrono::steady_clock::now();
        f();
        outsec += std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
    };
    const bool camclay = (variant == ECamClay);
    const std::string prefix = camclay ? "embankment" : "embankment_elastic";

    TPZGeoMesh *gmesh = CreateGeoMesh();
    if (output && camclay) timed([&] { mcc::WriteMeshCSV(gmesh, "embankment_mesh"); });
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, variant, mode, mat);

    // analysis: non-symmetric skyline matrix and LU decomposition (no renumbering, see TPZPoroElastoPlasticUPAnalysis)
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetPredictor(false);
    analysis.SetVerbose(fVerbose);

    TResult res;
    res.fMode = mode;
    res.fR0 = InitialResidual(analysis, mat, res.fBaseReaction0);
    analysis.ResetCounters();

    // VTK file series of every converged state (series time = index of the state, see Run)
    std::unique_ptr<mcc::TVTKSeries> vtk;
    if (output && fWriteVTK)
        timed([&] {
            vtk = std::make_unique<mcc::TVTKSeries>(mphys, mat, "vtk/" + prefix, prefix, "state",
                                                    std::vector<std::string>{"t", "lambda"});
        });
    const std::vector<bool> &ispressure = analysis.PressureEquations();
    res.fNEquations = mphys->NEquations();
    res.fNPressureNodes = std::count(ispressure.begin(), ispressure.end(), true);
    res.fNDisplacementNodes = (res.fNEquations - res.fNPressureNodes) / 3;

    // monitored points: top vertices of the face z = 0 at x = 0, 2, 4, 6 m and the elements of the zones pp1 and pp2
    std::array<int64_t, 4> top = {-1, -1, -1, -1};
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        if (std::fabs(x[1] - fHeight) > 1e-9 || std::fabs(x[2]) > 1e-9) continue;
        for (int k = 0; k < 4; ++k)
            if (std::fabs(x[0] - 2. * k) < 1e-9) top[k] = n;
    }
    for (int k = 0; k < 4; ++k)
        if (top[k] < 0) DebugStop();
    auto locate = [&](REAL x, REAL y) {
        TPZManVector<REAL, 3> xp = {x, y, 0.5 * fThickness}, qsi(3, 0.);
        TPZGeoEl *gel = mcc::LocatePoint(gmesh, EMatId, xp, qsi);
        if (!gel) DebugStop();
        return gel;
    };
    TPZGeoEl *zone1 = locate(0.5, 9.5), *zone2 = locate(1.5, 7.5);
    auto zonepressure = [&](TPZGeoEl *gel) { // mean of the pore pressure at the eight vertices
        REAL s = 0.;
        for (int i = 0; i < gel->NCornerNodes(); ++i) s += analysis.NodalValue(gel->NodeIndex(i), 1, 0);
        return s / gel->NCornerNodes();
    };
    if (output && camclay)
        timed([&] {
            std::ofstream out("embankment_monitor.csv");
            out << "name,x,y,z,xmin,xmax,ymin,ymax,zmin,zmax\n" << std::setprecision(12);
            for (int k = 0; k < 4; ++k) {
                TPZManVector<REAL, 3> x(3);
                gmesh->NodeVec()[top[k]].GetCoordinates(x);
                out << "s_x" << 2 * k << "," << x[0] << "," << x[1] << "," << x[2] << ",,,,,,\n";
            }
            for (auto zone : {std::make_pair("pp1", zone1), std::make_pair("pp2", zone2)}) {
                REAL lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
                for (int i = 0; i < zone.second->NCornerNodes(); ++i)
                    for (int c = 0; c < 3; ++c) {
                        lo[c] = std::min(lo[c], zone.second->NodePtr(i)->Coord(c));
                        hi[c] = std::max(hi[c], zone.second->NodePtr(i)->Coord(c));
                    }
                out << zone.first << "," << 0.5 * (lo[0] + hi[0]) << "," << 0.5 * (lo[1] + hi[1]) << ","
                    << 0.5 * (lo[2] + hi[2]) << "," << lo[0] << "," << hi[0] << "," << lo[1] << "," << hi[1] << ","
                    << lo[2] << "," << hi[2] << "\n";
            }
        });

    int stage = 0; // 0 undrained loading, 1 consolidation
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        if (stage == 1 && istep == 0) return; // the start of the consolidation is the end of the loading
        std::vector<REAL> row = {s.fTime, s.fLambda};
        // settlement, positive downwards (0 - u_y writes the initial value as 0 instead of -0)
        for (int k = 0; k < 4; ++k) row.push_back(0. - analysis.NodalValue(top[k], 0, 1));
        row.push_back(zonepressure(zone1));
        row.push_back(zonepressure(zone2));
        res.fHistory.push_back(row);
        if (vtk) timed([&] { vtk->Write(REAL(res.fHistory.size() - 1), {s.fTime, s.fLambda}); });
        if (output && camclay && stage == 1 && std::fabs(s.fTime - 1.e6) < 1e-6 * s.fTime)
            timed([&] { WriteState(analysis, mat, mphys, "t1e6", 1, res.fExcess1e6, res.fExcess1e6X); });
    };

    // undrained loading: increments of the load factor with Dt = 0
    std::vector<mcc::TAnalysis::TLoadState> stepsU, stepsC;
    for (int j = 1; j <= kNUndrained; ++j) stepsU.emplace_back(0., REAL(j) / kNUndrained, 0.);
    res.fConverged = analysis.Run(stepsU, monitor);
    res.fNUndrained = res.fHistory.size();
    CountTypesByLayer(mat, mphys, res.fLayerZ, res.fTypesUndrained, res.fTotalUndrained);
    res.fLogUndrained = analysis.StepLog();
    analysis.ClearStepLog();
    if (output && camclay && res.fConverged) {
        timed([&] { WriteState(analysis, mat, mphys, "undrained", 0, res.fExcessUndrained, res.fExcessUndrainedX); });
        res.fCvUndrained = ConsolidationCoefficient(mat, mphys, zone2);
    }

    // consolidation: t = 10^(2 + j/4), j = 0..24, from the end of the loading (t = 0, lambda = 1)
    if (res.fConverged) {
        stage = 1;
        for (int j = 0; j < kNConsolidation; ++j) stepsC.emplace_back(std::pow(10., 2. + j / 4.), 1., 0.);
        res.fConverged = analysis.Run(stepsC, monitor, mcc::TAnalysis::TLoadState(0., 1., 0.));
        CountTypesByLayer(mat, mphys, res.fLayerZ, res.fTypesFinal, res.fTotalFinal);
        res.fLogConsolidation = analysis.StepLog();
    }
    res.fNGlobalIterations = analysis.NGlobalIterations();
    res.fNBisections = analysis.NBisections();
    MatrixProfile(analysis, res, output && camclay);

    if (res.fConverged) {
        // reactions at the end of the consolidation (Lagrange nodes of the boundary faces), per metre of slab
        res.fReactionX0 = analysis.Reaction({ELeft}, 0, {EBottom});
        res.fReactionX20 = analysis.Reaction({ERight}, 0, {EBottom});
        res.fReactionBaseX = analysis.Reaction({EBottom}, 0);
        res.fReactionBaseY = analysis.Reaction({EBottom}, 1);
        // out-of-plane reactions of the whole faces z = 0 and z = 1 m (base line included): -/+ integral of sigma_zz
        res.fReactionZ0 = analysis.Reaction({EZ0}, 2);
        res.fReactionZ1 = analysis.Reaction({EZ1}, 2);
        mat->ForEachIntegrationPoint(mphys, [&](TPZCompEl *, int, const TPZVec<REAL> &, const TPZVec<REAL> &, REAL w,
                                                TPZElastoPlasticMem &mem) {
            res.fSigmaZZForce -= w * (mem.m_sigma.ZZ() - mem.m_elastoplastic_state.fpressure) / fThickness; // alpha_B = 1
        });
        PlaneStrainCheck(analysis, mphys, res);

        // heave of the top at the end (vertices of the face z = 0, Fig. 14d)
        std::vector<std::pair<REAL, REAL>> profile;
        for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
            TPZManVector<REAL, 3> x(3);
            gmesh->NodeVec()[n].GetCoordinates(x);
            if (std::fabs(x[1] - fHeight) < 1e-9 && std::fabs(x[2]) < 1e-9)
                profile.push_back({x[0], -analysis.NodalValue(n, 0, 1)});
        }
        std::sort(profile.begin(), profile.end());
        for (auto &pr : profile) {
            res.fHeave = std::min(res.fHeave, pr.second);
            if (pr.second < 0. && res.fHeaveStart < 0.) res.fHeaveStart = pr.first;
        }
        if (output && camclay) {
            REAL excess;
            std::array<REAL, 3> xexcess;
            timed([&] { WriteState(analysis, mat, mphys, "t1e8", 2, excess, xexcess); });
            res.fCvFinal = ConsolidationCoefficient(mat, mphys, zone2);
        }
    }

    if (output)
        timed([&] {
            // histories (Fig. 13) and convergence records (Sect. 6.7)
            mcc::WriteCSV(prefix + "_history.csv", {"t", "lambda", "s_x0", "s_x2", "s_x4", "s_x6", "pp1", "pp2"},
                          res.fHistory);
            std::vector<std::vector<REAL>> conv;
            int stg = 0;
            for (auto *log : {&res.fLogUndrained, &res.fLogConsolidation}) {
                for (size_t k = 0; k < log->size(); ++k) {
                    const auto &l = (*log)[k];
                    for (size_t i = 0; i < l.fResiduals.size(); ++i)
                        conv.push_back({REAL(stg), REAL(k + 1), l.fState.fTime, l.fState.fLambda, REAL(i + 1),
                                        l.fResiduals[i]});
                }
                stg++;
            }
            mcc::WriteCSV(prefix + "_convergence.csv", {"stage", "increment", "t", "lambda", "evaluation", "residual"},
                          conv);
            // plastic points per layer of integration points (Sect. 6.6)
            std::vector<std::vector<REAL>> lay;
            for (size_t l = 0; l < res.fLayerZ.size(); ++l) {
                std::vector<REAL> row = {REAL(l), res.fLayerZ[l]};
                for (int k = 0; k < 3; ++k) row.push_back(l < res.fTypesUndrained.size() ? res.fTypesUndrained[l][k] : 0.);
                for (int k = 0; k < 3; ++k) row.push_back(l < res.fTypesFinal.size() ? res.fTypesFinal[l][k] : 0.);
                lay.push_back(row);
            }
            mcc::WriteCSV(prefix + "_types.csv",
                          {"layer", "z", "elastic_undrained", "subcritical_undrained", "supercritical_undrained",
                           "elastic_t1e8", "subcritical_t1e8", "supercritical_t1e8"},
                          lay);
        });

    vtk.reset(); // the post-processing meshes refer to the meshes of the run
    mcc::DeleteMeshes(mphys);
    res.fSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    res.fOutputSeconds = outsec;
    return res;
}

inline void EmbankmentConsolidation::WriteSummary(const std::string &file) const {
    std::ofstream out(file);
    out << "quantity,value,v06_2d,unit\n" << std::setprecision(12);
    for (auto &s : fSummary) {
        out << s.fKey << "," << s.fValue << ",";
        if (std::isfinite(s.fRef2D)) out << s.fRef2D;
        out << "," << s.fUnit << "\n";
    }
}

inline void EmbankmentConsolidation::RunAll() {
    const auto clock0 = std::chrono::steady_clock::now();
    const REAL nan = std::nan("");
    const auto D = TPZPlasticStepModifiedCamClay::EConsistentTangent;
    fSummary.clear();
    auto types = [](const std::array<int, 3> &t) {
        return "[" + std::to_string(t[0]) + ", " + std::to_string(t[1]) + ", " + std::to_string(t[2]) + "]";
    };
    auto itslist = [](const std::vector<mcc::TAnalysis::TStepLog> &log) {
        std::string a = "[";
        for (size_t k = 0; k < log.size(); ++k) a += std::to_string(log[k].fResiduals.size()) + (k + 1 < log.size() ? "," : "]");
        return a;
    };
    auto reflist = [](const int *ref, int n) {
        std::string b = "[";
        for (int k = 0; k < n; ++k) b += std::to_string(ref[k]) + (k + 1 < n ? "," : "]");
        return b;
    };
    auto maxits = [](const std::vector<mcc::TAnalysis::TStepLog> &log) {
        size_t m = 0;
        for (auto &l : log) m = std::max(m, l.fResiduals.size());
        return int(m);
    };
    auto totalits = [](const std::vector<mcc::TAnalysis::TStepLog> &log) {
        size_t m = 0;
        for (auto &l : log) m += l.fResiduals.size();
        return int(m);
    };
    auto residuals = [](const mcc::TAnalysis::TStepLog &l) {
        std::ostringstream out;
        out << std::scientific << std::setprecision(1);
        for (REAL v : l.fResiduals) out << " " << v;
        return out.str();
    };
    const int iU = kNUndrained, iF = kNStates - 1; // rows of the end of the loading and of t = 1e8 s
    const auto &RU = kRefHistory[iU], &RF = kRefHistory[iF];

    TResult r;
    bool haveModel = false;
    if (fParts & EPartModel) {
        r = Run(ECamClay, D, true);
        const auto &h = r.fHistory;
        haveModel = r.fConverged && h.size() == size_t(kNStates);
        if (!haveModel) {
            std::cout << "Embankment: the analysis failed (" << h.size() << " converged states)" << std::endl;
        } else {
            int i6 = iU; // row of t = 1e6 s
            while (i6 < iF && std::fabs(h[i6][0] - 1.e6) > 1.) i6++;
            const auto &U = h[iU], &F = h[iF];
            std::cout << "Embankment on a Cam-Clay foundation (Sect. 6.6): 20 x 10 x 1 Hex20-Hex8 elements (slab of "
                      << fThickness << " m, u_z = 0 on z = 0 and z = " << fThickness << " m), " << r.fNDisplacementNodes
                      << " displacement nodes, " << r.fNPressureNodes << " pore pressure nodes, " << r.fNEquations
                      << " equations (" << r.fNFreeEquations << " after the elimination of the Dirichlet conditions)\n"
                      << "  2D model of v0.6: 20 x 10 Q8-Q4, 661 nodes, 231 pore pressure nodes, 1553 equations\n";
            std::cout << std::fixed << std::setprecision(1) << "run time " << r.fSeconds << " s (" << r.fOutputSeconds
                      << " s writing CSV and VTK files), " << r.fNGlobalIterations << " global iterations; skyline profile "
                      << r.fProfile << " entries, LU of " << std::scientific << std::setprecision(2) << r.fLUFlops
                      << " multiply-adds\n";
            std::cout << std::scientific << std::setprecision(1) << "initial state: largest residual at the free equations "
                      << r.fR0 << " kN (2D: 8.4e-13 kN/m, article v0.6 8e-13); " << std::fixed << std::setprecision(3)
                      << "vertical reaction of the base " << r.fBaseReaction0 << " kN (weight of the soil "
                      << fGammaSat * fWidth * fHeight * fThickness << ")\n";
            REAL sv, sh;
            GeostaticStress(0., sv, sh);
            TPZTensor<REAL> s0;
            s0.XX() = s0.ZZ() = sh;
            s0.YY() = sv;
            std::cout << std::setprecision(1) << "at the base: p'0 = " << mcc::MeanEffectiveStress(s0)
                      << " kPa, q0 = " << mcc::DeviatoricStress(s0) << " kPa (article 84, 69)\n";
            Summary("geostatic_p0_base", mcc::MeanEffectiveStress(s0), 84., "kPa");
            Summary("geostatic_q0_base", mcc::DeviatoricStress(s0), 69., "kPa");

            // Table 8
            const char *names[6] = {"settlement x = 0 (m)", "settlement x = 2 m (m)", "settlement x = 4 m (m)",
                                    "settlement x = 6 m (m)", "pore pressure pp1 (kPa)", "pore pressure pp2 (kPa)"};
            const char *keys[6] = {"s_x0", "s_x2", "s_x4", "s_x6", "pp1", "pp2"};
            const REAL flacU[6] = {0.140, 0.135, 0.055, -0.042, 18.1, 62.4}, flacF[6] = {0.193, 0.186, 0.104, 0.004, 5.1, 25.1};
            std::cout << "\nTable 8                 |       end of undrained loading      |              t = 1e8 s\n"
                      << "                        |    3D         2D v0.6       FLAC3D  |    3D         2D v0.6       FLAC3D\n";
            for (int i = 0; i < 6; ++i) {
                const int prec = (i < 4) ? 6 : 4;
                std::cout << std::left << std::setw(24) << names[i] << std::right << "|" << std::fixed
                          << std::setprecision(prec) << std::setw(11) << U[2 + i] << std::setw(12) << RU[1 + i]
                          << std::setprecision(3) << std::setw(10) << flacU[i] << "   |" << std::setprecision(prec)
                          << std::setw(11) << F[2 + i] << std::setw(12) << RF[1 + i] << std::setprecision(3)
                          << std::setw(10) << flacF[i] << "\n";
                Summary(std::string("table8_undrained_") + keys[i], U[2 + i], RU[1 + i], i < 4 ? "m" : "kPa");
                Summary(std::string("table8_t1e8_") + keys[i], F[2 + i], RF[1 + i], i < 4 ? "m" : "kPa");
            }
            // whole histories against the 2D model (Fig. 13)
            REAL ds = 0., dp = 0.;
            for (int k = 0; k < kNStates; ++k) {
                for (int c = 1; c <= 4; ++c) ds = std::max(ds, std::fabs(h[k][1 + c] - kRefHistory[k][c]));
                for (int c = 5; c <= 6; ++c) dp = std::max(dp, std::fabs(h[k][1 + c] - kRefHistory[k][c]));
            }
            std::cout << std::scientific << std::setprecision(2) << "largest difference from the 2D histories ("
                      << kNStates << " states): settlements " << ds << " m, pore pressures " << dp << " kPa\n";
            Summary("history_maxdiff_settlement", ds, 0., "m");
            Summary("history_maxdiff_pp", dp, 0., "kPa");
            // FLAC3D: undrained settlement at x = 0 (0.140 m) and heave at x = 6 m (-0.0416 m, from its history)
            const REAL dflac0 = 100. * (U[2] - 0.140) / 0.140, dflac6 = 100. * (U[5] + 0.0416) / 0.0416;
            std::cout << std::fixed << std::setprecision(1) << "difference from FLAC3D at the end of the loading: "
                      << "settlement x = 0 " << dflac0 << "% (v0.6: 9%), heave x = 6 m " << dflac6
                      << "% of 0.0416 m (v0.6: 7%)\n";
            Summary("flac3d_diff_settlement_x0_undrained", dflac0, 100. * (RU[1] - 0.140) / 0.140, "%");
            Summary("flac3d_diff_heave_x6_undrained", dflac6, 100. * (RU[4] + 0.0416) / 0.0416, "%");

            // reactions
            std::cout << std::fixed << std::setprecision(3) << "\nreactions at t = 1e8 s (kN, slab of 1 m): x = 0 "
                      << r.fReactionX0 << " (2D 898.088 kN/m), x = 20 m " << r.fReactionX20 << " (-795.210), base x "
                      << r.fReactionBaseX << " (-102.878), base y " << r.fReactionBaseY
                      << " (4800.000; 23 x 20 x 10 + 50 x 4 = 4800)\n"
                      << std::scientific << std::setprecision(1) << "sum of the horizontal reactions "
                      << r.fReactionX0 + r.fReactionX20 + r.fReactionBaseX << " kN (v0.6: zero to within 1e-5)\n"
                      << std::fixed << std::setprecision(3) << "out-of-plane reactions of the faces z = 0 and z = 1 m "
                      << "(all their nodes): " << r.fReactionZ0 << ", " << r.fReactionZ1
                      << " kN (-/+ the integral of the total sigma_zz over a face, the force that keeps the plane strain; "
                      << "integral over the integration points " << r.fSigmaZZForce
                      << " kN; 2D: 17165.275 kN/m from the stresses of data_aterro.pkl)\n";
            Summary("reaction_x0", r.fReactionX0, 898.088, "kN");
            Summary("reaction_x20", r.fReactionX20, -795.210, "kN");
            Summary("reaction_base_x", r.fReactionBaseX, -102.878, "kN");
            Summary("reaction_base_y", r.fReactionBaseY, 4800.000, "kN");
            Summary("reaction_sum_horizontal", r.fReactionX0 + r.fReactionX20 + r.fReactionBaseX, 9.4e-7, "kN");
            Summary("reaction_z0", r.fReactionZ0, 17165.275, "kN");
            Summary("reaction_z1", r.fReactionZ1, -17165.275, "kN");
            Summary("sigma_zz_force", r.fSigmaZZForce, 17165.275, "kN");
            Summary("initial_residual_max", r.fR0, 8.4e-13, "kN");
            Summary("initial_base_reaction", r.fBaseReaction0, 4600., "kN");
            std::cout << std::scientific << std::setprecision(1) << "plane strain at t = 1e8 s: largest |u(z=0) - u(z=1)| "
                      << r.fPlaneStrainU << " m, |p(z=0) - p(z=1)| " << r.fPlaneStrainP << " kPa, |u_z| at z = 0.5 m "
                      << r.fPlaneStrainUz << " m, |u(z=0.5) - u(z=0)| " << r.fPlaneStrainMid << " m\n";
            Summary("planestrain_du_faces", r.fPlaneStrainU, nan, "m");
            Summary("planestrain_dp_faces", r.fPlaneStrainP, nan, "kPa");
            Summary("planestrain_uz_mid", r.fPlaneStrainUz, nan, "m");
            Summary("planestrain_du_mid", r.fPlaneStrainMid, nan, "m");

            // plastic points
            std::cout << "\nintegration points [elastic, subcritical, supercritical]: end of loading "
                      << types(r.fTotalUndrained) << " of " << r.fTotalUndrained[0] + r.fTotalUndrained[1] + r.fTotalUndrained[2]
                      << " (2D [1800, 0, 0]), t = 1e8 s " << types(r.fTotalFinal) << " (2D [1628, 140, 32] of 1800)\n";
            for (size_t l = 0; l < r.fLayerZ.size(); ++l)
                std::cout << "  layer z = " << std::fixed << std::setprecision(4) << r.fLayerZ[l] << " m: end of loading "
                          << types(r.fTypesUndrained[l]) << ", t = 1e8 s " << types(r.fTypesFinal[l]) << "\n";
            Summary("points_total", r.fTotalFinal[0] + r.fTotalFinal[1] + r.fTotalFinal[2], 1800, "");
            Summary("points_elastic_t1e8", r.fTotalFinal[0], 1628, "");
            Summary("points_subcritical_t1e8", r.fTotalFinal[1], 140, "");
            Summary("points_supercritical_t1e8", r.fTotalFinal[2], 32, "");
            Summary("points_plastic_undrained", r.fTotalUndrained[1] + r.fTotalUndrained[2], 0, "");
            Summary("points_plastic_t1e8", r.fTotalFinal[1] + r.fTotalFinal[2], 172, "");
            for (size_t l = 0; l < r.fLayerZ.size(); ++l) {
                Summary("points_subcritical_t1e8_layer" + std::to_string(l), r.fTypesFinal[l][1], 140, "");
                Summary("points_supercritical_t1e8_layer" + std::to_string(l), r.fTypesFinal[l][2], 32, "");
            }

            // global iterations (Sect. 6.7)
            const REAL mU = mcc::MeanEvaluations(r.fLogUndrained), mC = mcc::MeanEvaluations(r.fLogConsolidation);
            std::cout << "evaluations of the residual per increment (v0.6: 4.3 undrained, 3.6 consolidation):\n"
                      << "  undrained     " << itslist(r.fLogUndrained) << " mean " << std::fixed << std::setprecision(2)
                      << mU << " (2D " << reflist(kRefItsUndrained, kNUndrained) << ")\n"
                      << "  consolidation " << itslist(r.fLogConsolidation) << " mean " << mC << " (2D "
                      << reflist(kRefItsConsolidation[0], kNConsolidation) << ")\n";
            int nbis = 0;
            for (auto *log : {&r.fLogUndrained, &r.fLogConsolidation})
                for (auto &l : *log) nbis += (l.fLevel > 0);
            std::cout << "  increments obtained by bisection: " << nbis << ", failed attempts " << r.fNBisections
                      << " (v0.6: none)\n";
            Summary("evaluations_undrained_mean", mU, 4.3, "");
            Summary("evaluations_consolidation_mean", mC, 3.56, "");
            Summary("evaluations_consolidation_max", maxits(r.fLogConsolidation), 5, "");
            Summary("global_iterations_total", r.fNGlobalIterations, 132, "");
            Summary("bisections", r.fNBisections, 0, "");
            std::cout << "  last undrained increment:" << residuals(r.fLogUndrained.back())
                      << "   (2D: 2.2e-02 6.6e-04 8.1e-07 2.0e-12)\n";
            for (size_t i = 0; i < r.fLogUndrained.back().fResiduals.size(); ++i)
                Summary("residual_last_undrained_" + std::to_string(i + 1), r.fLogUndrained.back().fResiduals[i],
                        i == 0 ? 2.2e-2 : i == 1 ? 6.6e-4 : i == 2 ? 8.1e-7 : i == 3 ? 2.0e-12 : nan, "");
            for (auto &l : r.fLogConsolidation)
                if (std::fabs(l.fState.fTime - 1.e6) < 1.) {
                    std::cout << "  step ending at t = 1e6 s:" << residuals(l) << "   (2D: 7.7e-05 5.6e-03 1.7e-05 2.0e-10)\n";
                    for (size_t i = 0; i < l.fResiduals.size(); ++i)
                        Summary("residual_step_t1e6_" + std::to_string(i + 1), l.fResiduals[i],
                                i == 0 ? 7.7e-5 : i == 1 ? 5.6e-3 : i == 2 ? 1.7e-5 : i == 3 ? 2.0e-10 : nan, "");
                }

            // Mandel-Cryer effect and consolidation rate near pp2
            int kmax = iU, kref = iU;
            for (int k = iU; k < kNStates; ++k) {
                if (h[k][7] > h[kmax][7]) kmax = k;
                if (kRefHistory[k][6] > kRefHistory[kref][6]) kref = k;
            }
            std::cout << std::fixed << std::setprecision(3) << "\nMandel-Cryer effect: pp2 rises from " << U[7] << " to "
                      << h[kmax][7] << " kPa at t = " << std::scientific << std::setprecision(2) << h[kmax][0]
                      << " s (2D " << std::fixed << std::setprecision(3) << RU[6] << " -> " << kRefHistory[kref][6]
                      << " kPa at " << std::scientific << std::setprecision(2) << kRefHistory[kref][0] << " s)\n";
            Summary("mandel_cryer_pp2_peak", h[kmax][7], kRefHistory[kref][6], "kPa");
            Summary("mandel_cryer_t_peak", h[kmax][0], kRefHistory[kref][0], "s");
            const REAL frac = (h[i6][2] - U[2]) / (F[2] - U[2]);
            const REAL fracref = (kRefHistory[i6][1] - RU[1]) / (RF[1] - RU[1]);
            std::cout << std::fixed << std::setprecision(3) << "t = 1e6 s: pp2 = " << h[i6][7] << " kPa (2D "
                      << kRefHistory[i6][6] << "), settlement x = 0 has developed " << std::setprecision(1) << 100. * frac
                      << "% of its consolidation part (2D " << 100. * fracref << "%)\n";
            Summary("pp2_t1e6", h[i6][7], kRefHistory[i6][6], "kPa");
            Summary("settlement_x0_share_t1e6", 100. * frac, 100. * fracref, "%");
            const REAL tcv = 2.5e5, dcv = 2.5; // time and depth of the estimate of the article
            std::cout << std::scientific << std::setprecision(3) << "consolidation coefficient of the skeleton in zone pp2, "
                      << "cv = k(K + 4G/3): " << r.fCvUndrained << " (end of loading) to " << r.fCvFinal
                      << " m2/s (t = 1e8 s) (2D 1.24e-06 to 2.58e-06); cv t/d2 at t = 2.5e5 s, d = 2.5 m: " << std::fixed
                      << std::setprecision(3) << r.fCvUndrained * tcv / (dcv * dcv) << " to "
                      << r.fCvFinal * tcv / (dcv * dcv) << " (v0.6 0.04-0.1)\n";
            Summary("cv_undrained", r.fCvUndrained, 1.24e-6, "m2/s");
            Summary("cv_t1e8", r.fCvFinal, 2.58e-6, "m2/s");
            Summary("cv_t_d2_undrained", r.fCvUndrained * tcv / (dcv * dcv), 0.050, "");
            Summary("cv_t_d2_t1e8", r.fCvFinal * tcv / (dcv * dcv), 0.103, "");

            // fields of Fig. 14
            std::cout << std::fixed << std::setprecision(2) << "largest excess pore pressure: end of loading "
                      << r.fExcessUndrained << " kPa at (" << r.fExcessUndrainedX[0] << ", " << r.fExcessUndrainedX[1]
                      << ", " << r.fExcessUndrainedX[2] << ") (2D 56.14 at (1, 9)); t = 1e6 s " << r.fExcess1e6
                      << " kPa at (" << r.fExcess1e6X[0] << ", " << r.fExcess1e6X[1] << ", " << r.fExcess1e6X[2]
                      << ") (2D 33.42 at (0, 8))\n";
            std::cout << std::setprecision(4) << "heave of the top at t = 1e8 s: up to " << -1e3 * r.fHeave
                      << " mm, from x = " << std::setprecision(1) << r.fHeaveStart << " m (2D 5.1132 mm from x = 9 m)\n";
            Summary("excess_max_undrained", r.fExcessUndrained, 56.14, "kPa");
            Summary("excess_max_undrained_x", r.fExcessUndrainedX[0], 1., "m");
            Summary("excess_max_undrained_y", r.fExcessUndrainedX[1], 9., "m");
            Summary("excess_max_t1e6", r.fExcess1e6, 33.42, "kPa");
            Summary("excess_max_t1e6_x", r.fExcess1e6X[0], 0., "m");
            Summary("excess_max_t1e6_y", r.fExcess1e6X[1], 8., "m");
            Summary("heave_max", -1e3 * r.fHeave, 5.1132, "mm");
            Summary("heave_start_x", r.fHeaveStart, 9., "m");
            Summary("nodes_displacement", r.fNDisplacementNodes, 661, "");
            Summary("nodes_pressure", r.fNPressureNodes, 231, "");
            Summary("equations", r.fNEquations, 1553, "");
            Summary("equations_free", r.fNFreeEquations, nan, "");
            Summary("skyline_profile", r.fProfile, nan, "");
            Summary("lu_multiply_adds", r.fLUFlops, nan, "");
            Summary("run_time_model", r.fSeconds, 15.7, "s");
            Summary("run_time_model_output", r.fOutputSeconds, nan, "s");
        }
        std::cout << std::flush;
    }

    // elastic variant (aterro_elastic.py)
    if (fParts & EPartElastic) {
        TResult e = Run(EElastic, D, true);
        if (!e.fConverged || e.fHistory.size() != size_t(kNStates)) {
            std::cout << "Embankment, elastic variant: the analysis failed" << std::endl;
        } else {
            REAL de = 0.;
            for (int k = 0; k < kNStates; ++k)
                for (int c = 1; c <= 4; ++c) de = std::max(de, std::fabs(e.fHistory[k][1 + c] - kRefElastic[k][c]));
            const auto &EF = e.fHistory.back();
            std::cout << std::fixed << std::setprecision(6) << "\nelastic foundation (pc = 1e7 kPa), t = 1e8 s: settlements "
                      << EF[2] << " " << EF[3] << " " << EF[4] << " " << EF[5] << " m (2D " << kRefElastic[iF][1] << " "
                      << kRefElastic[iF][2] << " " << kRefElastic[iF][3] << " " << kRefElastic[iF][4]
                      << "; v0.6 0.268 at x = 0); integration points " << types(e.fTotalFinal) << " (no yielding); "
                      << "largest difference from the 2D history " << std::scientific << std::setprecision(2) << de
                      << " m; " << std::fixed << std::setprecision(2) << "evaluations per step: undrained "
                      << mcc::MeanEvaluations(e.fLogUndrained) << ", consolidation " << mcc::MeanEvaluations(e.fLogConsolidation)
                      << "; run time " << std::setprecision(1) << e.fSeconds << " s\n";
            const char *keys[4] = {"s_x0", "s_x2", "s_x4", "s_x6"};
            for (int c = 0; c < 4; ++c) Summary(std::string("elastic_t1e8_") + keys[c], EF[2 + c], kRefElastic[iF][1 + c], "m");
            Summary("elastic_history_maxdiff", de, 0., "m");
            Summary("elastic_plastic_points", e.fTotalFinal[1] + e.fTotalFinal[2], 0, "");
            Summary("elastic_evaluations_undrained_mean", mcc::MeanEvaluations(e.fLogUndrained), nan, "");
            Summary("elastic_evaluations_consolidation_mean", mcc::MeanEvaluations(e.fLogConsolidation), nan, "");
            Summary("run_time_elastic", e.fSeconds, 14.5, "s");
            if (haveModel) {
                const REAL share = r.fHistory[iF][2] - EF[2];
                std::cout << std::setprecision(4) << "plastic share of the final settlement at x = 0: " << share
                          << " m (2D " << RF[1] - kRefElastic[iF][1] << ", v0.6 0.008)\n";
                Summary("plastic_share_x0", share, RF[1] - kRefElastic[iF][1], "m");
            }
        }
        std::cout << std::flush;
    }

    // comparison of the tangent operators (Table 10, embankment column; tangentes() of gen_data.py)
    if (fParts & EPartTangents) {
        std::vector<TResult> runs;
        for (auto mode : fTangentModes) {
            runs.push_back(Run(ECamClay, mode, false));
            const auto &t = runs.back();
            std::cout << "\ntangent " << TPZPlasticStepModifiedCamClay::TangentModeName(mode) << ": "
                      << (t.fConverged ? "converged" : "FAILED") << ", undrained " << itslist(t.fLogUndrained)
                      << ", consolidation " << itslist(t.fLogConsolidation) << ", " << t.fNGlobalIterations
                      << " global iterations, " << t.fNBisections << " failed attempts, " << std::fixed
                      << std::setprecision(1) << t.fSeconds << " s" << std::endl;
        }
        const TResult *ref = nullptr; // run with the consistent tangent
        for (auto &t : runs)
            if (t.fMode == D && t.fConverged) ref = &t;
        std::cout << "\nTable 10, embankment: evaluations of the residual per consolidation step (25 steps), "
                  << "mean (largest), and per undrained increment (10 increments)\n"
                  << "operator  undrained  consolidation    total  global its  failed  time (s) | 2D v0.6: consolidation  "
                     "total  time (s, Python)\n";
        std::ofstream tab("embankment_tangents.csv");
        tab << "mode,undrained_mean,undrained_max,consolidation_mean,consolidation_max,evaluations_total,"
               "global_iterations,failed_attempts,seconds,seconds_per_iteration,converged,final_s_x0,final_pp2,"
               "maxdiff_settlement_vs_D,maxdiff_pp_vs_D,maxreldiff_undrained_residuals_vs_D,"
               "v06_consolidation_mean,v06_consolidation_max,v06_global_iterations,v06_seconds_python\n"
            << std::setprecision(12);
        std::ofstream its("embankment_tangents_iterations.csv");
        its << "mode,stage,increment,t,lambda,evaluations,level\n" << std::setprecision(12);
        std::ofstream conv("embankment_tangents_convergence.csv");
        conv << "mode,stage,increment,t,lambda,evaluation,residual\n" << std::setprecision(12);
        for (auto &t : runs) {
            const std::string name = TPZPlasticStepModifiedCamClay::TangentModeName(t.fMode);
            const int row = RefRow(t.fMode);
            REAL dsx = nan, dpp = nan, dres = nan;
            if (ref && t.fConverged && t.fHistory.size() == ref->fHistory.size()) {
                dsx = dpp = 0.;
                for (size_t k = 0; k < t.fHistory.size(); ++k) {
                    for (int c = 2; c <= 5; ++c) dsx = std::max(dsx, std::fabs(t.fHistory[k][c] - ref->fHistory[k][c]));
                    for (int c = 6; c <= 7; ++c) dpp = std::max(dpp, std::fabs(t.fHistory[k][c] - ref->fHistory[k][c]));
                }
                dres = 0.;
                for (size_t k = 0; k < t.fLogUndrained.size() && k < ref->fLogUndrained.size(); ++k)
                    for (size_t i = 0; i < t.fLogUndrained[k].fResiduals.size() && i < ref->fLogUndrained[k].fResiduals.size(); ++i)
                        dres = std::max(dres, std::fabs(t.fLogUndrained[k].fResiduals[i] / ref->fLogUndrained[k].fResiduals[i] - 1.));
            }
            const REAL mU = mcc::MeanEvaluations(t.fLogUndrained), mC = mcc::MeanEvaluations(t.fLogConsolidation);
            const int tot = totalits(t.fLogUndrained) + totalits(t.fLogConsolidation);
            REAL refmean = nan;
            int refmax = -1;
            if (row >= 0) {
                refmean = 0.;
                for (int k = 0; k < kNConsolidation; ++k) {
                    refmean += kRefItsConsolidation[row][k];
                    refmax = std::max(refmax, kRefItsConsolidation[row][k]);
                }
                refmean /= kNConsolidation;
            }
            std::cout << std::left << std::setw(9) << name << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << mU << std::setw(10) << mC << " (" << maxits(t.fLogConsolidation) << ")"
                      << std::setw(9) << tot << std::setw(12) << t.fNGlobalIterations << std::setw(8) << t.fNBisections
                      << std::setw(10) << std::setprecision(1) << t.fSeconds << " | ";
            if (row >= 0)
                std::cout << std::setprecision(2) << std::setw(10) << refmean << " (" << refmax << ")" << std::setw(10)
                          << kRefItCount[row] << std::setw(10) << std::setprecision(1) << kRefTime[row];
            else
                std::cout << "      not run in v0.6";
            std::cout << "\n";
            auto num = [](REAL v) {
                std::ostringstream o;
                o << std::setprecision(12);
                if (std::isfinite(v)) o << v;
                return o.str();
            };
            const bool ok = t.fConverged && t.fHistory.size() == size_t(kNStates);
            tab << name << "," << mU << "," << maxits(t.fLogUndrained) << "," << mC << "," << maxits(t.fLogConsolidation)
                << "," << tot << "," << t.fNGlobalIterations << "," << t.fNBisections << "," << t.fSeconds << ","
                << t.fSeconds / std::max<int64_t>(1, t.fNGlobalIterations) << "," << (ok ? 1 : 0) << ","
                << num(ok ? t.fHistory.back()[2] : nan) << "," << num(ok ? t.fHistory.back()[7] : nan) << "," << num(dsx)
                << "," << num(dpp) << "," << num(dres) << "," << num(refmean) << "," << (row >= 0 ? std::to_string(refmax) : "")
                << "," << (row >= 0 ? std::to_string(kRefItCount[row]) : "") << ","
                << (row >= 0 ? num(kRefTime[row]) : "") << "\n";
            int stg = 0;
            for (auto *log : {&t.fLogUndrained, &t.fLogConsolidation}) {
                for (size_t k = 0; k < log->size(); ++k) {
                    const auto &l = (*log)[k];
                    its << name << "," << stg << "," << k + 1 << "," << l.fState.fTime << "," << l.fState.fLambda << ","
                        << l.fResiduals.size() << "," << l.fLevel << "\n";
                    for (size_t i = 0; i < l.fResiduals.size(); ++i)
                        conv << name << "," << stg << "," << k + 1 << "," << l.fState.fTime << "," << l.fState.fLambda
                             << "," << i + 1 << "," << l.fResiduals[i] << "\n";
                }
                stg++;
            }
            // consolidation steps whose number of evaluations differs from that with D, with both residual sequences
            if (ref && &t != ref && t.fConverged && t.fLogConsolidation.size() == ref->fLogConsolidation.size()) {
                for (size_t k = 0; k < t.fLogConsolidation.size(); ++k) {
                    const auto &a = t.fLogConsolidation[k], &b = ref->fLogConsolidation[k];
                    if (a.fResiduals.size() == b.fResiduals.size()) continue;
                    std::cout << "          consolidation step " << k + 1 << ": " << a.fResiduals.size() << " evaluations"
                              << residuals(a) << "; with D " << b.fResiduals.size() << residuals(b) << "\n";
                }
            }
            if (t.fMode == TPZPlasticStepModifiedCamClay::ETransposedTangent) {
                std::cout << "          undrained loading with D^T: " << itslist(t.fLogUndrained)
                          << ", largest relative change of the residuals from D " << std::scientific << std::setprecision(1)
                          << dres << " (2D: 0, all the points are elastic and the elastic tangent is symmetric)\n";
                Summary("transposed_undrained_mean", mU, 4.3, "");
                Summary("transposed_undrained_reldiff_residuals", dres, 0., "");
            }
        }
        std::cout << "(time: run time of the C++ code, one core, including the assembly of the reactions and the "
                     "post-processing of the history; 2D v0.6 times are those of the Python code)\n";
    }
    if (!fSummary.empty()) {
        Summary("run_time_total", std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count(), nan, "s");
        WriteSummary("embankment_summary.csv");
    }
}
