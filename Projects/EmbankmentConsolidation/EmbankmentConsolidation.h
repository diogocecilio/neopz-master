/**
 * @file EmbankmentConsolidation.h
 * @brief Sect. 6.5 of the article: embankment loading on a Modified Cam-Clay foundation (FLAC3D example),
 * a plane strain u-p model with geostatic initial state, undrained loading and consolidation
 * (Figs. 10 to 12, Table 7 and the global iterations of Sect. 6.6).
 */
#pragma once

#include "MCCPaperTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

/**
 * @ingroup mccpaper
 * @brief Embankment on a Cam-Clay foundation (function aterro() of gen_data.py and aterro_elastic.py).
 *
 * Half of the problem, the domain \f$[0,20]\times[0,10]\f$ m, in plane strain with 20 x 10 Q8-Q4
 * elements (serendipity displacement, linear pore pressure) and 3 x 3 Gauss points (Fig. 10).
 *
 * Material (Table 1): \f$M=0.888\f$, \f$\lambda=0.161\f$, \f$\kappa=0.062\f$, \f$v_\lambda=2.858\f$,
 * porous elasticity with constant Poisson ratio \f$\nu=0.3\f$ and uniform preconsolidation pressure
 * \f$p'_{c0}=160\f$ kPa. Unit weights \f$\gamma_{sat}=23\f$ and \f$\gamma_w=10\f$ kN/m\f$^3\f$, pore fluid with
 * \f$K_f=2\times10^5\f$ kPa and porosity \f$n=0.3\f$ (\f$\alpha_B=1\f$, \f$1/M_B=n/K_f\f$), mobility
 * \f$k=10^{-9}\f$ m\f$^2\f$/(kPa s).
 *
 * Initial state at the depth \f$d=10-y\f$: hydrostatic pore pressure \f$p_w=\gamma_w d\f$, effective stresses
 * \f$\sigma'_v=-(\gamma_{sat}-\gamma_w)d\f$ and \f$\sigma'_h=0.7(-\gamma_{sat}d)+\gamma_w d\f$ (the total horizontal
 * stress is 0.7 times the total vertical stress), \f$p'_c=160\f$ kPa and the specific volume
 * \f$v_0=v_\lambda-\lambda\ln p'_{c0}+\kappa\ln(p'_{c0}/p'_0)\f$ at each integration point.
 *
 * Boundary conditions: \f$u_x=0\f$ at x = 0 (symmetry) and x = 20 m, fixed and impermeable base, drained top
 * (\f$p_w=0\f$) and the strip load q = 50 kPa on \f$0\le x\le 4\f$ m of the top, scaled by the load factor.
 *
 * Loading: ten undrained increments of the load factor (\f$\Delta t=0\f$), followed by 25 consolidation steps
 * with four steps per decade from \f$t=10^2\f$ to \f$10^8\f$ s; no predictor. Monitoring: settlements of the
 * top at x = 0, 2, 4 and 6 m and the pore pressures pp1 and pp2, averages of the four vertices of the
 * elements centred at (0.5, 9.5) and (1.5, 7.5) m.
 *
 * Variants (EVariant): the elastic model of aterro_elastic.py uses the same model with \f$p_c=10^7\f$ kPa
 * in the state (no yielding) and \f$v_0\f$ computed with \f$p'_{c0}=160\f$ kPa, which isolates the plastic share
 * of the settlement; as at the end of aterro(), the undrained loading is also repeated with the transpose
 * \f$D^T\f$ of the consistent tangent.
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
        EBottom = -1, ///< base: \f$u_x=u_y=0\f$, impermeable
        ERight = -2,  ///< x = 20 m: \f$u_x=0\f$
        ELeft = -4,   ///< x = 0 (symmetry): \f$u_x=0\f$
        ELoad = -5,   ///< loaded strip of the top, \f$0\le x\le 4\f$ m: traction (0, -q) times the load factor
        EPTop = -13   ///< whole top (coincident lines): drained, \f$p_w=0\f$
    };

    /** @brief Models solved by Run */
    enum EVariant {
        ECamClay,   ///< Modified Cam-Clay foundation, undrained loading and consolidation (aterro())
        EElastic,   ///< no yielding (\f$p_c=10^7\f$ kPa), undrained loading and consolidation (aterro_elastic.py)
        ETransposed ///< Cam-Clay, undrained loading only, with the transposed tangent \f$D^T\f$ (end of aterro())
    };

    /** @name Loading (Sect. 6.5) */
    /** @{ */
    static constexpr int kNUndrained = 10;    ///< undrained increments of the load factor
    static constexpr int kNConsolidation = 25; ///< consolidation steps, \f$t=10^{2+j/4}\f$ s, j = 0..24
    static constexpr int kNStates = kNUndrained + kNConsolidation + 1; ///< monitored states (initial state included)
    /** @} */

    /** @brief Results of a run */
    struct TResult {
        /** @brief Monitored states: rows (t, lambda, s(x=0), s(x=2), s(x=4), s(x=6), pp1, pp2), the initial state first */
        std::vector<std::vector<REAL>> fHistory;
        size_t fNUndrained = 0;                                  ///< rows of the undrained stage (initial state included)
        std::array<int, 3> fTypesUndrained = {0, 0, 0};         ///< elastic, subcritical, supercritical points at the end of the loading
        std::array<int, 3> fTypesFinal = {0, 0, 0};             ///< elastic, subcritical, supercritical points at t = 1e8 s
        std::vector<mcc::TAnalysis::TStepLog> fLogUndrained;     ///< convergence records of the undrained increments
        std::vector<mcc::TAnalysis::TStepLog> fLogConsolidation; ///< convergence records of the consolidation steps
        REAL fR0 = 0.;                ///< largest initial residual at the free displacement equations (nodal basis)
        REAL fBaseReaction0 = 0.;     ///< initial vertical reaction of the base (weight of the soil)
        REAL fReactionX0 = 0.;        ///< horizontal reaction at x = 0 (base nodes excluded), t = 1e8 s
        REAL fReactionX20 = 0.;       ///< horizontal reaction at x = 20 m (base nodes excluded), t = 1e8 s
        REAL fReactionBaseX = 0.;     ///< horizontal reaction of the base, t = 1e8 s
        REAL fReactionBaseY = 0.;     ///< vertical reaction of the base, t = 1e8 s
        REAL fExcessUndrained = 0.;   ///< largest excess pore pressure at the end of the loading (Fig. 12a)
        std::array<REAL, 2> fExcessUndrainedX = {0., 0.}; ///< node where it occurs
        REAL fExcess1e6 = 0.;         ///< largest excess pore pressure at t = 1e6 s (Fig. 12b)
        std::array<REAL, 2> fExcess1e6X = {0., 0.};       ///< node where it occurs
        REAL fHeave = 0.;             ///< largest heave of the top at t = 1e8 s (negative settlement, Fig. 12d)
        REAL fHeaveStart = -1.;       ///< first x of the top with heave at t = 1e8 s (-1 if none)
        REAL fCvUndrained = 0.;       ///< consolidation coefficient of the skeleton in zone pp2, end of the loading
        REAL fCvFinal = 0.;           ///< consolidation coefficient of the skeleton in zone pp2, t = 1e8 s
        int64_t fNEquations = 0;      ///< number of equations
        int64_t fNDisplacementNodes = 0; ///< number of displacement nodes (vertices and mid-edge nodes)
        int64_t fNPressureNodes = 0;  ///< number of pore pressure nodes
        double fSeconds = 0.;         ///< run time
        bool fConverged = true;       ///< false if an increment failed
    };

    /** @name Data of the example (Table 1 and Sect. 6.5) */
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
    REAL fLoadWidth = 4.;   ///< loaded strip \f$0\le x\le 4\f$ m
    REAL fRatioTotal = 0.7; ///< ratio of the total horizontal to the total vertical geostatic stress
    int fNx = 20;           ///< elements along x
    int fNy = 10;           ///< elements along y
    /** @} */

    /**
     * @brief Write the VTK file series of every converged state (mcc::TVTKSeries) of the Cam-Clay and elastic runs,
     * in the directories vtk/embankment and vtk/embankment_elastic (the command line argument "novtk" disables it)
     */
    bool fWriteVTK = true;

    /**
     * @brief Effective geostatic stresses (tension positive) at the height y
     * @param y height (m); the depth is \f$d=10-y\f$
     * @param[out] sv vertical effective stress \f$-(\gamma_{sat}-\gamma_w)d\f$
     * @param[out] sh horizontal effective stress \f$0.7(-\gamma_{sat}d)+\gamma_w d\f$
     */
    void GeostaticStress(REAL y, REAL &sv, REAL &sh) const {
        const REAL d = fHeight - y;
        sv = -(fGammaSat - fGammaW) * d;
        sh = fRatioTotal * (-fGammaSat * d) + fGammaW * d;
    }

    /** @brief Hydrostatic pore pressure \f$\gamma_w(10-y)\f$ (water table at the surface) at the height y */
    REAL Hydrostatic(REAL y) const { return fGammaW * (fHeight - y); }

    /**
     * @brief Geometric mesh: 20 x 10 quadrilaterals with the boundary lines (coincident lines EPTop and ELoad
     * on the loaded part of the top)
     * @return the geometric mesh
     */
    TPZGeoMesh *CreateGeoMesh();

    /**
     * @brief Displacement, pore pressure and multiphysics meshes, material, boundary conditions and initial
     * state (memory of the integration points and pore pressure field)
     * @param gmesh geometric mesh
     * @param variant model (elastic variant: pc = 1e7 kPa; transposed variant: tangent \f$D^T\f$)
     * @param[out] mat the u-p material of the clay layer
     * @return the multiphysics mesh
     */
    TPZMultiphysicsCompMesh *CreateCompMesh(TPZGeoMesh *gmesh, EVariant variant, mcc::TPoroMaterial *&mat);

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
     * \f$c_v=k(K+4G/3)\f$ with \f$K=v_0p'/\kappa\f$ and \f$G=3K(1-2\nu)/(2(1+\nu))\f$ (Sect. 6.5)
     * @param mat the u-p material
     * @param mphys multiphysics mesh
     * @param gel volume element
     * @return mean of \f$c_v\f$ over the integration points of gel (m2/s)
     */
    REAL ConsolidationCoefficient(mcc::TPoroMaterial *mat, TPZMultiphysicsCompMesh *mphys, TPZGeoEl *gel) const;

    /**
     * @brief Writes the fields of a state (Fig. 12): nodal VTK (displacement and pore pressure), VTK of the
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
                    int step, REAL &excess, std::array<REAL, 2> &xexcess);

    /**
     * @brief Solves a model (analysis, structural matrix, solver, increments) and post-processes it
     * @param variant ECamClay and EElastic: undrained loading and consolidation, with CSV files (and the
     * fields of Fig. 12 for ECamClay); ETransposed: undrained loading only, without output files
     * @return monitored history, convergence records and post-processed quantities of the run
     *
     * With fWriteVTK, the ECamClay and EElastic runs write the VTK file series of the 36 converged states
     * (mcc::TVTKSeries) in vtk/embankment and vtk/embankment_elastic: the series time is the index of the
     * state (0 initial state, 1 to 10 undrained increments, 11 to 35 consolidation steps) and the file
     * \<name\>_states.csv gives the time t and the load factor of each index.
     */
    TResult Run(EVariant variant);

    /** @brief Runs the three models and prints the comparison with the article and the Python code */
    void RunAll();

    /** @name Reference values (Python transcription: data_aterro.pkl and aterro_elastic.pkl) */
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
    /** @brief Evaluations of the residual per increment of the undrained loading (also with \f$D^T\f$) */
    static constexpr int kRefItsUndrained[kNUndrained] = {5, 5, 5, 4, 4, 4, 4, 4, 4, 4};
    /** @brief Evaluations of the residual per consolidation step */
    static constexpr int kRefItsConsolidation[kNConsolidation] = {2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
                                                                  4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 4};
    /** @} */
};

inline TPZGeoMesh *EmbankmentConsolidation::CreateGeoMesh() {
    const REAL xload = fLoadWidth;
    return mcc::CreateRectangleMesh(0., 0., fWidth, fHeight, fNx, fNy, EMatId,
                                    [xload](int side, const TPZVec<REAL> &xmid) {
                                        switch (side) {
                                        case 0:
                                            return std::vector<int>{EBottom};
                                        case 1:
                                            return std::vector<int>{ERight};
                                        case 2: // drained top; strip load on the faces with midpoint x <= 4 m
                                            if (xmid[0] <= xload + 1e-9) return std::vector<int>{EPTop, ELoad};
                                            return std::vector<int>{EPTop};
                                        default:
                                            return std::vector<int>{ELeft};
                                        }
                                    });
}

inline TPZMultiphysicsCompMesh *EmbankmentConsolidation::CreateCompMesh(TPZGeoMesh *gmesh, EVariant variant,
                                                                       mcc::TPoroMaterial *&mat) {
    const std::set<int> bcids = {EBottom, ERight, ELeft, ELoad, EPTop};
    TPZCompMesh *cmeshU = mcc::CreateDisplacementMesh(gmesh, 2, EMatId, bcids);
    TPZCompMesh *cmeshP = mcc::CreatePressureMesh(gmesh, 2, EMatId, bcids);

    TPZMultiphysicsCompMesh *mphys = new TPZMultiphysicsCompMesh(gmesh);
    mphys->SetDimModel(2);
    mat = new mcc::TPoroMaterial(EMatId, TPZMatPoroElastoPlasticUPBase::EPlaneStrain);
    mcc::TPlastic model;
    model.SetModifiedCamClay(fM, fLambda, fKappa);
    model.SetPorousElasticity();
    model.SetPoissonRatio(fNu);
    model.SetTransposedTangent(variant == ETransposed);
    mat->SetPlasticModel(model);
    mat->SetBiot(1., fPorosity / fKf);
    mat->SetPermeability(fMobility);
    TPZManVector<REAL, 3> b(2, 0.), rhowg(2, 0.);
    b[1] = -fGammaSat;
    rhowg[1] = -fGammaW;
    mat->SetBodyForce(b);
    mat->SetFluidWeight(rhowg);
    mat->SetIntegrationOrder(4); // 3 x 3 Gauss points
    mphys->InsertMaterialObject(mat);

    using B = TPZMatPoroElastoPlasticUPBase;
    TPZFNMatrix<4, STATE> val1(2, 2, 0.);
    TPZManVector<STATE, 3> val2(2, 0.), zero(1, 0.);
    // fixed base
    mphys->InsertMaterialObject(mat->CreateBC(mat, EBottom, B::EDirichletU, val1, val2));
    // u_x = 0 at x = 0 (symmetry) and x = 20 m
    val1(0, 0) = 1.;
    mphys->InsertMaterialObject(mat->CreateBC(mat, ELeft, B::EDirichletUDirectional, val1, val2));
    mphys->InsertMaterialObject(mat->CreateBC(mat, ERight, B::EDirichletUDirectional, val1, val2));
    // strip load, scaled by the load factor
    val1.Zero();
    TPZManVector<STATE, 3> load(2, 0.);
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

inline void EmbankmentConsolidation::WriteState(mcc::TAnalysis &an, mcc::TPoroMaterial *mat,
                                                TPZMultiphysicsCompMesh *mphys, const std::string &tag, int step,
                                                REAL &excess, std::array<REAL, 2> &xexcess) {
    // native NeoPZ post-processing of the nodal fields and the integration points
    mcc::WriteNodalVTK(an, 2, "embankment_nodal.vtk", step);
    mcc::WriteGaussPointsVTK(mat, mphys, "embankment_gauss_" + tag + ".vtk");
    // vertex values: displacement, pore pressure and excess pore pressure (Fig. 12a, b and d)
    TPZGeoMesh *gmesh = mphys->Reference();
    std::vector<std::vector<REAL>> rows;
    excess = -1e300;
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        const REAL p = an.NodalValue(n, 1, 0);
        const REAL pex = p - Hydrostatic(x[1]);
        rows.push_back({x[0], x[1], an.NodalValue(n, 0, 0), an.NodalValue(n, 0, 1), p, pex});
        if (pex > excess) {
            excess = pex;
            xexcess = {x[0], x[1]};
        }
    }
    mcc::WriteCSV("embankment_nodal_" + tag + ".csv", {"x", "y", "ux", "uy", "p", "p_excess"}, rows);
    // integration points (Fig. 12c)
    rows.clear();
    for (auto &g : mcc::GaussPoints(mat, mphys))
        rows.push_back({g.fX[0], g.fX[1], mcc::MeanEffectiveStress(g.fSigma), mcc::DeviatoricStress(g.fSigma), g.fPc,
                        g.fV0, REAL(g.fType)});
    mcc::WriteCSV("embankment_gauss_" + tag + ".csv", {"x", "y", "p_eff", "q", "pc", "v0", "type"}, rows);
}

inline EmbankmentConsolidation::TResult EmbankmentConsolidation::Run(EVariant variant) {
    const auto clock0 = std::chrono::steady_clock::now();
    TPZGeoMesh *gmesh = CreateGeoMesh();
    mcc::TPoroMaterial *mat = nullptr;
    TPZMultiphysicsCompMesh *mphys = CreateCompMesh(gmesh, variant, mat);

    // analysis: non-symmetric skyline matrix and LU decomposition (no renumbering, see TPZPoroElastoPlasticUPAnalysis)
    mcc::TAnalysis analysis(mphys, mat);
    TPZSkylineNSymStructMatrix<STATE> skyl(mphys);
    skyl.SetNumThreads(0);
    analysis.SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);
    analysis.SetPredictor(false);

    TResult res;
    res.fR0 = InitialResidual(analysis, mat, res.fBaseReaction0);

    // VTK file series of every converged state (series time = index of the state, see Run)
    std::unique_ptr<mcc::TVTKSeries> vtk;
    if (fWriteVTK && variant != ETransposed) {
        const std::string name = (variant == EElastic) ? "embankment_elastic" : "embankment";
        vtk = std::make_unique<mcc::TVTKSeries>(mphys, mat, "vtk/" + name, name, "state",
                                                std::vector<std::string>{"t", "lambda"});
    }
    const std::vector<bool> &ispressure = analysis.PressureEquations();
    res.fNEquations = mphys->NEquations();
    res.fNPressureNodes = std::count(ispressure.begin(), ispressure.end(), true);
    res.fNDisplacementNodes = (res.fNEquations - res.fNPressureNodes) / 2;

    // monitored points: top vertices at x = 0, 2, 4, 6 m and the elements of the zones pp1 and pp2
    std::array<int64_t, 4> top = {-1, -1, -1, -1};
    for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
        TPZManVector<REAL, 3> x(3);
        gmesh->NodeVec()[n].GetCoordinates(x);
        if (std::fabs(x[1] - fHeight) > 1e-9) continue;
        for (int k = 0; k < 4; ++k)
            if (std::fabs(x[0] - 2. * k) < 1e-9) top[k] = n;
    }
    auto locate = [&](REAL x, REAL y) {
        TPZManVector<REAL, 3> xp = {x, y, 0.}, qsi(2, 0.);
        TPZGeoEl *gel = mcc::LocatePoint(gmesh, EMatId, xp, qsi);
        if (!gel) DebugStop();
        return gel;
    };
    TPZGeoEl *zone1 = locate(0.5, 9.5), *zone2 = locate(1.5, 7.5);
    auto zonepressure = [&](TPZGeoEl *gel) { // mean of the pore pressure at the four vertices
        REAL s = 0.;
        for (int i = 0; i < gel->NCornerNodes(); ++i) s += analysis.NodalValue(gel->NodeIndex(i), 1, 0);
        return s / gel->NCornerNodes();
    };

    int stage = 0; // 0 undrained loading, 1 consolidation
    auto monitor = [&](int istep, const mcc::TAnalysis::TLoadState &s) {
        if (stage == 1 && istep == 0) return; // the start of the consolidation is the end of the loading
        std::vector<REAL> row = {s.fTime, s.fLambda};
        // settlement, positive downwards (0 - u_y writes the initial value as 0 instead of -0)
        for (int k = 0; k < 4; ++k) row.push_back(0. - analysis.NodalValue(top[k], 0, 1));
        row.push_back(zonepressure(zone1));
        row.push_back(zonepressure(zone2));
        res.fHistory.push_back(row);
        if (vtk) vtk->Write(REAL(res.fHistory.size() - 1), {s.fTime, s.fLambda});
        if (variant == ECamClay && stage == 1 && std::fabs(s.fTime - 1.e6) < 1e-6 * s.fTime)
            WriteState(analysis, mat, mphys, "t1e6", 1, res.fExcess1e6, res.fExcess1e6X);
    };

    // undrained loading: increments of the load factor with Dt = 0
    std::vector<mcc::TAnalysis::TLoadState> stepsU, stepsC;
    for (int j = 1; j <= kNUndrained; ++j) stepsU.emplace_back(0., REAL(j) / kNUndrained, 0.);
    res.fConverged = analysis.Run(stepsU, monitor);
    res.fNUndrained = res.fHistory.size();
    res.fTypesUndrained = mcc::CountTypes(mat, mphys);
    res.fLogUndrained = analysis.StepLog();
    analysis.ClearStepLog();
    if (variant == ECamClay && res.fConverged) {
        WriteState(analysis, mat, mphys, "undrained", 0, res.fExcessUndrained, res.fExcessUndrainedX);
        res.fCvUndrained = ConsolidationCoefficient(mat, mphys, zone2);
    }

    // consolidation: t = 10^(2 + j/4), j = 0..24, from the end of the loading (t = 0, lambda = 1)
    if (variant != ETransposed && res.fConverged) {
        stage = 1;
        for (int j = 0; j < kNConsolidation; ++j) stepsC.emplace_back(std::pow(10., 2. + j / 4.), 1., 0.);
        res.fConverged = analysis.Run(stepsC, monitor, mcc::TAnalysis::TLoadState(0., 1., 0.));
        res.fTypesFinal = mcc::CountTypes(mat, mphys);
        res.fLogConsolidation = analysis.StepLog();

        // reactions at the end of the consolidation (Lagrange nodes of the boundary lines)
        res.fReactionX0 = analysis.Reaction({ELeft}, 0, {EBottom});
        res.fReactionX20 = analysis.Reaction({ERight}, 0, {EBottom});
        res.fReactionBaseX = analysis.Reaction({EBottom}, 0);
        res.fReactionBaseY = analysis.Reaction({EBottom}, 1);

        // heave of the top at the end (Fig. 12d)
        std::vector<std::pair<REAL, REAL>> profile;
        for (int64_t n = 0; n < gmesh->NNodes(); ++n) {
            TPZManVector<REAL, 3> x(3);
            gmesh->NodeVec()[n].GetCoordinates(x);
            if (std::fabs(x[1] - fHeight) < 1e-9) profile.push_back({x[0], -analysis.NodalValue(n, 0, 1)});
        }
        std::sort(profile.begin(), profile.end());
        for (auto &pr : profile) {
            res.fHeave = std::min(res.fHeave, pr.second);
            if (pr.second < 0. && res.fHeaveStart < 0.) res.fHeaveStart = pr.first;
        }
        if (variant == ECamClay) {
            REAL excess;
            std::array<REAL, 2> xexcess;
            WriteState(analysis, mat, mphys, "t1e8", 2, excess, xexcess);
            res.fCvFinal = ConsolidationCoefficient(mat, mphys, zone2);
        }

        // histories (Fig. 11) and convergence records (Sect. 6.6)
        const std::string prefix = (variant == EElastic) ? "embankment_elastic" : "embankment";
        mcc::WriteCSV(prefix + "_history.csv", {"t", "lambda", "s_x0", "s_x2", "s_x4", "s_x6", "pp1", "pp2"},
                      res.fHistory);
        std::vector<std::vector<REAL>> conv;
        int stg = 0;
        for (auto *log : {&res.fLogUndrained, &res.fLogConsolidation}) {
            for (size_t k = 0; k < log->size(); ++k) {
                const auto &l = (*log)[k];
                for (size_t i = 0; i < l.fResiduals.size(); ++i)
                    conv.push_back({REAL(stg), REAL(k + 1), l.fState.fTime, l.fState.fLambda, REAL(i + 1), l.fResiduals[i]});
            }
            stg++;
        }
        mcc::WriteCSV(prefix + "_convergence.csv", {"stage", "increment", "t", "lambda", "evaluation", "residual"}, conv);
    }

    vtk.reset(); // the post-processing meshes refer to the meshes of the run
    mcc::DeleteMeshes(mphys);
    res.fSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    return res;
}

inline void EmbankmentConsolidation::RunAll() {
    TResult r = Run(ECamClay);
    const auto &h = r.fHistory;
    if (!r.fConverged || h.size() != size_t(kNStates)) {
        std::cout << "Embankment: the analysis failed (" << h.size() << " converged states)" << std::endl;
        return;
    }
    const int iU = kNUndrained, iF = kNStates - 1; // rows of the end of the loading and of t = 1e8 s
    int i6 = iU;                                    // row of t = 1e6 s
    while (i6 < iF && std::fabs(h[i6][0] - 1.e6) > 1.) i6++;
    const auto &U = h[iU], &F = h[iF];
    const auto &RU = kRefHistory[iU], &RF = kRefHistory[iF];

    std::cout << "Embankment on a Cam-Clay foundation (Sect. 6.5): 20 x 10 Q8-Q4 elements, " << r.fNDisplacementNodes
              << " nodes, " << r.fNPressureNodes << " pore pressure nodes (article 661, 231), " << r.fNEquations
              << " equations, run time " << std::fixed << std::setprecision(1) << r.fSeconds << " s\n";
    std::cout << std::scientific << std::setprecision(1) << "initial state: largest residual at the free equations "
              << r.fR0 << " kN/m (Python 8.4e-13, article 8e-13); " << std::fixed << std::setprecision(3)
              << "vertical reaction of the base " << r.fBaseReaction0 << " kN/m (weight of the soil "
              << fGammaSat * fWidth * fHeight << ")\n";
    {
        REAL sv, sh;
        GeostaticStress(0., sv, sh);
        TPZTensor<REAL> s0;
        s0.XX() = s0.ZZ() = sh;
        s0.YY() = sv;
        std::cout << std::setprecision(1) << "at the base: p'0 = " << mcc::MeanEffectiveStress(s0)
                  << " kPa, q0 = " << mcc::DeviatoricStress(s0) << " kPa (article 84, 69)\n";
    }

    // Table 7
    const char *names[6] = {"settlement x = 0 (m)", "settlement x = 2 m (m)", "settlement x = 4 m (m)",
                            "settlement x = 6 m (m)", "pore pressure pp1 (kPa)", "pore pressure pp2 (kPa)"};
    const REAL artU[6] = {0.153, 0.152, 0.067, -0.039, 33.1, 55.7}, artF[6] = {0.275, 0.269, 0.164, 0.023, 5.0, 25.1};
    const REAL flacU[6] = {0.140, 0.135, 0.055, -0.042, 18.1, 62.4}, flacF[6] = {0.193, 0.186, 0.104, 0.004, 5.1, 25.1};
    std::cout << "\nTable 7                 |        end of undrained loading         |                t = 1e8 s\n"
              << "                        |  this work      Python  article  FLAC3D |  this work      Python  article  FLAC3D\n";
    for (int i = 0; i < 6; ++i) {
        const int prec = (i < 4) ? 6 : 4, pref = (i < 4) ? 3 : 1;
        std::cout << std::left << std::setw(24) << names[i] << std::right << "|" << std::fixed << std::setprecision(prec)
                  << std::setw(11) << U[2 + i] << std::setw(12) << RU[1 + i] << std::setprecision(pref) << std::setw(9)
                  << artU[i] << std::setw(8) << flacU[i] << " |" << std::setprecision(prec) << std::setw(11) << F[2 + i]
                  << std::setw(12) << RF[1 + i] << std::setprecision(pref) << std::setw(9) << artF[i] << std::setw(8)
                  << flacF[i] << "\n";
    }
    // whole histories against Python (Fig. 11)
    REAL ds = 0., dp = 0.;
    for (int k = 0; k < kNStates; ++k) {
        for (int c = 1; c <= 4; ++c) ds = std::max(ds, std::fabs(h[k][1 + c] - kRefHistory[k][c]));
        for (int c = 5; c <= 6; ++c) dp = std::max(dp, std::fabs(h[k][1 + c] - kRefHistory[k][c]));
    }
    std::cout << std::scientific << std::setprecision(1) << "largest difference from the Python histories ("
              << kNStates << " states): settlements " << ds << " m, pore pressures " << dp << " kPa\n";
    // FLAC3D: undrained settlement at x = 0 (0.140 m) and heave at x = 6 m (-0.0416 m, from its history)
    std::cout << std::fixed << std::setprecision(1) << "difference from FLAC3D at the end of the loading: settlement x = 0 "
              << 100. * (U[2] - 0.140) / 0.140 << "% (article 9%), heave x = 6 m " << 100. * (U[5] + 0.0416) / 0.0416
              << "% of 0.0416 m (article 7%)\n";

    // reactions
    std::cout << std::fixed << std::setprecision(3) << "\nreactions at t = 1e8 s (kN/m): x = 0 " << r.fReactionX0
              << " (Python 898.088), x = 20 m " << r.fReactionX20 << " (-795.210), base x " << r.fReactionBaseX
              << " (-102.878), base y " << r.fReactionBaseY << " (4800.000; 23 x 20 x 10 + 50 x 4 = 4800)\n"
              << std::scientific << std::setprecision(1) << "sum of the horizontal reactions "
              << r.fReactionX0 + r.fReactionX20 + r.fReactionBaseX << " kN/m (article: zero to within 1e-5)\n";

    // plastic points
    auto types = [](const std::array<int, 3> &t) {
        return "[" + std::to_string(t[0]) + ", " + std::to_string(t[1]) + ", " + std::to_string(t[2]) + "]";
    };
    std::cout << "\nintegration points [elastic, subcritical, supercritical]: end of loading " << types(r.fTypesUndrained)
              << " (Python [1800, 0, 0]), t = 1e8 s " << types(r.fTypesFinal) << " (Python [1628, 140, 32])\n";

    // global iterations (Sect. 6.6)
    auto its = [](const std::vector<mcc::TAnalysis::TStepLog> &log, const int *ref, int nref) {
        std::string a = "[", b = "[";
        bool same = (int)log.size() == nref;
        for (size_t k = 0; k < log.size(); ++k) {
            a += std::to_string(log[k].fResiduals.size()) + (k + 1 < log.size() ? "," : "]");
            if (same && (int)log[k].fResiduals.size() != ref[k]) same = false;
        }
        for (int k = 0; k < nref; ++k) b += std::to_string(ref[k]) + (k + 1 < nref ? "," : "]");
        std::ostringstream out;
        out << a << " mean " << std::fixed << std::setprecision(2) << mcc::MeanEvaluations(log) << " (Python " << b
            << (same ? ", identical)" : ", DIFFERENT)");
        return out.str();
    };
    std::cout << "evaluations of the residual per increment (article: 4.3 undrained, 3.6 consolidation):\n  undrained     "
              << its(r.fLogUndrained, kRefItsUndrained, kNUndrained) << "\n  consolidation "
              << its(r.fLogConsolidation, kRefItsConsolidation, kNConsolidation) << "\n";
    int nbis = 0;
    for (auto *log : {&r.fLogUndrained, &r.fLogConsolidation})
        for (auto &l : *log) nbis += (l.fLevel > 0);
    std::cout << "  increments obtained by bisection: " << nbis << " (article: none)\n";
    auto residuals = [](const mcc::TAnalysis::TStepLog &l) {
        std::ostringstream out;
        out << std::scientific << std::setprecision(1);
        for (REAL v : l.fResiduals) out << " " << v;
        return out.str();
    };
    std::cout << "  last undrained increment:" << residuals(r.fLogUndrained.back())
              << "   (article 2.2e-02 6.6e-04 8.1e-07 2.0e-12)\n";
    for (auto &l : r.fLogConsolidation)
        if (std::fabs(l.fState.fTime - 1.e6) < 1.)
            std::cout << "  step ending at t = 1e6 s:" << residuals(l) << "   (article 7.7e-05 5.6e-03 1.7e-05 2.0e-10)\n";

    // Mandel-Cryer effect and consolidation rate near pp2
    int kmax = iU, kref = iU;
    for (int k = iU; k < kNStates; ++k) {
        if (h[k][7] > h[kmax][7]) kmax = k;
        if (kRefHistory[k][6] > kRefHistory[kref][6]) kref = k;
    }
    std::cout << std::fixed << std::setprecision(3) << "\nMandel-Cryer effect: pp2 rises from " << U[7] << " to "
              << h[kmax][7] << " kPa at t = " << std::scientific << std::setprecision(2) << h[kmax][0] << " s (Python "
              << std::fixed << std::setprecision(3) << RU[6] << " -> " << kRefHistory[kref][6] << " kPa at "
              << std::scientific << std::setprecision(2) << kRefHistory[kref][0]
              << " s; article 55.7 -> 58.8 kPa at 3.2e5 s)\n";
    const REAL frac = (h[i6][2] - U[2]) / (F[2] - U[2]);
    const REAL fracref = (kRefHistory[i6][1] - RU[1]) / (RF[1] - RU[1]);
    std::cout << std::fixed << std::setprecision(3) << "t = 1e6 s: pp2 = " << h[i6][7] << " kPa (Python "
              << kRefHistory[i6][6] << "; article: still at its undrained value), settlement x = 0 has developed "
              << std::setprecision(1) << 100. * frac << "% of its consolidation part (Python " << 100. * fracref
              << "%, article 31%)\n";
    const REAL tcv = 2.5e5, dcv = 2.5; // time and depth of the estimate of the article
    std::cout << std::scientific << std::setprecision(2) << "consolidation coefficient of the skeleton in zone pp2, "
              << "cv = k(K + 4G/3): " << r.fCvUndrained << " (end of loading) to " << r.fCvFinal
              << " m2/s (t = 1e8 s) (Python 1.24e-06 to 2.58e-06; article 1-2.5e-6); cv t/d2 at t = 2.5e5 s, d = 2.5 m: "
              << std::fixed << std::setprecision(3) << r.fCvUndrained * tcv / (dcv * dcv) << " to "
              << r.fCvFinal * tcv / (dcv * dcv) << " (article 0.04-0.1)\n";

    // fields of Fig. 12
    std::cout << std::fixed << std::setprecision(2) << "largest excess pore pressure: end of loading " << r.fExcessUndrained
              << " kPa at (" << r.fExcessUndrainedX[0] << ", " << r.fExcessUndrainedX[1]
              << ") (Python 56.14 at (1, 9), article 56.1); t = 1e6 s " << r.fExcess1e6 << " kPa at (" << r.fExcess1e6X[0]
              << ", " << r.fExcess1e6X[1] << ") (Python 33.42 at (0, 8), article 33.4)\n";
    std::cout << std::setprecision(4) << "heave of the top at t = 1e8 s: up to " << -1e3 * r.fHeave << " mm, from x = "
              << std::setprecision(1) << r.fHeaveStart << " m (Python 5.1132 mm from x = 9 m; article: up to 5 mm beyond x = 9 m)\n";

    // transposed tangent in the undrained loading (end of aterro())
    TResult t = Run(ETransposed);
    if (!t.fConverged) {
        std::cout << "Embankment, transposed tangent: the analysis failed" << std::endl;
        return;
    }
    REAL dres = 0.;
    for (size_t k = 0; k < t.fLogUndrained.size() && k < r.fLogUndrained.size(); ++k)
        for (size_t i = 0; i < t.fLogUndrained[k].fResiduals.size() && i < r.fLogUndrained[k].fResiduals.size(); ++i)
            dres = std::max(dres, std::fabs(t.fLogUndrained[k].fResiduals[i] / r.fLogUndrained[k].fResiduals[i] - 1.));
    std::cout << "\nundrained loading with the transposed tangent D^T: " << its(t.fLogUndrained, kRefItsUndrained, kNUndrained)
              << "; largest relative change of the residuals from D " << std::scientific << std::setprecision(1) << dres
              << " (Python 0: all the points are elastic and the elastic tangent is symmetric); run time "
              << std::fixed << t.fSeconds << " s\n";

    // elastic variant (aterro_elastic.py)
    TResult e = Run(EElastic);
    if (!e.fConverged || e.fHistory.size() != size_t(kNStates)) {
        std::cout << "Embankment, elastic variant: the analysis failed" << std::endl;
        return;
    }
    REAL de = 0.;
    for (int k = 0; k < kNStates; ++k)
        for (int c = 1; c <= 4; ++c) de = std::max(de, std::fabs(e.fHistory[k][1 + c] - kRefElastic[k][c]));
    const auto &EF = e.fHistory.back();
    std::cout << std::fixed << std::setprecision(6) << "\nelastic foundation (pc = 1e7 kPa), t = 1e8 s: settlements "
              << EF[2] << " " << EF[3] << " " << EF[4] << " " << EF[5] << " m (Python " << kRefElastic[iF][1] << " "
              << kRefElastic[iF][2] << " " << kRefElastic[iF][3] << " " << kRefElastic[iF][4] << "; article 0.268 at x = 0)\n"
              << std::setprecision(4) << "plastic share of the final settlement at x = 0: " << F[2] - EF[2]
              << " m (Python " << RF[1] - kRefElastic[iF][1] << ", article 0.008); integration points "
              << types(e.fTypesFinal) << " (no yielding); largest difference from the Python history " << std::scientific
              << std::setprecision(1) << de << " m; run time " << std::fixed << e.fSeconds << " s\n";
}
