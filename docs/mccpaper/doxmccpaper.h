/**
 * \file
 * @brief Documentation page and group of the implementation of the Modified Cam-Clay u-p article.
 */

/**
 * @defgroup mccpaper Examples of the Modified Cam-Clay u-p article
 * @brief Projects that reproduce the numerical examples of the article "Return mapping for Modified
 * Cam-Clay plasticity in rotated Haigh-Westergaard space with consistent tangent operator and coupled
 * u-p consolidation" (D. Lira Cecilio), one executable per example (folder Projects).
 */

/**
 * @page mccpaperpage Modified Cam-Clay in rotated Haigh-Westergaard space and coupled u-p consolidation
 *
 * @section mccpaper_overview Overview
 *
 * The NeoPZ implementation of the article follows the structure of the Wolfram Language packages
 * camclay-perf.m (constitutive model) and poro-camclay-fem.m (u-p finite elements) and of their Python
 * transcription, and reuses the native classes of the library:
 *
 * | Article / WL routine | NeoPZ |
 * |---|---|
 * | HardeningCC, PhiCC, ResCC, JacCC, dResdTrialCC, ProjectHWCC, GradCC (Sect. 4, eqs. (11)-(22), (A.1)) | TPZYCModifiedCamClayRHW |
 * | TrialStressCC, ProjectStressCC, ComputedDep (Algorithm 1, eqs. (8)-(10), (16)) | TPZPlasticStepModifiedCamClay (spectral decomposition with TPZTensor::EigenSystem) |
 * | ComputeBN, ContributePorous, ContributePlasticity (eqs. (25)-(28)) | TPZMatPoroElastoPlasticUP (multiphysics material with memory, TPZMatWithMem) |
 * | SolveStepUP, AdvanceUP, IterativeProcessUP, ReactionByMarker (Sect. 5.4) | TPZPoroElastoPlasticUPAnalysis (derived from TPZLinearAnalysis) |
 * | SubdivideQuadMesh, BoxMesh3D, QuarterCylinderMesh, LocatePoint, GaussPointWeights, TriaxialPointCC, closed forms (Appendix B) | Projects/Common/MCCPaperTools.h (namespace mcc) |
 *
 * Main conventions: tension positive; Voigt order of TPZTensor (XX, XY, XZ, YY, YZ, ZZ), the same of the
 * article; engineering shear strains; the consistent tangent is assembled column by column (the
 * transposition of Sect. 3 of the article is built in, TPZPlasticStepModifiedCamClay::SetTransposedTangent
 * reproduces the wrong operator for the comparisons of Sects. 4.5 and 6.6).
 *
 * The finite element spaces are native NeoPZ H1 spaces built with TPZMultiphysicsCompMesh: quadratic
 * displacement and linear pore pressure. With the face and volume connects of quadrilaterals and hexahedra
 * reduced to order 1, the hierarchical space of order 2 is exactly the serendipity Q8/Hex20 space of the
 * article (Taylor-Hood pairs Q8-Q4 and Hex20-Hex8). The Gauss rules are chosen by the material
 * (TPZMatPoroElastoPlasticUP::SetIntegrationOrder: 2 x 2 reduced or 3 x 3 full). The non-symmetric
 * monolithic systems are solved with TPZSkylineNSymStructMatrix and the LU decomposition, and the equations
 * of the Dirichlet conditions are eliminated with TPZEquationFilter.
 *
 * @section mccpaper_examples Examples (folder Projects)
 *
 * | Project | Article |
 * |---|---|
 * | YieldSurfaceProjection | Figs. 1 and 2: yield surface and closest-point projection in the meridian plane |
 * | TaylorTest | Sect. 4.5, Fig. 3: Taylor test of the consistent tangent |
 * | RS2Triaxial | Sect. 6.1, Fig. 4, Table 2: drained triaxial tests of the RS2 manual |
 * | FrozenBulkModulus | Sect. 6.1, Table 3: exact integration of the porous law versus frozen bulk modulus |
 * | FLAC3DTriaxial | Sect. 6.2, Fig. 5, Table 4: triaxial tests of FLAC3D with one u-p element |
 * | TerzaghiConsolidation | Sect. 6.3, Fig. 6, Table 5: Terzaghi consolidation (Q8-Q4 and Hex20-Hex8) |
 * | AbaqusTriaxialConsolidation | Sects. 6.4 and 6.6, Figs. 7-9, Tables 6 and 8: Abaqus benchmark 1.15.2 |
 * | EmbankmentConsolidation | Sect. 6.5, Figs. 10-12, Table 7: embankment on a Cam-Clay foundation |
 *
 * Each project is a class written in the style of the NeoPZ examples (geometric mesh, computational
 * meshes, analysis with structural matrix and solver, incremental solution and post-processing). The
 * projects are built with -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON.
 *
 * @section mccpaper_paraview Viewing the solution in ParaView
 *
 * EmbankmentConsolidation and AbaqusTriaxialConsolidation write every converged state (every increment
 * or time step) as VTK file series with mcc::TVTKSeries, one directory per run under vtk/ (the command
 * line argument novtk skips them):
 *  - \<prefix\>_nodal.vtk.series: displacement and pore pressure at the nodes, written by the native graph
 *    mesh of the multiphysics mesh (TPZAnalysis::DefineGraphMesh once, SetStep and PostProcess per state);
 *  - \<prefix\>_intpoints.vtk.series: the variables of the integration points of TPZMatPoroElastoPlasticUP
 *    (TPZMatPoroElastoPlasticUP::ESolutionVar: p', q, p_c, type of response, volumetric strain, specific
 *    volume, effective and total stresses, principal effective stresses), projected element by element on a
 *    discontinuous mesh by TPZPostProcAnalysis, as in the footing example of NeoPZ (SetPostProcessVariables
 *    once, TransferSolution and PostProcess per state);
 *  - \<prefix\>_gausspoints.vtk.series: the integration points as a point cloud (mcc::WriteGaussPointsVTK);
 *  - \<prefix\>_states.csv: the time, load factor or displacement of each state.
 *
 * The .series files (JSON, "file-series-version" 1.0) give ParaView the time of each file. TPZPostProcAnalysis
 * supports a TPZMultiphysicsCompMesh whose material is a combined-space material with memory: each
 * post-processing element (TPZCompElPostProc) refers to the multiphysics element, builds the vector of
 * material data of the atomic spaces as TPZMultiphysicsCompEl::CalcStiff does, sets the memory index of each
 * integration point and calls TPZMatCombinedSpacesT::Solution. Its order is n-1 for n x n (x n) Gauss points,
 * so the projection is the Lagrange extrapolation of the integration point values (the interpolation of
 * mcc::StressAtPoint). The READMEs of the two projects describe the files and the steps in ParaView.
 */
