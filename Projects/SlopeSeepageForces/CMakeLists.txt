# Slope stability under the seepage forces of a rapid drawdown (Ceron et al., IJNAMG 2025, Section 4.3): parametric
# slope meshes, steady anisotropic seepage, seepage force fields (FE and analytical), kinematic limit analysis and FEM
# gravity-increase stability (SlopeAnalysis.h); commands documented in main.cpp
add_pz_target(NAME SlopeSeepageForces
              SOURCES main.cpp SeepageCommands.cpp LimitAnalysisCommands.cpp FEMCommands.cpp SelfTests.cpp
                      Commands.h FigureCases.h ResumableCSV.h
                      SlopeGeometry.h DelaunayMesher.h AnisotropicDarcy.h SeepageFE.h SeepageForceField.h
                      AnalyticalSeepage.h LimitAnalysis.h FEMStability.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
