# Slope stability under the seepage forces of a rapid drawdown (Ceron et al., IJNAMG 2025, Section 4.3): parametric
# slope meshes, steady anisotropic seepage, seepage force field and FEM gravity-increase stability (SlopeAnalysis.h)
add_pz_target(NAME SlopeSeepageForces
              SOURCES main.cpp SlopeGeometry.h DelaunayMesher.h AnisotropicDarcy.h SeepageFE.h SeepageForceField.h
                      FEMStability.h LimitAnalysis.h AnalyticalSeepage.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
