# 3D slope factor of safety (tetrahedra or hexahedra) by gravity increase and strength reduction (RHW Mohr-Coulomb)
add_pz_target(NAME SlopeMohrCoulomb3D SOURCES main.cpp SlopeAnalysis3D.h REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
