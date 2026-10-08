# Slope factor of safety by gravity increase and strength reduction (RHW Mohr-Coulomb)
add_pz_target(NAME SlopeMohrCoulomb SOURCES main.cpp SlopeAnalysis.h SlopeModel.h REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
