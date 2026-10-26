# Unsaturated slope under rainfall (Rahardjo et al. 2007): coupled u - p_w - p_a Mohr-Coulomb analysis
# (TPZMatUnsaturatedUPwPa) and factor of safety (strength reduction and gravity increase) in time
add_pz_target(NAME SlopeRainfall SOURCES main.cpp UnsatSlope.h REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_compile_definitions(SlopeRainfall PRIVATE SLOPERAINFALL_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
