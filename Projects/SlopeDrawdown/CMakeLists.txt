# Rapid drawdown of the reservoir in front of the slope of SlopeMohrCoulomb: u-p poro-elastoplastic Mohr-Coulomb
# analysis and factor of safety (gravity increase and strength reduction) with the seepage forces
add_pz_target(NAME SlopeDrawdown SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
