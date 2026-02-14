# - Config file for the PZ package
# It defines the following variables
#  PZ_INCLUDE_DIRS - include directories for using PZ
#  PZ_LIBRARIES    - PZ library to link against



include(/opt/neopz/lib/cmake/neopz/NeoPZTargets.cmake)
include(/opt/neopz/lib/cmake/neopz/add_pz_target.cmake)
include(/opt/neopz/lib/cmake/neopz/check_pz_opt.cmake)
set(PZ_INCLUDE_DIRS "/opt/neopz/pz/include")

## Compute paths
set(PZ_BRANCH "copilot/verify-transfer-solution-function")
set(PZ_REVISION "a0ef19e")
set(PZ_REVISION_DATE "Sat Feb 14 14:06:32 2026")

# These are IMPORTED targets created by PZTargets.cmake
if(NOT TARGET NeoPZ::pz)
  message(FATAL_ERROR "Could not find PZ libs!")
endif()
set(PZ_LIBRARIES NeoPZ::pz)
message(STATUS "PZ_INCLUDE_DIRS: ${PZ_INCLUDE_DIRS}")
message(STATUS "Link to: ${PZ_LIBRARIES}")

set_property(GLOBAL PROPERTY PZ_REAL_TYPE double)
set_property(GLOBAL PROPERTY PZ_STATE_TYPE double)
set_property(GLOBAL PROPERTY PZ_BUILD_PLASTICITY_MATERIALS OFF)
set_property(GLOBAL PROPERTY PZ_USING_BOOST OFF)
set_property(GLOBAL PROPERTY PZ_USING_TBB OFF)
set_property(GLOBAL PROPERTY PZ_USING_LAPACK OFF)
set_property(GLOBAL PROPERTY PZ_USING_PAPI OFF)
set_property(GLOBAL PROPERTY PZ_USING_MKL OFF)
set_property(GLOBAL PROPERTY PZ_LOG OFF)
