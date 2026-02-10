# Install script for directory: /home/runner/work/neopz-master/neopz-master/Material

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/opt/neopz")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Material" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Material/TPZKarhunenLoeveMat.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMaterial.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMaterialT.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatBase.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZBndCond.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZBndCondT.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZBndCondBase.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMaterialData.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMaterialDataT.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatTypes.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatError.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatWithMem.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatLoadCases.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatGeneralisedEigenVal.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatSingleSpace.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatErrorSingleSpace.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatInterfaceSingleSpace.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatTransientSingleSpace.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatCombinedSpaces.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatErrorCombinedSpaces.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatInterfaceCombinedSpaces.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZNullMaterial.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZNullMaterialCS.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZLagrangeMultiplier.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZLagrangeMultiplierCS.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatTypes.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatKLCov2D.h"
    "/home/runner/work/neopz-master/neopz-master/Material/TPZMatKLKernel.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/Elasticity/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/ConsLaw/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/BlackOil/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/Projection/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/Poisson/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/Electromagnetics/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Material/DarcyFlow/cmake_install.cmake")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build/Material/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
