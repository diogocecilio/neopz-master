# Install script for directory: /home/runner/work/neopz-master/neopz-master/Pre

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Pre" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Pre/MMeshType.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZExtendGridDimension.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGenSpecialGrid.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZMHMeshControl.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZReadGIDGrid.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzhyperplane.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZAcademicGeoMesh.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGMSHReadMesh.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGeoMeshBuilder.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZMHMixedHybridMeshControl.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/doxpre.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzidentifyrefpattern.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzreadtetgen.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZAnalyticSolution.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGenGrid2D.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGmshReader.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZMHMixedMeshChannelControl.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzbuildmultiphysicsmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzpargrid.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/tpzhierarquicalgrid.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZGenGrid3D.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZHybridizeHDiv.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/TPZMHMixedMeshControl.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzdatafi.h"
    "/home/runner/work/neopz-master/neopz-master/Pre/pzreadmesh.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build_test/Pre/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
