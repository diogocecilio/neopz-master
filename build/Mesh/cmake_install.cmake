# Install script for directory: /home/runner/work/neopz-master/neopz-master/Mesh

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Mesh" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Mesh/tpzcompmeshreferred.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/doxmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzconnect.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgeoelbc.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompMeshTools.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZGeoMeshTools.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcheckgeom.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcreateapproxspace.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgeoel.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzintel.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZGeoElement.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZInterfaceEl.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcheckmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzelchdivbound2.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgeoelrefless.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzinterpolationspace.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElDisc.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZGeoElement.h.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZMultiphysicsCompMesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcheckrestraint.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzelchdiv.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgeoelrefless.h.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzmultiphysicscompel.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/tpzgeoelmapped.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZMultiphysicsInterfaceEl.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzelctemp.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgeoelside.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzmultiphysicselement.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElHCurl.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/tpzgeoelrefpattern.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZOneShapeRestraint.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcompel.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzelementgroup.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzreducedspace.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElHDivCollapsed.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/tpzgeoelrefpattern.h.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcompelwithmem.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzelmat.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZElementMatrixT.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzgnode.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElLagrange.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZGeoElSideAncestors.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzcondensedcompel.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzsubcmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZGeoElSidePartition.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/pzflowcmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/tpzagglomeratemesh.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZAgglomerateEl.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElH1.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElKernelHDiv.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZCompElKernelHDiv3D.h"
    "/home/runner/work/neopz-master/neopz-master/Mesh/TPZHCurlEquationFilter.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build/Mesh/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
