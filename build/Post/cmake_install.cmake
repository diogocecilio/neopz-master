# Install script for directory: /home/runner/work/neopz-master/neopz-master/Post

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
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Post" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Post/TPZDrawStyle.h"
    "/home/runner/work/neopz-master/neopz-master/Post/TPZVTKGeoMesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzdxmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphel1d.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphelq2dd.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphnode.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzv3dmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/tpzgraphelprismmapped.h"
    "/home/runner/work/neopz-master/neopz-master/Post/tpzgraphelt3d.h"
    "/home/runner/work/neopz-master/neopz-master/Post/TPZMeshSolution.h"
    "/home/runner/work/neopz-master/neopz-master/Post/doxpost.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphel1dd.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphelq3dd.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzmvmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pztrigraph.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzvisualmatrix.h"
    "/home/runner/work/neopz-master/neopz-master/Post/tpzgraphelpyramidmapped.h"
    "/home/runner/work/neopz-master/neopz-master/Post/TPZProjectEllipse.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphel.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphelq2d.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgraphmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pztrigraphd.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzvtkmesh.h"
    "/home/runner/work/neopz-master/neopz-master/Post/tpzgraphelt2dmapped.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzgradientreconstruction.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzpostprocanalysis.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzcompelpostproc.h"
    "/home/runner/work/neopz-master/neopz-master/Post/pzpostprocmat.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build/Post/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
