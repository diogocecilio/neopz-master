# Install script for directory: /home/runner/work/neopz-master/neopz-master/Refine

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Refine" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Refine/TPZRefCube.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/TPZRefLinear.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/TPZRefPattern.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/TPZRefPatternDataBase.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/TPZRefPatternTools.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/doxrefine.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzrefpoint.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzrefprism.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzrefpyram.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzrefquad.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzreftetrahedra.h"
    "/home/runner/work/neopz-master/neopz-master/Refine/pzreftriangle.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/runner/work/neopz-master/neopz-master/build/Refine/RefPatterns/cmake_install.cmake")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build/Refine/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
