# Install script for directory: /home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/External/FAD/TinyFad/Specializations" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadone.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadtwo.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadthree.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadfour.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadfive.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadsix.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadseven.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadeight.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadnine.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadten.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadeleven.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadtwelve.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadthirteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadfourteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadfifteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadsixteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadseventeen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadeighteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadnineteen.h"
    "/home/runner/work/neopz-master/neopz-master/Common/FAD/TinyFad/Specializations/tinyfadtwenty.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build_test/Common/FAD/TinyFad/Specializations/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
