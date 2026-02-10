# Install script for directory: /home/runner/work/neopz-master/neopz-master/Save

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Save" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Save/TPZBFileStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZChunkTranslator.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZContBufferedStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZGeneralFStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZRestoredInstance.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/pzmd5stream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZChunkInTranslation.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZCircBufferedStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZFileStream.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZPersistenceManager.h"
    "/home/runner/work/neopz-master/neopz-master/Save/TPZSavable.h"
    "/home/runner/work/neopz-master/neopz-master/Save/doxsave.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build_test/Save/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
