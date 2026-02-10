# Install script for directory: /home/runner/work/neopz-master/neopz-master/Refine/RefPatterns

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/pz/include/Refine/RefPatterns" TYPE FILE FILES
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Face_Side_8.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Face_Side_8_5quad.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4_5_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4_5_5_6_6_7_7.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_4_6_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_5_5_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_5_5_7.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_5_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_5_6_7.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Quad_Rib_Side_4_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Face_Side_6.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_3.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_3_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_3_4_permuted.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_3_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_3_5_permuted.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_4_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_OnlyTriang_Side_4_5_permuted.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/2D_Triang_Rib_Side_3_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa2Tetrahedras.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Face_20.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Face_20_22_25.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_08.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_10.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_12.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_14.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_16.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_16_17_19.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_16_18_18.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_17.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_17_17_19_26.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_17_18.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_16_18.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_8_10.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_Rib_Side_8_9_11_16_17_19.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_directional_1face.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_directional_2faces.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_directional_3faces.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Hexa_directional_3faces_cross.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_1RaisedEdge_Directional.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_2RaisedEdges.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_2lados_lateral.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_5_7.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_BaseEdge_Directional.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_BaseFace.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_Rib_Side_10.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Piram_Rib_Side_5.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Prism_Rib_Side_10.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Prism_Rib_Side_10_11_12_14.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Prism_Rib_Side_12.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Prism_Rib_Side_15.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Prism_Rib_Side_6_7_12_13.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_1lado.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_2ladosopostos.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_2nodes_000.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_2nodes_003.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_5ribs.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_Rib_Side_4.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_UniformTetra.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/3D_Tetra_hexas.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/PrismNode.rpt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/CMakeLists.txt"
    "/home/runner/work/neopz-master/neopz-master/Refine/RefPatterns/VERY_IMPORTANT!!!!!.txt"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/runner/work/neopz-master/neopz-master/build/Refine/RefPatterns/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
