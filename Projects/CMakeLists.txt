# @file neopz/Projects/CMakeLists.txt
# Examples of the article "Return mapping for Modified Cam-Clay plasticity in rotated Haigh-Westergaard
# space with consistent tangent operator and coupled u-p consolidation" (one executable per example).

cmake_minimum_required(VERSION 3.14)

project(MCCPaperExamples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package if not in neopz buildtree (NeoPZConfig.cmake provides add_pz_target)
if(NOT CMAKE_IS_PZ_BUILDTREE)
    find_package(NeoPZ REQUIRED)
endif()

# header-only utilities shared by the examples
add_library(mcc_paper_common INTERFACE)
target_include_directories(mcc_paper_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Common)

# add_mcc_example(<name> <sources...>): one executable per example, linked to pz
function(add_mcc_example name)
    add_pz_target(NAME ${name} SOURCES ${ARGN} REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
    target_link_libraries(${name} PRIVATE mcc_paper_common)
    set_target_properties(${name} PROPERTIES FOLDER "MCCPaper")
endfunction()

# section, table and figure numbers of the article v0.6 (sections and tables are kept in v0.7; the figures after
# Fig. 3 may be renumbered when the finite element meshes are added)
add_subdirectory(YieldSurfaceProjection)       # Figs. 1-2
add_subdirectory(TaylorTest)                   # Sect. 4.5, Fig. 3; Table 10 (Taylor slopes of the operators)
add_subdirectory(RS2Triaxial)                  # Sect. 6.1, Fig. 4, Table 2
add_subdirectory(FrozenBulkModulus)            # Sect. 6.1, Table 3
add_subdirectory(IntegrationSchemes)           # Sect. 6.2, Fig. 5, Table 4
add_subdirectory(FLAC3DTriaxial)               # Sect. 6.3, Fig. 6, Table 5
add_subdirectory(TerzaghiConsolidation)        # Sect. 6.4, Fig. 7, Table 6
add_subdirectory(AbaqusTriaxialConsolidation)  # Sects. 6.5 and 6.7, Figs. 8-10, Tables 7, 9 and 10
add_subdirectory(EmbankmentConsolidation)      # Sects. 6.6 and 6.7, Figs. 11-13, Tables 8 and 10
