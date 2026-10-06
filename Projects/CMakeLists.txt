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

add_subdirectory(YieldSurfaceProjection)       # Figs. 1-2
add_subdirectory(TaylorTest)                   # Sect. 4.5, Fig. 3
add_subdirectory(RS2Triaxial)                  # Sect. 6.1, Fig. 4, Table 2
add_subdirectory(FrozenBulkModulus)            # Sect. 6.1, Table 3
add_subdirectory(FLAC3DTriaxial)               # Sect. 6.2, Fig. 5, Table 4
add_subdirectory(TerzaghiConsolidation)        # Sect. 6.3, Fig. 6, Table 5
add_subdirectory(AbaqusTriaxialConsolidation)  # Sects. 6.4 and 6.6, Figs. 7-9, Tables 6 and 8
add_subdirectory(EmbankmentConsolidation)      # Sect. 6.5, Figs. 10-12, Table 7
