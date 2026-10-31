# @file neopz/Projects/CMakeLists.txt
# Unsaturated soil (u - p_w - p_a, Mohr-Coulomb of Lira Cecilio, J. Eng. Math. 157:10, 2026): one test project per
# part of the verification plan of docs/biot_nao_saturado/formulacao/formulacao.tex. Each project is a Catch2
# executable registered in CTest (run "ctest" in the build folder of Projects; slow tests are tagged [.lento]).

cmake_minimum_required(VERSION 3.14)

project(UnsaturatedSoilTests LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_IS_PZ_BUILDTREE)
    find_package(NeoPZ REQUIRED)
endif()

# Catch2 (v2, single header)
if(NOT TARGET Catch2::Catch2)
    find_package(Catch2 2 QUIET)
endif()
if(NOT TARGET Catch2::Catch2)
    include(FetchContent)
    FetchContent_Declare(Catch2 GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v2.13.10)
    FetchContent_MakeAvailable(Catch2)
endif()

enable_testing()

# header-only code shared by the tests (model, meshes, analysis, factors of safety, reference solutions)
add_library(unsat_common INTERFACE)
target_include_directories(unsat_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/UnsatCommon)

# add_unsat_test(<name> <sources...>): Catch2 executable linked to pz, registered in CTest
function(add_unsat_test name)
    add_pz_target(NAME ${name} SOURCES ${ARGN} REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
    target_link_libraries(${name} PRIVATE unsat_common Catch2::Catch2)
    target_compile_definitions(${name} PRIVATE UNSAT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(${name} PROPERTIES FOLDER "UnsaturatedSoil")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_subdirectory(UnsatMaterialPoint) # ponto material: leis hidraulicas, Bishop, Mohr-Coulomb PV2 e tangente
add_subdirectory(UnsatWaterPhase)    # teste 2: agua com ar passivo e esqueleto rigido (Richards)
add_subdirectory(UnsatAirPhase)      # teste 3: agua e ar com esqueleto rigido
add_subdirectory(UnsatCoupled)       # testes 0 e 4: Jacobiano, Terzaghi, Mandel, Liakopoulos
add_subdirectory(UnsatSRM)           # teste 1: reducao de resistencia em solo seco
add_subdirectory(UnsatGIM)           # teste 1: aumento de gravidade em solo seco
add_subdirectory(LearningTherzagui)  # projeto de estudo do autor: coluna de Terzaghi
add_subdirectory(LearningMohrCoulomb) # projeto de estudo do autor: Mohr-Coulomb em sapata corrida e talude
add_subdirectory(LearningSlopeWater)  # projeto de estudo do autor: talude saturado, transiente e FS no tempo
