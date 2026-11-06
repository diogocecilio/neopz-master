# @file neopz/Projects/CMakeLists.txt
# Verificacao do modelo nao saturado (u - p_w - p_a, Mohr-Coulomb associado, SRM, campos de Karhunen-Loeve) contra
# exemplos da literatura: um projeto por verificacao (V1 a V6, ver README.md). Cada projeto tem um executavel (o caso
# completo, para rodar no PC ou no cluster) e um executavel de testes Catch2 registrado no CTest (rode "ctest" na pasta
# de build de Projects; os testes demorados tem a marca [.lento] e rodam com "ctest -C lento" ou "<teste> [.lento]").

cmake_minimum_required(VERSION 3.14)

project(UnsaturatedSoilVerification LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_IS_PZ_BUILDTREE)
    find_package(NeoPZ REQUIRED)
endif()

# Catch2 (v2, um cabecalho so)
if(NOT TARGET Catch2::Catch2)
    find_package(Catch2 2 QUIET)
endif()
if(NOT TARGET Catch2::Catch2)
    include(FetchContent)
    FetchContent_Declare(Catch2 GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v2.13.10)
    FetchContent_MakeAvailable(Catch2)
endif()

enable_testing()

# codigo compartilhado (somente cabecalhos): modelo, malhas, analise, fatores de seguranca, campos de Karhunen-Loeve
add_library(unsat_common INTERFACE)
target_include_directories(unsat_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Comum)

# add_verif_test(<nome> <fontes...>): executavel Catch2 ligado a pz, registrado no CTest; VERIF_DIR = a pasta do projeto
function(add_verif_test name)
    add_pz_target(NAME ${name} SOURCES ${ARGN} REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
    target_link_libraries(${name} PRIVATE unsat_common Catch2::Catch2)
    target_compile_definitions(${name} PRIVATE VERIF_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
                               COMUM_DIR="${PROJECT_SOURCE_DIR}/Comum")
    set_target_properties(${name} PROPERTIES FOLDER "Verificacao")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# add_verif_exe(<nome> <fontes...>): o executavel do caso completo de um projeto
function(add_verif_exe name)
    add_pz_target(NAME ${name} SOURCES ${ARGN} REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
    target_link_libraries(${name} PRIVATE unsat_common)
    target_compile_definitions(${name} PRIVATE VERIF_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
                               COMUM_DIR="${PROJECT_SOURCE_DIR}/Comum")
    set_target_properties(${name} PROPERTIES FOLDER "Verificacao")
endfunction()

add_subdirectory(V1_LiakopoulosTrifasico) # Laloui et al. (2003): coluna de Liakopoulos com ar ativo (variantes i e ii)
add_subdirectory(V2_TaludesCho)           # Liao et al. (2025) / Cho (2010): SRM, campos KL e Monte Carlo em taludes secos
add_subdirectory(V3_SuperficieLivre)      # Sharma et al. (2021): percolacao com superficie livre e face de percolacao
add_subdirectory(V4_TaludeLiu)            # Liu et al. (2024): talude nao saturado sob chuva, campos de c', phi' e k_s
add_subdirectory(V5_EncostaDuasCamadas)   # Abbasov et al. (2024): encosta em duas camadas sob chuva, FS local
add_subdirectory(V6_MediaLocalKL)         # Yuan et al. (2021): reducao de variancia na media do campo KL num arco
