# Exemplo 1-1 da Secao 6.2 de Mika, Hughes, Schillinger, Wriggers e Hiemstra (2021, CMAME 379:113730): campo aleatorio
# com nucleo exponencial no meio cilindro vazado, por Karhunen-Loeve de Galerkin com TPZKarhunenLoeveMat e
# CalcStiffB/CalcStiffC do NeoPZ (KLCilindroMika, main.cpp). KLTaludeCeron (talude.cpp): a mesma KL de Galerkin em 2D com
# o nucleo exponencial separavel de Vargas Ceron et al. (2025), no retangulo e no talude, contra a forma fechada do
# retangulo de ../SlopeRandomField/KLSeparavel.h (ver README.md). Precisam do Eigen3 (problema generalizado de
# autovalores); sem ele, o projeto fica de fora e os outros configuram normalmente.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "KLCilindroMika: Eigen3 nao encontrado; o projeto fica de fora")
    return()
endif()
find_package(Threads REQUIRED)
add_pz_target(NAME KLCilindroMika SOURCES main.cpp)
target_link_libraries(KLCilindroMika PRIVATE Eigen3::Eigen Threads::Threads)
add_pz_target(NAME KLTaludeCeron SOURCES talude.cpp)
target_link_libraries(KLTaludeCeron PRIVATE Eigen3::Eigen Threads::Threads)
# malhas do talude (Gmsh 4.1) usadas por padrao
target_compile_definitions(KLTaludeCeron PRIVATE
    MALHAS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas")
