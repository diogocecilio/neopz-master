# Exemplo 1-1 da Secao 6.2 de Mika, Hughes, Schillinger, Wriggers e Hiemstra (2021, CMAME 379:113730): campo aleatorio
# com nucleo exponencial no meio cilindro vazado, por Karhunen-Loeve de Galerkin com TPZKarhunenLoeveMat e
# CalcStiffB/CalcStiffC do NeoPZ (ver README.md). Precisa do Eigen3 (problema generalizado de autovalores); sem ele, o
# projeto fica de fora e os outros configuram normalmente.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "KLCilindroMika: Eigen3 nao encontrado; o projeto fica de fora")
    return()
endif()
find_package(Threads REQUIRED)
add_pz_target(NAME KLCilindroMika SOURCES main.cpp)
target_link_libraries(KLCilindroMika PRIVATE Eigen3::Eigen Threads::Threads)
