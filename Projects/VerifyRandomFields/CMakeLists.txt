# Verificacao do campo aleatorio de coesao do talude de Vargas Ceron et al. (2025): KL de Galerkin do NeoPZ numa malha
# uniforme fina, forma fechada do retangulo (../SlopeRandomField/KLSeparavel.h, a do SlopeRandomField) e Nystrom, com os
# mesmos numeros aleatorios, comparados em todos os pontos de integracao da malha de Galerkin (ver README.md). Precisa do
# Eigen3 (autopares densos); sem ele, o projeto fica de fora e os outros configuram normalmente.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "VerifyRandomFields: Eigen3 nao encontrado; o projeto fica de fora")
    return()
endif()
find_package(Threads REQUIRED)
add_pz_target(NAME VerifyRandomFields SOURCES main.cpp ../SlopeRandomField/KLSeparavel.h)
target_link_libraries(VerifyRandomFields PRIVATE Eigen3::Eigen Threads::Threads)
# diretorio do projeto: as malhas uniformes do talude (malhas/) usadas por padrao
target_compile_definitions(VerifyRandomFields PRIVATE VRF_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
