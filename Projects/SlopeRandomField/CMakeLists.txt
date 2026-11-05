# Talude de Rahardjo et al. (2007) sob chuva com campos aleatorios de Karhunen-Loeve de c', phi' e k (Vargas Ceron et al.
# 2025): Monte Carlo com fator de seguranca por reducao de resistencia. Uma realizacao por execucao (ver README.md).
add_pz_target(NAME SlopeRandomField SOURCES main.cpp KLSeparavel.h REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
# o main.cpp incluido (LearningMohrCoulomb2/cluster/fonte/main.cpp) usa MESH_DIR nas funcoes do projeto original
target_compile_definitions(SlopeRandomField PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/meshes")
