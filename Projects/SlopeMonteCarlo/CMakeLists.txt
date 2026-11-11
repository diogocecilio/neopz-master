add_pz_target(NAME SlopeMonteCarlo SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeMonteCarlo PRIVATE slope_common) # Projects/Common: SlopeAnalysis e RandomFields
# malhas do Gmsh do talude do Ceron (as mesmas do relatorio), lidas a partir de qualquer pasta de build
target_compile_definitions(SlopeMonteCarlo PRIVATE
                           MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas")

# As Secoes 6.1 a 6.5 do Ceron no cluster (fermi): o mesmo modelo, com as pastas e os arquivos da campanha ceron_sat
add_pz_target(NAME SlopeMonteCarloCluster SOURCES cluster.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeMonteCarloCluster PRIVATE slope_common)
target_compile_definitions(SlopeMonteCarloCluster PRIVATE
                           MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas")
