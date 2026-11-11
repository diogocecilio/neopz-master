add_pz_target(NAME SlopeRandomFieldRefinement SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeRandomFieldRefinement PRIVATE slope_common) # Projects/Common: SlopeAnalysis
# pasta das malhas do Gmsh do talude do Ceron (as mesmas do relatorio), lida a partir de qualquer pasta de build
target_compile_definitions(SlopeRandomFieldRefinement PRIVATE
                           MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas"
                           PROJECT_DIR="${CMAKE_CURRENT_SOURCE_DIR}") # malhas/ deste projeto
