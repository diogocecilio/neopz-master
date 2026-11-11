# Campo de percolacao semi-analitico de Vargas Ceron et al. (2025) no GIM do rebaixamento rapido, comparado com o MEF
add_pz_target(NAME SlopeSemiAnalytical SOURCES main.cpp SemiAnalyticalSeepage.cpp SemiAnalyticalSeepage.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeSemiAnalytical PRIVATE slope_common)
target_compile_definitions(SlopeSemiAnalytical PRIVATE
    MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas")
