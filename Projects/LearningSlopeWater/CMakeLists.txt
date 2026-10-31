# Projeto de estudo do autor: talude saturado, transiente com elevacao do nivel d'agua e FS ao longo do tempo
add_pz_target(NAME LearningSlopeWater SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
# malhas do Gmsh do projeto LearningMohrCoulomb
target_compile_definitions(LearningSlopeWater PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb/meshes")
