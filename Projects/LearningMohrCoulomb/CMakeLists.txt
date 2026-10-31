# Projeto de estudo do autor: Mohr-Coulomb em sapata corrida e talude, escrito passo a passo
add_pz_target(NAME LearningMohrCoulomb SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
# pasta das malhas do Gmsh (.geo e .msh), lida a partir de qualquer pasta de build
target_compile_definitions(LearningMohrCoulomb PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/meshes")
