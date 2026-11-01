# Projeto de estudo do autor: dessaturacao da coluna de areia de Liakopoulos (Secao 10 do resumo nao saturado),
# escrito passo a passo como o LearningMohrCoulomb
add_pz_target(NAME LearningLiakopoulos SOURCES main.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS
              # pos-processamento em Python, copiado para a pasta de build (rodar ao lado dos CSVs)
              FILES plot_liakopoulos.py coluna_liakopoulos_1d.py data/callari2009_fig3_fig4_digitalizado.csv)
# pasta da malha do Gmsh (.geo e .msh), lida a partir de qualquer pasta de build
target_compile_definitions(LearningLiakopoulos PRIVATE MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/meshes")
