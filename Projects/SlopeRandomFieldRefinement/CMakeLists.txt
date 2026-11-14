# UnsatModel e UnsatAnalysis: classes genericas no material (u-p ou u-p_w-p_a, qualquer lei plastica)
add_pz_target(NAME SlopeRandomFieldRefinement
              SOURCES main.cpp UnsatModel.cpp UnsatModel.h UnsatAnalysis.cpp UnsatAnalysis.h HPMaterial.h
                      LayerGMesh.cpp LayerGMesh.h DrawdownTransient.cpp DrawdownTransient.h AnimationFiles.h PointLocator.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeRandomFieldRefinement PRIVATE slope_common) # Projects/Common: SlopeAnalysis
# pasta das malhas do Gmsh do talude do Ceron (as mesmas do relatorio), lida a partir de qualquer pasta de build
target_compile_definitions(SlopeRandomFieldRefinement PRIVATE
                           MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../LearningMohrCoulomb2/cluster/malhas"
                           PROJECT_DIR="${CMAKE_CURRENT_SOURCE_DIR}") # malhas/ deste projeto

# Chuva de Rahardjo et al. (2007), deterministico (rahardjo.cpp: copia do main.cpp com a TriGMesh): os tres solos da
# Tabela 1, transiente de 48 h com 24 h de chuva e o FS (SRM) de hora em hora com refinamento adaptativo
add_pz_target(NAME SlopeRahardjo
              SOURCES rahardjo.cpp UnsatModel.cpp UnsatModel.h UnsatAnalysis.cpp UnsatAnalysis.h HPMaterial.h
                      RainTransient.cpp RainTransient.h AnimationFiles.h PointLocator.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)

# Monte Carlo do rebaixamento rapido do Ceron com o transiente (Secoes 6.1 a 6.5 do artigo): o TransientDrawDown do
# main.cpp com campos aleatorios e Gamma (GIM) no fim do rebaixamento
add_pz_target(NAME SlopeDrawdownMonteCarlo
              SOURCES rebaixamento_mc.cpp UnsatModel.cpp UnsatModel.h UnsatAnalysis.cpp UnsatAnalysis.h HPMaterial.h
                      RandomField.cpp RandomField.h KLSeparavel.h DrawdownTransient.cpp DrawdownTransient.h
                      AnimationFiles.h PointLocator.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)

# Monte Carlo com chuva: campos aleatorios de c', phi' e k (KL analitico), transiente numa malha uniforme e o fator
# Gamma numa malha refinada em h, com as pressoes transferidas entre as duas (as malhas: TriGMesh ou em camadas,
# LayerGMesh). Nao usa o Common
add_pz_target(NAME SlopeRainMonteCarlo
              SOURCES montecarlo.cpp UnsatModel.cpp UnsatModel.h UnsatAnalysis.cpp UnsatAnalysis.h HPMaterial.h
                      RandomField.cpp RandomField.h KLSeparavel.h RainTransient.cpp RainTransient.h PointLocator.h
                      LayerGMesh.cpp LayerGMesh.h AnimationFiles.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
