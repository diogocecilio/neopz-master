# Reproducao de Cho (2016), Eng. Geol. 211:184-197: talude sob chuva de 3 k_s, fluxo so da agua contra agua-ar, FS pelo
# SRM ao longo do tempo. Usa as classes do SlopeRandomFieldRefinement (UnsatModel, UnsatAnalysis, HPMaterial,
# LayerGMesh, PointLocator, RainTransient), compiladas aqui de la
set(SRF ${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement)
add_pz_target(NAME Cho2016
              SOURCES main.cpp ${SRF}/UnsatModel.cpp ${SRF}/UnsatModel.h ${SRF}/UnsatAnalysis.cpp ${SRF}/UnsatAnalysis.h
                      ${SRF}/HPMaterial.h ${SRF}/LayerGMesh.cpp ${SRF}/LayerGMesh.h ${SRF}/PointLocator.h
                      ${SRF}/RainTransient.cpp ${SRF}/RainTransient.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_include_directories(Cho2016 PRIVATE ${SRF})

# a coluna vertical 1D (coluna.cpp): convergencia de malha da frente sp e tp, separando a malha da fisica
add_pz_target(NAME Cho2016Coluna
              SOURCES coluna.cpp ${SRF}/UnsatModel.cpp ${SRF}/UnsatModel.h ${SRF}/UnsatAnalysis.cpp ${SRF}/UnsatAnalysis.h
                      ${SRF}/HPMaterial.h ${SRF}/RainTransient.cpp ${SRF}/RainTransient.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_include_directories(Cho2016Coluna PRIVATE ${SRF})
