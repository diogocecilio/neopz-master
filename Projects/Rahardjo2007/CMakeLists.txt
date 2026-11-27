# Reproducao da Serie B de Rahardjo et al. (2007), J. Geotech. Geoenviron. Eng. 133(12):1532-1543: talude homogeneo a
# 45 graus sob chuva de 24 h, tres solos e sete intensidades por solo, modelo (iii) (esqueleto rigido, ar passivo), FS
# pelo SRM ao longo do tempo. Usa as classes do SlopeRandomFieldRefinement (UnsatModel, UnsatAnalysis, HPMaterial,
# LayerGMesh, PointLocator, RainTransient), compiladas aqui de la
set(SRF ${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement)
add_pz_target(NAME Rahardjo2007
              SOURCES main.cpp ${SRF}/UnsatModel.cpp ${SRF}/UnsatModel.h ${SRF}/UnsatAnalysis.cpp ${SRF}/UnsatAnalysis.h
                      ${SRF}/HPMaterial.h ${SRF}/LayerGMesh.cpp ${SRF}/LayerGMesh.h ${SRF}/PointLocator.h
                      ${SRF}/RainTransient.cpp ${SRF}/RainTransient.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_include_directories(Rahardjo2007 PRIVATE ${SRF})
