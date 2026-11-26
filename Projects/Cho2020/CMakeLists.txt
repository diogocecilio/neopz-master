# Reproducao de Cho (2020), Landslides 17:79-91: talude infinito com rocha rasa sob chuva constante, na coluna normal
# ao talude do Kang2020. Usa as classes do SlopeRandomFieldRefinement e o talude_infinito.h do Kang2020
set(SRF ${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement)
add_pz_target(NAME Cho2020
              SOURCES main.cpp ${SRF}/UnsatModel.cpp ${SRF}/UnsatModel.h ${SRF}/HPMaterial.h
                      ${SRF}/RainTransient.cpp ${SRF}/RainTransient.h ${SRF}/UnsatAnalysis.cpp ${SRF}/UnsatAnalysis.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_include_directories(Cho2020 PRIVATE ${SRF} ${CMAKE_CURRENT_SOURCE_DIR}/../Kang2020)
