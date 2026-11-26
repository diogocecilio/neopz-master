# Reproducao de Kang, Cho, Kim e Go (2020), Water 12:812: talude infinito raso sob a chuva de 26-27/07/2011 (Umyeon),
# so a agua (ar passivo) contra o bifasico agua-ar (ar ativo), escrita passo a passo. O transiente usa as classes ja
# validadas do SlopeRandomFieldRefinement (UnsatModel, HPMaterial, RainTransient), compiladas aqui de la
set(SRF ${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement)
add_pz_target(NAME Kang2020
              SOURCES main.cpp kang_dados.h talude_infinito.h
                      ${SRF}/UnsatModel.cpp ${SRF}/UnsatModel.h ${SRF}/HPMaterial.h ${SRF}/RainTransient.cpp
                      ${SRF}/RainTransient.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_include_directories(Kang2020 PRIVATE ${SRF})
