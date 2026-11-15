# Cho (2010), Exemplo 2: talude c-phi seco com campos aleatorios correlacionados (c' e phi'), FS pelo SRM com refinamento
# adaptativo, Monte Carlo direto e amostragem por importancia. Usa as classes do SlopeRandomFieldRefinement (UnsatModel,
# UnsatAnalysis, RandomField, ImportanceSampling.h)
set(SRFR ${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement)
add_pz_target(NAME Cho2010
              SOURCES cho2010.cpp ${SRFR}/UnsatModel.cpp ${SRFR}/UnsatModel.h ${SRFR}/UnsatAnalysis.cpp
                      ${SRFR}/UnsatAnalysis.h ${SRFR}/HPMaterial.h ${SRFR}/RandomField.cpp ${SRFR}/RandomField.h
                      ${SRFR}/KLSeparavel.h ${SRFR}/ImportanceSampling.h
              REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
if(TARGET Cho2010)
    target_include_directories(Cho2010 PRIVATE ${SRFR})
endif()
