# Biblioteca comum dos projetos de talude: malhas e refinamento hp, materiais, caminho de equilibrio (GIM e SRM, por
# arc-length ou bisseccao), o transiente da chuva, os campos aleatorios (Karhunen-Loeve) e a classe SlopeAnalysis
add_library(slope_common STATIC
    SlopeMesh.cpp SlopeModel.cpp SlopeSolver.cpp SlopeRain.cpp SlopeFields.cpp SlopeAnalysis.cpp
    SlopeMesh.h SlopeModel.h SlopeSolver.h SlopeRain.h SlopeFields.h KLSeparavel.h SlopeAnalysis.h)
target_include_directories(slope_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(CMAKE_IS_PZ_BUILDTREE)
    target_link_libraries(slope_common PUBLIC pz)
else()
    target_link_libraries(slope_common PUBLIC NeoPZ::pz)
    target_include_directories(slope_common PUBLIC ${PZ_INCLUDE_DIRS})
endif()

# Teste dos quatro metodos do fator no talude seco de H = 10 m (Gamma = 1,777 e FS = 1,203 do artigo)
add_pz_target(NAME SlopeCommonMetodos SOURCES test_metodos.cpp REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_link_libraries(SlopeCommonMetodos PRIVATE slope_common)
target_compile_definitions(SlopeCommonMetodos PRIVATE
    MESH_TALUDE10="${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRandomFieldRefinement/malhas/talude_h10_grossa.msh")
