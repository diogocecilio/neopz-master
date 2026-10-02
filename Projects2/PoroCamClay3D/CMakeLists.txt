# Exemplos u-p 3D com Cam-Clay modificado (TPZModifiedCamClay); requer BUILD_PLASTICITY_MATERIALS=ON.
#   AterroCamClay          aterro sobre fundação Cam-Clay (exemplo do FLAC3D)
#   TriaxialAbaqusCamClay  Abaqus Benchmarks 1.15.2, adensamento de um corpo de prova triaxial
# Os dois usam a estrutura nativa do NeoPZ (TPZMultiphysicsCompMesh, TPZMatPoroElastoPlastic3DMem,
# TPZLinearAnalysis, DefineGraphMesh/PostProcess, TPZPostProcAnalysis).
# Os executáveis *Serendipity (pasta serendipity/) são a referência com o elemento serendipity de 20 nós
# do Python/Abaqus, montado fora do TPZCompMesh (o H1 de ordem 2 do NeoPZ no hexaedro tem 27 funções).
if(BUILD_PLASTICITY_MATERIALS)
    set(nativo PoroCamClayNativo.h PoroCamClayNativo.cpp malhas.h malhas.cpp)
    add_executable(AterroCamClay aterro.cpp ${nativo})
    target_link_libraries(AterroCamClay pz)
    add_executable(TriaxialAbaqusCamClay triaxial_abaqus.cpp ${nativo})
    target_link_libraries(TriaxialAbaqusCamClay pz)

    set(serendipity serendipity/TPZPoroCamClayUP.h serendipity/TPZPoroCamClayUP.cpp malhas.h malhas.cpp)
    add_executable(AterroCamClaySerendipity serendipity/aterro_serendipity.cpp ${serendipity})
    target_include_directories(AterroCamClaySerendipity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(AterroCamClaySerendipity pz)
    add_executable(TriaxialAbaqusCamClaySerendipity serendipity/triaxial_serendipity.cpp ${serendipity})
    target_include_directories(TriaxialAbaqusCamClaySerendipity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(TriaxialAbaqusCamClaySerendipity pz)
endif()
