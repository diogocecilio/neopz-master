# Verificação do Cam-Clay modificado (TPZModifiedCamClay); requer BUILD_PLASTICITY_MATERIALS=ON.
# Depois de rodar o executável, compara_python.py compara os CSV gerados com o código Python de referência.
if(BUILD_PLASTICITY_MATERIALS)
    add_executable(PlasticityTestsCamClay main.cpp)
    target_link_libraries(PlasticityTestsCamClay pz)
endif()
