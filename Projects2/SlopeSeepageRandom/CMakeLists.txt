# Talude com forças de percolação e variabilidade espacial (Vargas Ceron et al., IJNAMG 2025): campos aleatórios
# por Karhunen-Loève, Darcy desacoplado, elastoplasticidade (Mohr-Coulomb e Cam-Clay modificado) e Monte Carlo.
# Requer BUILD_PLASTICITY_MATERIALS=ON e USING_LAPACK=ON (autoproblema KL).
if(BUILD_PLASTICITY_MATERIALS)
    add_executable(SlopeSeepageRandom main.cpp
        SlopeGeometry.h SlopeGeometry.cpp
        KLRandomField.h KLRandomField.cpp
        SeepageProblem.h SeepageProblem.cpp
        SlopeStability.h SlopeStability.cpp)
    target_link_libraries(SlopeSeepageRandom pz)
endif()
