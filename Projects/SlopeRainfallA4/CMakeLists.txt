# Caso A4 (contraexemplos a Rahardjo et al. 2007) com a sucção inicial sem limite: programa autocontido, com malha
# geométrica, malhas computacionais, material, solver, análise, chuva, FS e saída VTK para o ParaView explícitos.
# Usa UnsatSlope.h e as malhas do projeto SlopeRainfall.
add_pz_target(NAME SlopeRainfallA4 SOURCES main.cpp ../SlopeRainfall/UnsatSlope.h REQUIRED PZ_BUILD_PLASTICITY_MATERIALS)
target_compile_definitions(SlopeRainfallA4 PRIVATE
                           SLOPERAINFALL_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../SlopeRainfall/meshes")
