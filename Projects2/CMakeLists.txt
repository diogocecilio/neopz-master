cmake_minimum_required(VERSION 3.8...3.15)

project (NeoPZExamples LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


add_subdirectory(Elasticity2D)

add_subdirectory(geomecprobdeter)

add_subdirectory(eigensolver)


add_subdirectory(eigensolverkl)

add_subdirectory(randomfem)

add_subdirectory(ElasticityRandomField)

add_subdirectory(GeoMecRandFields)

add_subdirectory(GeoMecDeterm)

add_subdirectory(GeoMecRandFieldsMonteCarlo)

add_subdirectory(PoroElastic)

add_subdirectory(PoroElastoPlastic)

add_subdirectory(PlasticityTests)
