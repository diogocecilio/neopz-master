cmake_minimum_required(VERSION 3.8...3.15)

project (NeoPZExamples LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_subdirectory(Elasticity2D)

add_subdirectory(Plasticity2D)

add_subdirectory(geomecproblems)
