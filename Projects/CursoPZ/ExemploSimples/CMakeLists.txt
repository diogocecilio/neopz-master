file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(ExemploSimples ${main_SRCS})

target_link_libraries(ExemploSimples pz)