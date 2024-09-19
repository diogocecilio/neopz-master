file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main1_SRCS main.cpp  ${sources} ${headers})

include_directories('src')
include_directories("/usr/include/eigen3")
add_executable(geomecproblems ${main1_SRCS})
target_link_libraries(geomecproblems pz)
