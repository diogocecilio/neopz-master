cmake_minimum_required(VERSION 3.8...3.15)


file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main1_SRCS main.cpp  ${sources} ${headers})

include_directories('src')
include_directories("/usr/include/eigen3")
add_executable(geomecprobdeter ${main1_SRCS})
target_link_libraries(geomecprobdeter pz)


