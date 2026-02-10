cmake_minimum_required(VERSION 3.8...3.15)


file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main1_SRCS main.cpp  ${sources} ${headers})

include_directories('src')
add_executable(hprefine ${main1_SRCS})
target_link_libraries(hprefine pz)



