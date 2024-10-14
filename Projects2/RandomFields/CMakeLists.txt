file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main1_SRCS randomfield.cpp  ${sources} ${headers})

include_directories('src')
include_directories("/usr/include/eigen3")
add_executable(randomfield ${main1_SRCS})
target_link_libraries(randomfield pz)


