file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main1_SRCS main.cpp  ${sources} ${headers})

include_directories('src')
add_executable(PlasticityTestsMohrCoulomb ${main1_SRCS})
target_link_libraries(PlasticityTestsMohrCoulomb pz)

