file(GLOB headers *.h)
file(GLOB sources *.cpp)

file (GLOB to_remove_from_sources main.cpp mainGUI.cpp poroelastoplastic.cpp mainPV.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})

set(main1_SRCS elasticity.cpp  ${sources} ${headers})

include_directories('src')
include_directories("/usr/include/eigen3")
add_executable(elasticity ${main1_SRCS})
target_link_libraries(elasticity pz)


