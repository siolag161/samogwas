cmake_minimum_required(VERSION 2.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -fPIC -O5 -fopenmp -mtune=corei7-avx") 

#set(PROJECT_DIR ..)
set(LIB_DIR "${PROJECT_SOURCE_DIR}/dependencies/" )
set(Boost_INCLUDE_DIR "${LIB_DIR}/boost_1_55_0/target/include/")
set(Boost_LIBRARY_DIR "${LIB_DIR}/boost_1_55_0/target/lib/")

set(ProBT_INCLUDE_DIR "${LIB_DIR}/probt_spl_2_3_99/include/")
set(ProBT_LIBRARY_DIR "${LIB_DIR}/probt_spl_2_3_99/lib/" )

add_subdirectory(tests)
add_subdirectory(projects)

