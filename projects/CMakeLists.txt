cmake_minimum_required(VERSION 2.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O5 -fPIC -fopenmp -mtune=corei7-avx") 

set(LIB_DIR "${PROJECT_SOURCE_DIR}/dependencies/" )
set(Boost_INCLUDE_DIR "${LIB_DIR}/boost_1_55_0/target/include/")
set(Boost_LIBRARY_DIR "${LIB_DIR}/boost_1_55_0/target/lib/")

set(ProBT_INCLUDE_DIR "${LIB_DIR}/probt_spl_2_3_0/include/")
set(ProBT_LIBRARY_DIR "${LIB_DIR}/probt_spl_2_3_0/lib/" )

#add_subdirectory(clustering/cast)
#add_subdirectory(clustering/dbscan)
add_subdirectory(fltm)
add_subdirectory(clustering/louvain)
add_subdirectory(data_generation)
add_subdirectory(clustering/clustering_algo)
add_subdirectory(gwas)
add_subdirectory(clustering_comparison)
