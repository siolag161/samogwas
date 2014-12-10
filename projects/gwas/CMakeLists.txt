set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#set(PROJECT_DIR ..)
set(LIB_DIR "${PROJECT_SOURCE_DIR}/dependencies/" )
set(Boost_INCLUDE_DIR "${LIB_DIR}/boost_1_55_0/target/include/")
set(Boost_LIBRARY_DIR "${LIB_DIR}/boost_1_55_0/target/lib/")

set(ProBT_INCLUDE_DIR "${LIB_DIR}/probt_spl_2_3_0/include/")
set(ProBT_LIBRARY_DIR "${LIB_DIR}/probt_spl_2_3_0/lib/" )

set(LOUVAIN_INCL "${PROJECT_SOURCE_DIR}/apps/include/" )
set(LOUVAIN_SRC "${PROJECT_SOURCE_DIR}/apps/src/clustering/louvain" )

set(PROJECT_COMMON "${PROJECT_SOURCE_DIR}/projects/common" )

include_directories( 
  "${ProBT_INCLUDE_DIR}"
  "${Boost_INCLUDE_DIR}"
  "${PROJECT_SOURCE_DIR}/core/include" 
  "${PROJECT_COMMON}" 
  "${LOUVAIN_INCL}"
)

link_directories(${Boost_LIBRARY_DIR} ${ProBT_LIBRARY_DIR})

file(GLOB_RECURSE SOURCE "${PROJECT_SOURCE_DIR}/core/src/*.cpp" )
file(GLOB LOUVAIN_SRC  ${LOUVAIN_SRC}/*.cpp  )

add_executable( gwas_tools main.cpp gwas_common.cpp ${SOURCE}  ${LOUVAIN_SRC} )
target_link_libraries( gwas_tools boost_program_options boost_filesystem boost_system boost_regex spl  )



add_executable( gwas main_gwas.cpp gwas_common.cpp ${SOURCE}  ${LOUVAIN_SRC} )
target_link_libraries( gwas boost_program_options boost_filesystem boost_system boost_regex spl  )
