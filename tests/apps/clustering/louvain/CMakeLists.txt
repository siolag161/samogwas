

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)


set(LIB_DIR "${PROJECT_SOURCE_DIR}/dependencies/" )
set(Boost_INCLUDE_DIR "${LIB_DIR}/boost_1_55_0/target/include/")
set(Boost_LIBRARY_DIR "${LIB_DIR}/boost_1_55_0/target/lib/")

set(LOUVAIN_INCL "${PROJECT_SOURCE_DIR}/apps/include/" )
set(LOUVAIN_SRC "${PROJECT_SOURCE_DIR}/apps/src/clustering/louvain" )


include_directories( 
  "${ProBT_INCLUDE_DIR}"
  "${Boost_INCLUDE_DIR}"
  "${PROJECT_SOURCE_DIR}/core/include" 
  "${PROJECT_SOURCE_DIR}/tests/common" 
  "${LOUVAIN_INCL}"
)
MESSAGE(${LOUVAIN_SRC})

link_directories(${Boost_LIBRARY_DIR} ${ProBT_LIBRARY_DIR})
file(GLOB TESTS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp  )
file(GLOB LOUVAIN_SRC  ${LOUVAIN_SRC}/*.cpp  )

add_executable( test_louvain ${TESTS} ${LOUVAIN_SRC} )
target_link_libraries( test_louvain boost_program_options boost_filesystem boost_system boost_regex boost_unit_test_framework spl  )

