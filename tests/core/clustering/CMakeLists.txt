set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)



include_directories( 
  "${ProBT_INCLUDE_DIR}"
  "${Boost_INCLUDE_DIR}"
  "${PROJECT_SOURCE_DIR}/core/include" 
  "${PROJECT_SOURCE_DIR}/tests/common" 
)

link_directories(${Boost_LIBRARY_DIR} ${ProBT_LIBRARY_DIR})
file(GLOB TESTS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp  )

add_executable( test_clustering ${TESTS}  )
target_link_libraries( test_clustering boost_program_options boost_filesystem boost_system boost_regex boost_unit_test_framework spl  )

