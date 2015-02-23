set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

#set(PROJECT_DIR ..)
set(LIB_DIR "${PROJECT_SOURCE_DIR}/dependencies/" )
set(Boost_INCLUDE_DIR "${LIB_DIR}/boost_1_55_0/target/include/")
set(Boost_LIBRARY_DIR "${LIB_DIR}/boost_1_55_0/target/lib/")



set(ProBT_INCLUDE_DIR "${LIB_DIR}/probt_spl_2_3_99/include/")
set(ProBT_LIBRARY_DIR "${LIB_DIR}/probt_spl_2_3_99/lib/" )


MESSAGE(${ProBT_LIBRARY_DIR})

include_directories( 
  "${ProBT_INCLUDE_DIR}"
  "${Boost_INCLUDE_DIR}"
  "${PROJECT_SOURCE_DIR}/core/include" 
  "${PROJECT_SOURCE_DIR}/tests/common" 
)

link_directories(${Boost_LIBRARY_DIR} ${ProBT_LIBRARY_DIR})
file(GLOB TESTS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp  )

file(GLOB_RECURSE SOURCE "${PROJECT_SOURCE_DIR}/core/src/*.cpp" )

add_executable( test_em ${CMAKE_CURRENT_SOURCE_DIR}/test_runner.cpp ${CMAKE_CURRENT_SOURCE_DIR}/test_probayes.cpp ${SOURCE} )
target_link_libraries( test_em boost_program_options boost_filesystem boost_system boost_regex boost_unit_test_framework spl  )

