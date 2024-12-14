set(CMAKE_CXX_STANDARD 20)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(matrix_test matrix_test.cpp)
target_link_libraries(matrix_test PRIVATE GTest::GTest GTest::Main)
add_test(NAME ConstructorAccessTest
         COMMAND matrix_test --gtest_color=yes --gtest_output=xml:${CMAKE_BINARY_DIR}/test-results/matrix_test.xml)

# TODO: add test on classes without some constructors, operators etc.