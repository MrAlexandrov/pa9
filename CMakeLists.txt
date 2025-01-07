cmake_minimum_required(VERSION 3.14)
project(main CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
    include
    third_party/matrix/include
)

add_executable(${PROJECT_NAME} main.cpp)

enable_testing()
add_subdirectory(third_party/matrix/tests)

set(CTEST_BINARY_DIRECTORY ${CMAKE_BINARY_DIR})
