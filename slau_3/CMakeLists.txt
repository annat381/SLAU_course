cmake_minimum_required(VERSION 3.14)
project(QRsolve VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")


add_library(MATRIX STATIC
    src/Vector.cpp
    src/Matrix.cpp
    src/QR.cpp
)

target_include_directories(MATRIX PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable(run_tests
    tests/testQR.cpp
)

target_link_libraries(run_tests PRIVATE
    MATRIX
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(run_tests)
