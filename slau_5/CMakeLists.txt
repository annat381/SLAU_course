cmake_minimum_required(VERSION 3.10)
project(ChebyshevSolver VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(${CMAKE_SOURCE_DIR}/include)


include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)


add_executable(test_chebyshev tests/test_chebyshev.cpp)
target_link_libraries(test_chebyshev PRIVATE GTest::gtest_main)


enable_testing()
add_test(NAME ChebyshevTests COMMAND test_chebyshev)
