cmake_minimum_required(VERSION 3.10)
project(IterativeMethods VERSION 1.0)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


add_library(iterative_methods_lib
    src/vector.cpp
    src/csr_matrix.cpp
    src/iterative_methods.cpp
)


target_include_directories(iterative_methods_lib
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)


add_executable(tests tests/test_methods.cpp)


target_link_libraries(tests
    iterative_methods_lib
    GTest::gtest_main
)


include(GoogleTest)
gtest_discover_tests(tests)


