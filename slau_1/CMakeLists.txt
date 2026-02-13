cmake_minimum_required(VERSION 3.10)
project(tridiagonal_solver)

set(CMAKE_CXX_STANDARD 17)

add_library(solver src/solver.cpp)
target_include_directories(solver PUBLIC include)

add_executable(run_tests tests/main.cpp)
target_link_libraries(run_tests solver)
