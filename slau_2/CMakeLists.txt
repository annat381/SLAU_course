cmake_minimum_required(VERSION 3.10)
project(csr_matrix)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_executable(csr_matrix test/main.cpp)
