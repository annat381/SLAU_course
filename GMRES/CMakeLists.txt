cmake_minimum_required(VERSION 3.10)
project(GMRES_Solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(include)

# Создаем исполняемый файл из тестов
add_executable(test_gmres tests/test_gmres.cpp)


