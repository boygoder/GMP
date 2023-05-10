cmake_minimum_required(VERSION 3.20.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
project(GMP CXX)
find_package(GMP REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(test_gmp test_gmp.cpp)
target_link_libraries(test_gmp PRIVATE gmpxx gmp)

add_executable(test_Newton test_Newton.cpp)
target_link_libraries(test_Newton PRIVATE gmpxx gmp)

add_executable(test_OrthoPoly test_OrthoPoly.cpp)
target_link_libraries(test_OrthoPoly PRIVATE gmpxx gmp)

add_executable(test_Prufer test_Prufer.cpp)
target_link_libraries(test_Prufer PRIVATE gmpxx gmp)

add_executable(test_RungeKutta test_RungeKutta.cpp)
target_link_libraries(test_RungeKutta PRIVATE gmpxx gmp)

add_executable(test_Gauss test_Gauss.cpp)
target_link_libraries(test_Gauss PRIVATE gmpxx gmp)


