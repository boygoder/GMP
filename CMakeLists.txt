cmake_minimum_required(VERSION 3.20.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
project(GMP CXX)
find_package(GMP REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# 构建静态库 g++ -c *.cpp
add_library(gmptools STATIC gmptools.h gmptools.cpp)
target_link_libraries(gmptools PUBLIC gmpxx gmp)

add_library(polynomial STATIC polynomial.h polynomial.cpp)
target_link_libraries(polynomial PUBLIC gmptools)

add_library(OrthogonalPolynomials STATIC OrthogonalPolynomials.h OrthogonalPolynomials.cpp)
target_link_libraries(OrthogonalPolynomials PUBLIC polynomial)

add_library(PruferTransformer STATIC PruferTransformer.h PruferTransformer.cpp)
target_link_libraries(PruferTransformer PUBLIC polynomial)

add_library(RungeKutta STATIC RungeKutta.h RungeKutta.cpp)
target_link_libraries(RungeKutta PUBLIC polynomial)

add_library(GaussianPoint STATIC GaussianPoint.h GaussianPoint.cpp)
target_link_libraries(GaussianPoint PUBLIC gmptools)

add_library(NewtonMethod STATIC NewtonMethod.h NewtonMethod.cpp)
target_link_libraries(NewtonMethod PUBLIC polynomial)

add_library(GaussIntegralTableGenerator STATIC GaussIntegralTableGenerator.h GaussIntegralTableGenerator.cpp)
target_link_libraries(GaussIntegralTableGenerator PUBLIC GaussianPoint OrthogonalPolynomials PruferTransformer RungeKutta NewtonMethod)

# 生成可执行文件 g++ -std=c++17 main.cpp -o main -lgmpxx -lgmp
add_executable(test_gmp test_gmp.cpp)
target_link_libraries(test_gmp OrthogonalPolynomials)

add_executable(test_polynomial test_polynomial.cpp )
target_link_libraries(test_polynomial polynomial)

add_executable(test_Newton test_Newton.cpp)
target_link_libraries(test_Newton OrthogonalPolynomials RungeKutta PruferTransformer NewtonMethod)

add_executable(test_OrthoPoly test_OrthoPoly.cpp)
target_link_libraries(test_OrthoPoly OrthogonalPolynomials)

add_executable(test_Prufer test_Prufer.cpp)
target_link_libraries(test_Prufer OrthogonalPolynomials PruferTransformer)

add_executable(test_RungeKutta test_RungeKutta.cpp)
target_link_libraries(test_RungeKutta OrthogonalPolynomials RungeKutta PruferTransformer)

add_executable(test_Gauss test_Gauss.cpp)
target_link_libraries(test_Gauss GaussIntegralTableGenerator)

add_executable(main main.cpp)
target_link_libraries(main GaussIntegralTableGenerator)
