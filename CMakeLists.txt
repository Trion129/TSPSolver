cmake_minimum_required(VERSION 3.16)
project(TSPSolver)

set(CMAKE_CXX_STANDARD 17)

add_executable(TSPSolver lkh-solver.cpp tour-tree.hpp)
set(CGAL_ROOT /cygdrive/c/dev/CGAL-5.0.2)

find_package(Boost 1.60.0)
find_package(CGAL CONFIG REQUIRED)
target_link_libraries(TSPSolver PRIVATE CGAL::CGAL)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()