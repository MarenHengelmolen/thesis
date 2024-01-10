cmake_minimum_required(VERSION 3.20)
project(histogram)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-O3" )


#CGAL
set (CGAL_DIR "-I/usr/include/CGAL")
set(CGAL_DATA_DIR ".")
find_package(CGAL QUIET COMPONENTS )
if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )
    message(STATUS ${CGAL_LIBRARIES})
    message(STATUS ${CGAL_3RD_PARTY_LIBRARIES})
else()
    message(SEND_ERROR "CGAL library is required")
    return()
endif()

find_package(Eigen3 3.4.0)
include(CGAL_Eigen3_support)

add_executable(histogram src/main.cpp src/definitions.h src/input.cpp src/input.h src/buildings.cpp src/buildings.h
        src/output.cpp
)



