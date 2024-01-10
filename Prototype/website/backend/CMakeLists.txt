cmake_minimum_required(VERSION 3.20)
project(backend)

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

add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/../../../val3dity/ build)


add_executable(backend src/main.cpp src/definitions.h src/input.cpp src/input.h src/output.cpp src/output.h src/buildings.cpp src/buildings.h src/geometry.cpp src/geometry.h
        src/mesh3.cpp src/mesh3.h
        src/mesh2.cpp
        src/mesh2.h
        src/domain.cpp
        src/domain.h
        src/meshtools.cpp
        src/meshtools.h
        src/RoI.cpp
        src/RoI.h
        src/refinementboxes.cpp
        src/refinementboxes.h
        src/refinementboxes.h
        src/refinementboxes.h
        src/parameters.cpp
        src/parameters.h)

target_link_libraries(backend LINK_PUBLIC val3dity)
