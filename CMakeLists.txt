cmake_minimum_required(VERSION 3.5)
project(LeastSquares3DTransform CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/leastSqs3Dparams.yaml
				DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include(ExternalProject)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(eigen-mvn
    GIT_REPOSITORY https://github.com/jlong29/eigen-mvn.git
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-mvn
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/eigen-mvn-install}
)

ExternalProject_Get_Property(eigen-mvn install_dir)

message("eigen-mvn install directory is " ${install_dir})

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED)

# SYSTEM LIBRARIES
pkg_check_modules(YAMLCPP REQUIRED yaml-cpp)
if(${YAMLCPP_FOUND})
	 MESSAGE(STATUS "Found YAMLCPP")
	 MESSAGE(STATUS "YAMLCPP INCLUDES: ${YAMLCPP_INCLUDE_DIRS}")
	 MESSAGE(STATUS "YAMLCPP LIBRARIES: ${YAMLCPP_LIBRARIES}")
else()
	 MESSAGE(STATUS "Could not locate YAMLCPP.")
endif()

# LOCAL LIBRARIES
add_library(utilities
		utilities.h
		utilities.cpp
)

target_include_directories(
	utilities
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		$ENV{EIGEN_ROOT_DIR}
		${YAMLCPP_INCLUDE_DIRS}
)

target_link_libraries(utilities ${YAMLCPP_LIBRARIES})

# Case1: Basic Recreation of Python example
add_executable(leastSq3D1 leastSq3D1.cpp)
add_dependencies(leastSq3D1 eigen-mvn)
set_target_properties(leastSq3D1 PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(leastSq3D1 PRIVATE -O3 -Wall -Wextra)
target_include_directories(
	leastSq3D1
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		$ENV{EIGEN_ROOT_DIR}
		${YAMLCPP_INCLUDE_DIRS}
)
target_link_libraries(leastSq3D1 utilities eigen_mvn ${YAMLCPP_LIBRARIES})
