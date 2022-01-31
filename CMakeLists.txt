cmake_minimum_required(VERSION 3.5)
project(LeastSquares3DTransform CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/leastSqs3Dparams.yaml
				DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include(ExternalProject)

ExternalProject_Add(eigen-mvn
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen-mvn
    GIT_REPOSITORY https://github.com/jlong29/eigen-mvn.git
)

# Needed to find headers and library
ExternalProject_Get_Property(eigen-mvn install_dir)
message("eigen-mvn install directory is " ${install_dir})

add_library(eigen_mvn STATIC IMPORTED)

# Needed to find the library
# NOTE: /src/eigen-mvn-build are default settings for how cmake organizes downloands and builds
# same for /src/eigen-mvn for headers below
set_property(TARGET eigen_mvn PROPERTY IMPORTED_LOCATION ${install_dir}/src/eigen-mvn-build/libeigen_mvn.a)

# Needed to trigger build of library
add_dependencies(eigen_mvn eigen-mvn)

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
set_target_properties(leastSq3D1 PROPERTIES LINKER_LANGUAGE CXX)
target_compile_options(leastSq3D1 PRIVATE -O3 -Wall -Wextra)
target_include_directories(
	leastSq3D1
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
		$ENV{EIGEN_ROOT_DIR}
		${install_dir}/src/eigen-mvn
		${YAMLCPP_INCLUDE_DIRS}
)
target_link_libraries(leastSq3D1 utilities eigen_mvn ${YAMLCPP_LIBRARIES})
