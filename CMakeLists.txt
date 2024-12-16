cmake_minimum_required(VERSION 3.10)

project(GeoJiang)


# Sourcefiles
set(SOURCES
    src/mathj.cpp
    src/main.cpp
    src/molecule.cpp
    src/atomType.cpp
    src/forceField.cpp
    src/internal.cpp
    src/hessian.cpp
)

set(CLAPACK
    liblapack.a
    libblas.a
    libf2c.a
)

# Test Matrix class
enable_testing()
set(TEST1
    test/test_matrix.cpp
    src/mathj.cpp)

include_directories(include)
include_directories(src)
link_directories(libs)

# set(CMAKE_VERBOSE_MAKEFILE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE debug)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

set(CMAKE_CXX_STANDARD 20)

add_executable(GeoJiang ${SOURCES})
target_link_libraries(GeoJiang PRIVATE ${CLAPACK})

add_executable(test_matrix ${TEST1})
target_link_libraries(test_matrix PRIVATE ${CLAPACK})
add_test(NAME test_matrix COMMAND test_matrix)

set_target_properties(GeoJiang PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})


