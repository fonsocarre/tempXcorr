CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(PIVlib)
# add directories to the end of the directory include paths
include_directories("include" ${CMAKE_CURRENT_SOURCE_DIR})

# for the stitchUtils lib
add_subdirectory("StitchUtils")
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()
SET( CMAKE_CXX_FLAGS -std=c++11 )

#target_link_libraries(stitchUtils hdf5 hdf5_cpp)
SET(SRC ${CMAKE_CURRENT_LIST_DIR}/src)
ADD_SUBDIRECTORY(${SRC})
file(GLOB PIVlib_SRC ".cpp")

#This tells CMake to fib.c and name it fibonacci
#add_library(PIVlib ${PIVlib_SRC})
