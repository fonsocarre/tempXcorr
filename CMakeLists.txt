#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(txcorr)

# add directories to the end of the directory include paths
include_directories("include" ${CMAKE_CURRENT_SOURCE_DIR} PIVlib/include)
# for hsf5 path
# Your-external "mylib"
#add_library( hdf5 SHARED IMPORTED )
# You can define two import-locations: one for debug and one for release.
#set_target_properties( hdf5 PROPERTIES IMPORTED_LOCATION /usr/local/lib )
FIND_LIBRARY(RESULT libhdf5.a PATHS /usr/local/lib)
# for the PIV lib
add_subdirectory("PIVlib")


# Boost options
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()
SET( CMAKE_CXX_FLAGS -std=c++11 )
#set(BOOST_ROOT "/usr/local/boost")
#find_package( Boost lexical_cast property_tree)
#include_directories( ${Boost_INCLUDE_DIRS} )

target_link_libraries(PIVlib)
add_subdirectory(src)

set_target_properties( txcorr PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
