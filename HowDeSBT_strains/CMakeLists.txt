project(HowDeSBT)

cmake_minimum_required(VERSION 3.10.0)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(ZLIB)

if (MAKE_RELEASE AND NOT APPLE)
    set(STATIC_BUILD "-static")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${STATIC_BUILD}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")

include(ExternalProject)
ExternalProject_add(CROAR
    PREFIX CROAR
    GIT_REPOSITORY https://github.com/RoaringBitmap/CRoaring
    CMAKE_ARGS -DROARING_BUILD_STATIC=ON
    INSTALL_COMMAND "")

ExternalProject_add(SDSL
    PREFIX SDSL
    GIT_REPOSITORY https://github.com/simongog/sdsl-lite
    INSTALL_COMMAND "")

ExternalProject_add(ZLIB
    PREFIX ZLIB
    GIT_REPOSITORY https://github.com/madler/zlib
    INSTALL_COMMAND "")

If (APPLE)
    set(DYN_LIB "*.dylib")
else()
    set(DYN_LIB "*.so*")
endif()

add_custom_target(rm_dyn
    COMMAND rm ${PROJECT_BINARY_DIR}/ZLIB/src/ZLIB-build/${DYN_LIB}
    DEPENDS ZLIB
)

if (!STRAIN_ONLY)
ExternalProject_add(JELLY
    PREFIX JELLY
    URL https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/JELLY/src/JELLY/configure --prefix=${PROJECT_BINARY_DIR}/jelly_build 
    BUILD_COMMAND make
    INSTALL_COMMAND make install)


set(JELLY_BUILD_DIR ${PROJECT_BINARY_DIR}/jelly_build)
set(JELLY_INCLUDE ${JELLY_BUILD_DIR}/include/jellyfish-2.3.0/)
set(JELLY_LINK ${JELLY_BUILD_DIR}/lib)
endif()

set(CROAR_DIR ${PROJECT_BINARY_DIR}/CROAR/src/CROAR)
set(ZLIB_DIR ${PROJECT_BINARY_DIR}/ZLIB/src/ZLIB)
set(CROAR_BUILD_DIR ${PROJECT_BINARY_DIR}/CROAR/src/CROAR-build/src)
set(SDSL_BUILD_DIR ${PROJECT_BINARY_DIR}/SDSL/src/SDSL-build)
set(ZLIB_BUILD ${PROJECT_BINARY_DIR}/ZLIB/src/ZLIB-build)

set(CROAR_INCLUDE ${CROAR_DIR}/include/)
set(SDSL_INCLUDE ${SDSL_BUILD_DIR}/include)
set(ZLIB_INCLUDE ${ZLIB_DIR} ${ZLIB_BUILD})

set(CROAR_LINK ${CROAR_BUILD_DIR})
set(SDSL_LINK ${SDSL_BUILD_DIR}/lib)
set(ZLIB_LINK ${ZLIB_BUILD})

include_directories(${CROAR_INCLUDE} ${SDSL_INCLUDE} ${JELLY_INCLUDE} ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/thirdparty ${ZLIB_INCLUDE})

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

if (STRAIN_ONLY)
    link_directories(${CROAR_LINK} ${SDSL_LINK} ${ZLIB_LINK})
    link_libraries(roaring sdsl pthread z)
else()
    link_directories(${CROAR_LINK} ${SDSL_LINK} ${JELLY_LINK} ${ZLIB_LINK})
    link_libraries(roaring sdsl jellyfish-2.0 pthread z)
    set(CMAKE_INSTALL_RPATH ${CROAR_LINK} ${JELLY_LINK})
endif()

if (STRAIN_ONLY)
    set(PROJECT_FILES howdesbt.cc cmd_make_bf_qgram.cc cmd_query_qgram.cc cmd_cluster.cc cmd_build_sbt.cc
                      cmd_version.cc query.cc cmd_distance.cc bloom_tree.cc bloom_filter.cc bit_vector.cc
                      file_manager.cc bit_utilities.cc utilities.cc support.cc)
    add_definitions(-DonlyQgrams)             
elseif (CLASSICAL)
    set(PROJECT_FILES howdesbt.cc cmd_make_bf.cc cmd_cluster.cc cmd_build_sbt.cc cmd_query.cc 
                      cmd_version.cc query.cc bloom_tree.cc bloom_filter.cc bit_vector.cc 
                      file_manager.cc bit_utilities.cc utilities.cc support.cc)
    add_definitions(-DincludeJellyQuery)
elseif (BOTH)
    set(PROJECT_FILES howdesbt.cc cmd_make_bf.cc cmd_make_bf_qgram.cc cmd_query_qgram.cc cmd_cluster.cc
                      cmd_query.cc cmd_build_sbt.cc cmd_version.cc query.cc bloom_tree.cc bloom_filter.cc
                      bit_vector.cc file_manager.cc bit_utilities.cc utilities.cc support.cc)
    add_definitions(-DincludeJellyQuery -DonlyQgrams)
elseif (FULL)
    file (GLOB_RECURSE PROJECT_FILES *.cc)
    add_definitions(-DincludeAll -DonlyQgrams -DincludeJellyQuery)
endif()


add_executable(howdesbt ${PROJECT_FILES})
if (STRAIN_ONLY)
    add_dependencies(howdesbt CROAR SDSL ZLIB rm_dyn)
else()
    add_dependencies(howdesbt CROAR SDSL JELLY)
endif()




