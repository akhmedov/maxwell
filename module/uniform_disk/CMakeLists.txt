#
#  CMakeLists.txt
#  uniform_disk.module.maxwell
#
#  Created by Rolan Akhmedov on 27.02.19
#  Copyright © 2019 Rolan Akhmedov. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "module:uniform_disk")

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

add_dependencies(${PROJECT_NAME} generate_build_version core)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++14 -pipe -Wall -Wextra -Wpedantic -Wshadow \
    -Wundef -Wno-missing-field-initializers -fPIC"
)

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(GNUMP_DIR ${ROOT_DIR}/gnump)
set(CORE_DIR ${ROOT_DIR}/core)

add_library (${PROJECT_NAME} SHARED
    ${SOURCE_DIR}/kerr_amendment.cpp
    ${SOURCE_DIR}/uniform_disk_current.cpp
    ${SOURCE_DIR}/updisk_meandr.cpp
)

execute_process(COMMAND bash -c "mysql_config --variable=pkgincludedir" OUTPUT_VARIABLE MYSQL_INC_DIR)
target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${GNUMP_DIR}/include>
    $<BUILD_INTERFACE:${CORE_DIR}/include>
    $<BUILD_INTERFACE:${MYSQL_INC_DIR}>
)

find_library(CORE_LIB core HINTS ${ROOT_DIR}/build/core)
find_library(GMP_LIB gmp HINTS ${GNUMP_DIR}/lib)
execute_process(COMMAND bash -c "mysql_config --variable=pkglibdir" OUTPUT_VARIABLE MYSQL_LIB_DIR)
string(STRIP ${MYSQL_LIB_DIR} MYSQL_LIB_DIR) # remove widespace in the end of MYSQL_LIB_DIR variable
find_library(MYSQL_LIB mysqlclient HINTS ${MYSQL_LIB_DIR})

target_link_libraries(${PROJECT_NAME} PUBLIC
    ${CORE_LIB}
    ${GMP_LIB}
    ${MYSQL_LIB}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${VERSION_MAJOR}
#    PUBLIC_HEADER "${PUBLIC_HEADERS}"
)


