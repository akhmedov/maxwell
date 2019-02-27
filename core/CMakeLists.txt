#
#  CMakeLists.txt
#  core.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "core")

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++14 -pipe -Wall -Wextra -Wpedantic -Wshadow \
    -Wundef -Wno-missing-field-initializers -fPIC"
)

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
get_filename_component(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(GNUMP_DIR ${ROOT_DIR}/gnump)

add_library(${PROJECT_NAME} SHARED
    ${SOURCE_DIR}/abstract_field.cpp
    ${SOURCE_DIR}/function.cpp
    ${SOURCE_DIR}/integral.cpp
    ${SOURCE_DIR}/linear_duhamel.cpp
    ${SOURCE_DIR}/linear_field.cpp
    ${SOURCE_DIR}/linear_medium.cpp
    ${SOURCE_DIR}/linear_source.cpp
    ${SOURCE_DIR}/logger.cpp
    ${SOURCE_DIR}/mysql_connect.cpp
    ${SOURCE_DIR}/noise.cpp
    ${SOURCE_DIR}/nonlinear_field.cpp
    ${SOURCE_DIR}/nonlinear_medium.cpp
    ${SOURCE_DIR}/phys_math.cpp
)

add_dependencies(${PROJECT_NAME} generate_build_version)

execute_process(COMMAND bash -c "mysql_config --variable=pkgincludedir" OUTPUT_VARIABLE MYSQL_INC_DIR)
target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${GNUMP_DIR}/include>
    $<BUILD_INTERFACE:${MYSQL_INC_DIR}>
)

find_library(GMP_LIB gmp HINTS ${GNUMP_DIR}/lib)
execute_process(COMMAND bash -c "mysql_config --variable=pkglibdir" OUTPUT_VARIABLE MYSQL_LIB_DIR)
string(STRIP ${MYSQL_LIB_DIR} MYSQL_LIB_DIR) # remove widespace in the end of MYSQL_LIB_DIR variable
find_library(MYSQL_LIB mysqlclient HINTS ${MYSQL_LIB_DIR})

target_link_libraries(${PROJECT_NAME} PUBLIC
    ${GMP_LIB}
    ${MYSQL_LIB}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${VERSION_MAJOR}
#    PUBLIC_HEADER "${PUBLIC_HEADERS}"
)

# include(GNUInstallDirs)

# install(TARGETS ${PROJECT_NAME}
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
# )

add_subdirectory(test)