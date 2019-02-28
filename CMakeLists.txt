#
#  CMakeLists.txt
#  maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

# Variables
set(VERSION_MAJOR 0)
set(VERSION_MINOR 9)
set(VERSION_PATCH 0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wundef -Wno-missing-field-initializers")

# Build type
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, setting to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# GnuMP and MySQL libs
set(GNUMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gnump)
set(GNUMP_INCLUDE_DIR ${GNUMP_DIR}/include)
set(GNUMP_LIB_DIR ${GNUMP_DIR}/lib)
find_library(GNUMP_LIB gmp
    HINTS ${GNUMP_LIB_DIR}
    NO_DEFAULT_PATH
)

execute_process(COMMAND bash -c "mysql_config --variable=pkgincludedir" 
    OUTPUT_VARIABLE MYSQL_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(COMMAND bash -c "mysql_config --variable=pkglibdir" 
    OUTPUT_VARIABLE MYSQL_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_library(MYSQL_LIB mysqlclient
    HINTS ${MYSQL_LIB_DIR}
    NO_DEFAULT_PATH
)

# Create build version on every build from the git HEAD commit hash
add_custom_command(
    OUTPUT build_version
    COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}
        -DVERSION_BUILD:INTERNAL=`git rev-parse --short HEAD`
)

add_custom_target(generate_build_version ALL
    DEPENDS build_version
)

add_subdirectory(core)
add_subdirectory(module)
add_subdirectory(example)
add_subdirectory(interface)
