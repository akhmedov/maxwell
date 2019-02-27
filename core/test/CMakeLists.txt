#
#  CMakeLists.txt
#  test.core.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "core-test")

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++14 -pipe -Wall -Wextra -Wpedantic -Wshadow \
    -Wundef -Wno-missing-field-initializers -fPIC"
)

set(CORE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(GNUMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../gnump)

include_directories(
    ${CORE_DIR}/include
    ${GNUMP_DIR}/include
)
find_library(GMP_LIB gmp HINTS ${GNUMP_DIR}/lib)
link_libraries(
    core
    ${GMP_LIB}
    mysqlclient
)

file(GLOB sources RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")
foreach(source ${sources})
    string(REPLACE ".cpp" "" test ${source})

    add_executable(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${source})
endforeach()