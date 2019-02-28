#
#  CMakeLists.txt
#  test.uniform_disk.module.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "uniform_disk-test")
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES
    ${SOURCE_DIR}/anal_numerical_I1.cpp
    ${SOURCE_DIR}/anal_numerical_I2.cpp
    ${SOURCE_DIR}/anal_numerical_I3.cpp
    ${SOURCE_DIR}/anal_numerical_I4.cpp
    ${SOURCE_DIR}/anal_numerical_I5.cpp
    ${SOURCE_DIR}/I1_time_partder.cpp
    ${SOURCE_DIR}/I2_time_partder.cpp
)

project(${PROJECT_NAME} LANGUAGES CXX)

link_libraries(
    uniform_disk
)

foreach(SOURCE ${SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE} NAME)
    string(REPLACE ".cpp" "" TEST ${SOURCE_NAME})

    add_executable(${TEST} ${SOURCE})
    add_test(${TEST} ${TEST})
endforeach()
