#
#  CMakeLists.txt
#  example.maxwell
#
#  Created by Oleh Zahrychanskyi on 28.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME example)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES
    # ${SOURCE_DIR}/plot_awgn_power.cpp
    ${SOURCE_DIR}/plot_emp_shape.cpp
    # ${SOURCE_DIR}/plot_emp_surface.cpp
    # ${SOURCE_DIR}/plot_kerr_source.cpp
    ${SOURCE_DIR}/plot_energy_slyse.cpp
    ${SOURCE_DIR}/plot_energy_distribution.cpp
    ${SOURCE_DIR}/dataset_samesnr.cpp
)

project(${PROJECT_NAME} LANGUAGES CXX)

link_libraries(
    maxwell
    uniform_disk
)

foreach(SOURCE ${SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE} NAME)
    string(REPLACE ".cpp" "" TEST ${SOURCE_NAME})

    add_executable(${TEST} ${SOURCE})
endforeach()

#enable_testing()
#add_subdirectory(test)
