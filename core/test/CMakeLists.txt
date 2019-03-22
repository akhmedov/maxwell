#
#  CMakeLists.txt
#  test.core.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "core-test")

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES
    ${SOURCE_DIR}/bessel_perp.cpp
    ${SOURCE_DIR}/binom_prod.cpp
    ${SOURCE_DIR}/dataset_metadata.cpp
    ${SOURCE_DIR}/dataset_series.cpp
    ${SOURCE_DIR}/invers_sqrt.cpp
    ${SOURCE_DIR}/monte_carlo_improper.cpp
    ${SOURCE_DIR}/monte_carlo_integral.cpp
    ${SOURCE_DIR}/monte_carlo_vector.cpp
    ${SOURCE_DIR}/mysql_driver.cpp
    ${SOURCE_DIR}/next_prime.cpp
    ${SOURCE_DIR}/point_converter.cpp
    ${SOURCE_DIR}/simpson2d.cpp
    ${SOURCE_DIR}/simpson2d_line2.cpp
    ${SOURCE_DIR}/simpson2d_line.cpp
    ${SOURCE_DIR}/simpson3d.cpp
    ${SOURCE_DIR}/simpson_runge_2d.cpp
    ${SOURCE_DIR}/simpson_runge.cpp
    ${SOURCE_DIR}/vector_to_matrix.cpp
    ${SOURCE_DIR}/yacobi_pol_compare.cpp
    ${SOURCE_DIR}/yacobi_pol_property.cpp
)

project(${PROJECT_NAME} LANGUAGES CXX)

link_libraries(
    maxwell
)

foreach(SOURCE ${SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE} NAME)
    string(REPLACE ".cpp" "" TEST ${SOURCE_NAME})

    add_executable(${TEST} ${SOURCE})
    add_test(${TEST} ${TEST})
endforeach()
