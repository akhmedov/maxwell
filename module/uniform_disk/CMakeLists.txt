#
#  CMakeLists.txt
#  uniform_disk.module.maxwell
#
#  Created by Rolan Akhmedov on 27.02.19
#  Copyright © 2019 Rolan Akhmedov. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME uniform_disk)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SOURCES
    ${SOURCE_DIR}/kerr_amendment.cpp
    ${SOURCE_DIR}/uniform_disk_current.cpp
    ${SOURCE_DIR}/updisk_meandr.cpp
)

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

add_library(${PROJECT_NAME} SHARED ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    maxwell
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${VERSION_MAJOR}
#    PUBLIC_HEADER "${PUBLIC_HEADERS}"
)

enable_testing()
add_subdirectory(test)