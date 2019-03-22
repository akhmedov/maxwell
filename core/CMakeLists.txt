#
#  CMakeLists.txt
#  core.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME maxwell)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SOURCES
    ${SOURCE_DIR}/abstract_field.cpp
    ${SOURCE_DIR}/dataset.cpp
    ${SOURCE_DIR}/function.cpp
    ${SOURCE_DIR}/gnu_plot.cpp
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
    ${SOURCE_DIR}/module_manager.cpp
    ${SOURCE_DIR}/space_point.cpp
)

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

add_library(${PROJECT_NAME} SHARED ${SOURCES})
# add_dependencies(${PROJECT_NAME} generate_build_version)
target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PUBLIC ${GNUMP_INCLUDE_DIR}
    PUBLIC ${MYSQL_INCLUDE_DIR}
)
target_link_libraries(${PROJECT_NAME}
    PUBLIC ${GNUMP_LIB}
    PRIVATE ${MYSQL_LIB}
)
if (${UNIX})
    target_link_libraries(${PROJECT_NAME} INTERFACE pthread)
endif()
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

enable_testing()
add_subdirectory(test)
