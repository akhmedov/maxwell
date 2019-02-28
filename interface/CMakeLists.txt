#
#  CMakeLists.txt
#  core.maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME interface)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SOURCES
    ${SOURCE_DIR}/dataset.cpp
    ${SOURCE_DIR}/gnu_plot.cpp
    ${SOURCE_DIR}/main.cpp
)

project(${PROJECT_NAME}
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES CXX
)

add_executable(${PROJECT_NAME} ${SOURCES})
add_dependencies(${PROJECT_NAME} generate_build_version)
target_include_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(${PROJECT_NAME} PRIVATE
    maxwell
    uniform_disk
)

# include(GNUInstallDirs)

# install(TARGETS ${PROJECT_NAME}
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
# )

#enable_testing()
#add_subdirectory(test)
