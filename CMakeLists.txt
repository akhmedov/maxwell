#
#  CMakeLists.txt
#  maxwell
#
#  Created by Oleh Zahrychanskyi on 27.02.19
#  Copyright © 2019 Oleh Zahrychanskyi. All rights reserved.
#

cmake_minimum_required(VERSION 3.0)

# Versioning
set(VERSION_MAJOR 0)
set(VERSION_MINOR 9)
set(VERSION_PATCH 0)

# Create build version on every build from the git HEAD commit hash
add_custom_command(
    OUTPUT build_version
    COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}
        -DVERSION_BUILD:INTERNAL=`git rev-parse --short HEAD`
)

add_custom_target(generate_build_version ALL
    DEPENDS build_version
)

# Build type
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_subdirectory(core)
add_subdirectory(module/uniform_disk)
