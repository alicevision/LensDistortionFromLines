CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(cameraLensCalibration C CXX)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
endif ()

#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FIND_PACKAGE(PNG REQUIRED)

file(
        GLOB_RECURSE
        files_header
        src/*.h
)
file(
        GLOB_RECURSE
        files_cpp
        src/*.cpp
)

add_executable(cameraLensCalibration ${files_header} ${files_cpp})
target_link_libraries(cameraLensCalibration
  ${PNG_LIBRARIES}
)
target_include_directories(cameraLensCalibration PUBLIC
  ${PNG_INCLUDE_DIRS}
)
