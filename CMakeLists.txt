cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(rbmd-gpu)

#cmake build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "cmake build type." FORCE)
	message(STATUS "Set CMAKE_BUILD_TYPE \"Release\".")
endif()

#option
include("${CMAKE_CURRENT_LIST_DIR}/cmake/option.cmake")

#hip setting
include("${CMAKE_CURRENT_LIST_DIR}/cmake/hip_setting.cmake")

#sub directory
add_subdirectory(./src)
add_subdirectory(./tools)

set(src 
./src/main.cpp
./src/system/md_system.cpp
./src/executioner/executioner.cpp
./src/reader/mmap_reader.cpp
./src/reader/atomic_reader.cpp
./src/reader/potential_reader.cpp
./src/reader/structure_reader.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
#./src/op/near_force/direct_truncation/direct_truncation_op.cpp
#./src/op/near_force/direct_truncation/rocm/direct_truncation_op.hip.cu
./src/parser/command_line.cpp
./src/ensemble/NVEensemble.cpp
./src/ensemble/NVPensemble.cpp
./src/ensemble/NVTensemble.cpp
./src/controller/position_controller/default_position_controller.cpp
./src/controller/velocity_controller/default_velocity_controller.cpp)

#include
include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/reader")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/model")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/system")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/executioner")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/op")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/ensemble")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/controller/position_controller")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/controller/velocity_controller")


#tool include
get_property(var GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
foreach(tool_include_path ${var})
	include_directories(${tool_include_path})
endforeach()

add_executable(${PROJECT_NAME} ${src})

#tool link
get_property(var GLOBAL PROPERTY TOOLS_LINK_PATH)
foreach(tool_link_path ${var})
	target_link_directories(${PROJECT_NAME} PUBLIC ${tool_link_path})
endforeach()

#target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/src)

target_link_libraries(${PROJECT_NAME} PRIVATE 
jsoncpp 
spdlog
hip::device
device_rocm) 

#set target property
include("./cmake/set_target_property.cmake")

#print property
#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS ${PROJECT_NAME} PROPERTIES COMPILE_DEFINITIONS)