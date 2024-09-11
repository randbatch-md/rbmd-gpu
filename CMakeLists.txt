cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(rbmd-gpu)

#cmake build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "cmake build type." FORCE)
	message(STATUS "Set CMAKE_BUILD_TYPE \"Release\".")
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cmake)
get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)

#option
include("${cmake_dir}/option.cmake")

#hip setting
include("${cmake_dir}/hip_setting.cmake")

#sub directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tools)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src)
include_directories(${CMAKE_CURRENT_LIST_DIR}/src/common)


set(src
./src/main.cpp
./src/executioner/executioner.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
./src/parser/command_line.cpp
./src/common/unit_factor.cpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/executioner")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/common")

file(GLOB_RECURSE all_include_dirs
    "${CMAKE_CURRENT_LIST_DIR}/src/*.h"
)

set(include_dirs)
foreach(header_file ${all_include_dirs})
    get_filename_component(dir ${header_file} DIRECTORY)
    list(APPEND include_dirs ${dir})
endforeach()
list(REMOVE_DUPLICATES include_dirs)

include_directories(${include_dirs})

# 链接找到的所有静态库到可执行文件
file(GLOB tools_path ${CMAKE_CURRENT_LIST_DIR}/tools/*.tar.gz)

foreach(tool_path ${tools_path})
	get_filename_component(toolname ${tool_path} NAME)
	get_filename_component(tooldir ${tool_path} PATH)
	string(REGEX REPLACE "\\.tar\\.gz$" "" toolname_without_ext ${toolname})
	set(TOOL_INSTALL_PATH ${tooldir}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/${toolname_without_ext})
	list(APPEND CMAKE_PREFIX_PATH ${TOOL_INSTALL_PATH})
endforeach()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED jsoncpp IMPORTED_TARGET GLOBAL)
pkg_check_modules(SPDLOG REQUIRED spdlog IMPORTED_TARGET GLOBAL)
pkg_check_modules(CXXOPTS REQUIRED cxxopts IMPORTED_TARGET GLOBAL)






# 添加可执行文件
add_executable(${PROJECT_NAME} ${src})

# 设置库文件的路径
set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/lib)

# 定义要查找的库名列表
set(LIBRARY_NAMES
		data_manager
		reader
		neighbor_list
		force
		simulate_pipeline
)

target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRARY_NAMES} PkgConfig::JSONCPP PkgConfig::SPDLOG PkgConfig::CXXOPTS)


#set(src
#./src/main.cpp
#./src/system/md_system.cpp
#./src/executioner/executioner.cpp
#./src/reader/mmap_reader.cpp
#./src/reader/atomic_reader.cpp
#./src/reader/potential_reader.cpp
#./src/reader/structure_reader.cpp
#./src/Utilities/string_util.cpp
#./src/application/application.cpp
#./src/application/md_application.cpp
##./src/op/near_force/direct_truncation/direct_truncation_op.cpp
##./src/op/near_force/direct_truncation/rocm/direct_truncation_op.hip.cu
#./src/parser/command_line.cpp)


#./src/ensemble/MDensemble.cpp
#./src/ensemble/NVEensemble.cpp
#./src/ensemble/NVPensemble.cpp
#./src/ensemble/NVTensemble.cpp
#./src/controller/position_controller/default_position_controller.cpp
#./src/controller/velocity_controller/default_velocity_controller.cpp
#./src/controller/shake_controller/shake_controller.cpp

#)

#include
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/reader")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/system")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/executioner")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/op")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/data_manager/include")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/ensemble")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/controller/position_controller")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/controller/velocity_controller")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/controller/shake_controller")


#tool include
#get_property(var GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
#foreach(tool_include_path ${var})
#	include_directories(${tool_include_path})
#endforeach()
#
#add_executable(${PROJECT_NAME} ${src})
#
##tool link
#get_property(var GLOBAL PROPERTY TOOLS_LINK_PATH)
#foreach(tool_link_path ${var})
#	target_link_directories(${PROJECT_NAME} PUBLIC ${tool_link_path})
#endforeach()
#
##target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/src)
#
#target_link_libraries(${PROJECT_NAME} PRIVATE 
#jsoncpp 
#spdlog
#hip::device) 
#
##set target property
#include("${cmake_dir}/compile_definitions.cmake")

#print property
#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS ${PROJECT_NAME} PROPERTIES COMPILE_DEFINITIONS)