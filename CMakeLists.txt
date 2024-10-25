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

file(GLOB_RECURSE output_srcs ./src/output/*.cpp ./src/output/*.cu)
set(src
./src/main.cpp
./src/executioner/executioner.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
./src/parser/command_line.cpp
./src/common/unit_factor.cpp
./src/common/RBEPSample.h
./src/common/erf_table.h
${output_srcs})

include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/executioner")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/common")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/output")


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

# �����ҵ������о�̬�⵽��ִ���ļ�
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






# ��ӿ�ִ���ļ�
add_executable(${PROJECT_NAME} ${src})

# ���ÿ��ļ���·��
set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/lib)

# ����Ҫ���ҵĿ����б�
set(LIBRARY_NAMES
		data_manager
		reader
		neighbor_list
		force
		simulate_pipeline
)

target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRARY_NAMES} PkgConfig::JSONCPP PkgConfig::SPDLOG PkgConfig::CXXOPTS)