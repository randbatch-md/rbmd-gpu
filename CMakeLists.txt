cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
include("${cmake_dir}/gpu_setting.cmake")

#sub directory
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tools)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src)
include_directories(${CMAKE_CURRENT_LIST_DIR}/src/common)

file(GLOB_RECURSE output_srcs ./src/output/*.cpp ./src/output/*.cu)
set(src
./src/main.cpp
#./src/executioner/executioner.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
./src/parser/command_line.cpp
./src/common/unit_factor.cpp
./src/common/RBEPSample.h
./src/common/erf_table.h
${output_srcs}
		src/force/include/maceload.h
		src/force/src/maceload.cpp
		src/data_manager/include/scheduler/mace_memory_scheduler.h
		src/data_manager/include/force_field/mace_force_field_data.h
		src/data_manager/src/mace_memory_scheduler.cpp
		src/data_manager/src/mace_memory_scheduler.cpp
		src/data_manager/src/mace_memory_scheduler.cpp
		src/neighbor_list/src/op/rocm/mace_neighbor_list_op.hip.cu
		src/neighbor_list/src/op/mace_neighbor_list_op.h
		src/neighbor_list/src/neighbor_list_builder/mace_neighbor_list_builder.cpp
		src/neighbor_list/include/neighbor_list_builder/mace_neighbor_list_builder.h)

include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
#include_directories("${CMAKE_CURRENT_LIST_DIR}/src/executioner")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/common")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/output")
if(USE_CUDA)
	set_source_files_properties(${src} PROPERTIES LANGUAGE CUDA)
endif()
#torch
#find_package(hipblas REQUIRED)
#find_package(hipfft REQUIRED)
#find_package(hiprand REQUIRED)
#find_package(hipsparse REQUIRED)
#find_package(hipsolver REQUIRED)
#add gflags and libtorch
#set(CMAKE_PREFIX_PATH "/home/amd/zkh/gflags" "/home/amd/libtorch" ${CMAKE_PREFIX_PATH})
#find_package(gflags REQUIRED)
#et(TORCH_DIR "/home/amd/libtorch")
#find_package(Torch REQUIRED)
#set(torch_INCLUDE_DIRS "${TORCH_DIR}/include" "${TORCH_DIR}/include/torch/csrc/api/include" "${TORCH_DIR}/include/torch/csrc/autograd" "${TORCH_DIR}/include/torch/csrc/cuda")
#set(Torch_LINK_DIRS "${TORCH_DIR}/lib")
#include_directories(${torch_INCLUDE_DIRS})
#link_directories(${torch_LINK_DIRS})
#list(APPEND depend_names ${TORCH_LIBRARIES})
#list(APPEND depend_names PkgConfig::JSONCPP)
#list(APPEND include_dirs "${TORCH_DIR}/include")
#list(APPEND include_dirs "${TORCH_DIR}/include/torch/csrc/api/include")
#list(APPEND include_dirs "${TORCH_DIR}/include/torch/csrc/autograd")
#list(APPEND include_dirs "${TORCH_DIR}/include/torch/csrc/cuda")


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

#
file(GLOB tools_path ${CMAKE_CURRENT_LIST_DIR}/tools/*.tar.gz)

foreach(tool_path ${tools_path})
	get_filename_component(toolname ${tool_path} NAME)
	get_filename_component(tooldir ${tool_path} PATH)
	string(REGEX REPLACE "\\.tar\\.gz$" "" toolname_without_ext ${toolname})
	set(TOOL_INSTALL_PATH ${tooldir}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/${toolname_without_ext})
	list(APPEND CMAKE_PREFIX_PATH ${TOOL_INSTALL_PATH})
endforeach()

find_package(PkgConfig REQUIRED)
pkg_check_modules(SPDLOG REQUIRED spdlog IMPORTED_TARGET GLOBAL)
pkg_check_modules(CXXOPTS REQUIRED cxxopts IMPORTED_TARGET GLOBAL)


add_executable(${PROJECT_NAME} ${src})

set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/lib)

set(LIBRARY_NAMES
		data_manager
		reader
		neighbor_list
		force
		simulate_pipeline
)
if(USE_CUDA)
	target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRARY_NAMES}  PkgConfig::SPDLOG PkgConfig::CXXOPTS)