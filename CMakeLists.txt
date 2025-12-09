cmake_minimum_required(VERSION 3.16)

execute_process(
    COMMAND getconf GNU_LIBC_VERSION
    OUTPUT_VARIABLE GLIBC_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
string(REGEX MATCH "[0-9]+\\.[0-9]+" GLIBC_VERSION_NUM "${GLIBC_VERSION}")
message(STATUS "Found GLIBC version: ${GLIBC_VERSION_NUM}")
    
    # GLIBC 2.19 及以上通常支持新的 CXX11 ABI
if(GLIBC_VERSION_NUM VERSION_GREATER_EQUAL "2.19")
        add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
        message(STATUS "Using new CXX11 ABI (GLIBC >= 2.19)")
else()
        add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
        message(STATUS "Using old ABI (GLIBC < 2.19)")
endif()

set_property(DIRECTORY . PROPERTY USE_MACE ${USE_MACE})
 
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
message(STATUS "Global USE_MACE setting: ${USE_MACE}")

file(GLOB_RECURSE output_srcs ./src/output/*.cpp ./src/output/*.cu)
set(add_files)
if(USE_MACE)
list(APPEND add_files
         src/force/include/maceload.h
                src/force/src/maceload.cpp
                src/data_manager/include/scheduler/mace_memory_scheduler.h
                src/data_manager/include/force_field/mace_force_field_data.h
                src/data_manager/src/mace_memory_scheduler.cpp
                src/neighbor_list/src/op/rocm/mace_neighbor_list_op.hip.cu
                src/neighbor_list/src/op/mace_neighbor_list_op.h
                src/neighbor_list/src/neighbor_list_builder/mace_neighbor_list_builder.cpp
                src/neighbor_list/include/neighbor_list_builder/mace_neighbor_list_builder.h
           )
endif()
set(src
./src/main.cpp
#./src/executioner/executioner.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
./src/parser/command_line.cpp
./src/common/unit_factor.cpp
./src/common/thermo_tui_view.cpp
./src/common/RBEPSample.h
./src/common/erf_table.h
${output_srcs}
${add_files})

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


file(GLOB_RECURSE all_include_dirs
    "${CMAKE_CURRENT_LIST_DIR}/src/*.h"
)

set(include_dirs)
foreach(header_file ${all_include_dirs})
    get_filename_component(dir ${header_file} DIRECTORY)
    list(APPEND include_dirs ${dir})
endforeach()
list(REMOVE_DUPLICATES include_dirs)
if(USE_MACE)
include_directories(${include_dirs}
      "${TORCH_DIR}/include"
      "${TORCH_DIR}/include/torch/csrc/api/include"
      "${glog_DIR}/include" )
else()
include_directories(${include_dirs})
endif()
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
pkg_check_modules(FTXUI REQUIRED ftxui IMPORTED_TARGET GLOBAL)


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
if(USE_MACE)
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MACE=1)
else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MACE=0)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRARY_NAMES}  PkgConfig::SPDLOG PkgConfig::CXXOPTS PkgConfig::FTXUI)
