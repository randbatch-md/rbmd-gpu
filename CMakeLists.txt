cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

project(rbmd-gpu)

#option
include("./cmake/option.cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "cmake build type." FORCE)
	message(STATUS "not defined CMAKE_BUILD_TYPE.")
endif()

if(USE_ROCM)
  if(NOT DEFINED ROCM_PATH )
    set (ROCM_PATH "/opt/rocm"  CACHE STRING "Default ROCM installation directory." FORCE )
  endif ()
  message(STATUS "ROCM_PATH: ${ROCM_PATH}")

  if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
        set(HIP_PATH "${ROCM_PATH}/hip" CACHE PATH "Path to which HIP has been installed" FORCE)
    else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed" FORCE)
    endif()
  endif()
  message(STATUS "HIP_PATH: ${HIP_PATH}")

  set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
  set (HIP_HIPCC_FLAGS -fno-gpu-rdc; --std=c++14) # --amdgpu-target=gfx906

  # find_package(ROCM REQUIRED)
  find_package(HIP REQUIRED)

  if(HIP_FOUND)
    message(STATUS "Found HIP: " ${HIP_VERSION})
    message(STATUS "HIP include directories: ${HIP_INCLUDE_DIRS}")
	message(STATUS "HIP libraries: ${HIP_LIBRARIES}")
  else()
    message(FATAL_ERROR "Could not find HIP. Ensure that HIP is either installed in ${ROCM_PATH}/hip or the variable HIP_PATH is set to point to the right location.")
  endif()

  include_directories(
    ${ROCM_PATH}/include
  )

  target_link_libraries(${PROJECT_NAME}
    hip::host
    hip::device
  )

  add_compile_definitions(__ROCM)
  add_compile_definitions(__UT_USE_ROCM)
  add_compile_definitions(__HIP_PLATFORM_HCC__)
endif()

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
./src/op/near_force/direct_truncation/direct_truncation_op.cpp
./src/parser/command_line.cpp)

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

target_link_libraries(${PROJECT_NAME} PRIVATE jsoncpp spdlog device) 

#set target property
include("./cmake/set_target_property.cmake")

#print property
#include(CMakePrintHelpers)
#cmake_print_properties(TARGETS ${PROJECT_NAME} PROPERTIES COMPILE_DEFINITIONS)