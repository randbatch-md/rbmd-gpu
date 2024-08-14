cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(controller)

file(GLOB_RECURSE force_controller_headers "${CMAKE_CURRENT_LIST_DIR}/force_controller/*.h" )
file(GLOB_RECURSE force_controller_sources "${CMAKE_CURRENT_LIST_DIR}/force_controller/*.cpp")
file(GLOB_RECURSE position_controller_headers "${CMAKE_CURRENT_LIST_DIR}/position_controller/*.h" )
file(GLOB_RECURSE position_controller_sources "${CMAKE_CURRENT_LIST_DIR}/position_controller/*.cpp")
file(GLOB_RECURSE temperature_controller_headers "${CMAKE_CURRENT_LIST_DIR}/temperature_controller/*.h" )
file(GLOB_RECURSE temperature_controller_sources "${CMAKE_CURRENT_LIST_DIR}/temperature_controller/*.cpp")
file(GLOB_RECURSE velocity_controller_headers "${CMAKE_CURRENT_LIST_DIR}/velocity_controller/*.h" )
file(GLOB_RECURSE velocity_controller_sources "${CMAKE_CURRENT_LIST_DIR}/velocity_controller/*.cpp")

add_library(${PROJECT_NAME} SHARED ${force_controller_headers} ${force_controller_sources})
add_library(${PROJECT_NAME} SHARED ${position_controller_headers} ${position_controller_sources})
add_library(${PROJECT_NAME} SHARED ${temperature_controller_headers} ${temperature_controller_sources})
add_library(${PROJECT_NAME} SHARED ${velocity_controller_headers} ${velocity_controller_sources})

target_include_directories( ${PROJECT_NAME} PUBLIC
	"${CMAKE_CURRENT_LIST_DIR}/force_controller"
	"${CMAKE_CURRENT_LIST_DIR}/position_controller"
	"${CMAKE_CURRENT_LIST_DIR}/temperature_controller"
	"${CMAKE_CURRENT_LIST_DIR}/velocity_controller"
)


