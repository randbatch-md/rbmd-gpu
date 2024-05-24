cmake_minimum_required(VERSION 3.16)

project(rbmd-gpu)

#option
include("./cmake/option.cmake")

add_subdirectory(./src)

set(src 
./src/main.cpp
./src/system/md_system.cpp
./src/reader/mmap_reader.cpp
./src/reader/atomic_reader.cpp
./src/reader/potential_reader.cpp
./src/reader/structure_reader.cpp
./src/Utilities/string_util.cpp
./src/application/application.cpp
./src/application/md_application.cpp
./src/parser/command_line.cpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/reader")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/model")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/Utilities")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/parser")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/application")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/system")

add_executable(${PROJECT_NAME} ${src})

#target property
set_property(TARGET ${PROJECT_NAME} PROPERTY TVAR "tval")
get_property(val TARGET ${PROJECT_NAME} PROPERTY TVAR)
message("val = ${val}")

if(${ENABLE_DOUBLE})
	set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_DEFINITIONS "USE_DOUBLE=${ENABLE_DOUBLE}")
endif()

if(${ENABLE_64BIT_IDS})
	set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_DEFINITIONS "USE_LONGLONG=${ENABLE_64BIT_IDS}")
endif()