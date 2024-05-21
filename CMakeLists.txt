cmake_minimum_required(VERSION 3.16)

project(rbmd-gpu)

set(src 
./main.cpp
./src/reader/mmap_reader.cpp
./src/reader/atomic_reader.cpp
./src/reader/potential_reader.cpp
./src/reader/structure_reader.cpp)

include_directories("${CMAKE_CURRENT_LIST_DIR}/src")
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/model")

add_executable(${PROJECT_NAME} ${src})
