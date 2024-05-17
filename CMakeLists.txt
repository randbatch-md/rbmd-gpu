cmake_minimum_required(VERSION 3.16)

project(rbmd-gpu)

set(src 
./main.cpp
./src/reader/conf_reader.cpp)

add_executable(${PROJECT_NAME} ${src})
