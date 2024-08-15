project(data_manager)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)
cpp_library(${PROJECT_NAME})

