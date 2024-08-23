project(simulate_pipeline)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)

get_property(spdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(spdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)


file(GLOB_RECURSE all_include_dirs
    "${CMAKE_CURRENT_LIST_DIR}/include/*.h"
    "${CMAKE_CURRENT_LIST_DIR}/src/*.h"
)

set(include_dirs)
foreach(header_file ${all_include_dirs})
    get_filename_component(dir ${header_file} DIRECTORY)
    list(APPEND include_dirs ${dir})
endforeach()
list(REMOVE_DUPLICATES include_dirs)

set(inclue_files "../data_manager/include" "../op/base" ${include_dirs} ${spdlog_include})
set(link_dirs ${LIBRARY_DIR} ${spdlog_link} )
set(depend_names "data_manager" "controller")

cpp_library(
name ${PROJECT_NAME}
depends_include ${inclue_files}
depends_link_dir ${link_dirs}
depends_name ${depend_names}
)