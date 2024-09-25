project(neighbor_list)

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

find_package(hipcub REQUIRED)
get_property(TOOLS_INCLUDE_PATH GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
get_property(TOOLS_LINK_PATH GLOBAL PROPERTY TOOLS_LINK_PATH)
set(inclue_files ${TOOLS_INCLUDE_PATH}  "../data_manager/include"  ${include_dirs} ${spdlog_include} ${hipcub_INCLUDE_DIRS})
set(link_dirs ${TOOLS_LINK_PATH}  ${spdlog_link} ${hipcub_LIBRARIES} )
#set(link_name "hip::device" "spdlog" "hip::hipcub")
set(depend_names "hip::hipcub"  "data_manager")



cpp_library(
name ${PROJECT_NAME}
depends_include ${inclue_files}
depends_link_dir ${link_dirs}
depends_name ${depend_names}
)