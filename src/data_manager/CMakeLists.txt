project(data_manager)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)
get_property(sppdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(sppdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)

get_property(TOOLS_INCLUDE_PATH GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
get_property(TOOLS_LINK_PATH GLOBAL PROPERTY TOOLS_LINK_PATH)

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

#set(include_dir ${TOOLS_INCLUDE_PATH}  ${sppdlog_include}    ${rocthrust_INCLUDE_DIRS} )  # ${jsoncpp_INCLUDE_DIRS}
#set(link_dir  ${TOOLS_LINK_PATH}  ${sppdlog_link}  ${rocthrust_LIBRARIES}   )  # ${jsoncpp_LIBRARIES}
#set(link_name "jsoncpp_lib" "spdlog" "roc::rocthrust")



# 链接找到的所有静态库到可执行文件
file(GLOB tools_path ${CMAKE_CURRENT_LIST_DIR}/../../tools/*.tar.gz)

foreach(tool_path ${tools_path})
	get_filename_component(toolname ${tool_path} NAME)
	get_filename_component(tooldir ${tool_path} PATH)
	string(REGEX REPLACE "\\.tar\\.gz$" "" toolname_without_ext ${toolname})
	set(TOOL_INSTALL_PATH ${tooldir}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/${toolname_without_ext})
	list(APPEND CMAKE_PREFIX_PATH ${TOOL_INSTALL_PATH})
endforeach()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED jsoncpp IMPORTED_TARGET GLOBAL)
pkg_check_modules(SPDLOG REQUIRED spdlog IMPORTED_TARGET GLOBAL)
#pkg_check_modules(CXXOPTS REQUIRED cxxopts IMPORTED_TARGET GLOBAL)

set(include_dir ${include_dirs} ${TOOLS_INCLUDE_PATH} ${rocthrust_INCLUDE_DIRS} ${SPDLOG_INCLUDE_DIRS} ${JSONCPP_INCLUDE_DIRS} "../common")  # ${jsoncpp_INCLUDE_DIRS}
set(link_dir  ${TOOLS_LINK_PATH} ${rocthrust_LIBRARIES} ${SPDLOG_LIBRARY_DIRS} ${JSONCPP_LIBRARY_DIRS})  # ${jsoncpp_LIBRARIES}
set(depends_name "roc::rocthrust" PkgConfig::JSONCPP PkgConfig::SPDLOG)

cpp_library(
name ${PROJECT_NAME}
depends_include ${include_dir}
depends_link_dir ${link_dir}
depends_name ${depends_name}
)