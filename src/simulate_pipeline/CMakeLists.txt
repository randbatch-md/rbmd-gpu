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


# 循环链接库
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

set(include_dir ${include_dirs} ${TOOLS_INCLUDE_PATH} ${rocthrust_INCLUDE_DIRS} ${SPDLOG_INCLUDE_DIRS} ${JSONCPP_INCLUDE_DIRS} ${hipcub_INCLUDE_DIRS} ${HIP_INCLUDE_DIRS} "../common" "../data_manager/include")  # ${jsoncpp_INCLUDE_DIRS}
set(link_dir  ${TOOLS_LINK_PATH} ${rocthrust_LIBRARIES} ${hipcub_LIBRARIES} ${HIP_LIBRARIES} ${SPDLOG_LIBRARY_DIRS} ${JSONCPP_LIBRARY_DIRS})  # ${jsoncpp_LIBRARIES}
set(depends_name "hip::device" "roc::rocthrust" PkgConfig::JSONCPP PkgConfig::SPDLOG)



#set(include_dir ${include_dirs} ${spdlog_include} ${rocthrust_INCLUDE_DIRS} ${hipcub_INCLUDE_DIRS} ${HIP_INCLUDE_DIRS} "../data_manager/include" "../common/")
#set(link_dirs ${LIBRARY_DIR} ${spdlog_link} ${rocthrust_LIBRARIES} ${hipcub_LIBRARIES} ${HIP_LIBRARIES})
#set(depend_names "data_manager" "controller" "roc::rocthrust")
#set(link_name "hip::hipcub" "roc::rocthrust" PkgConfig::JSONCPP PkgConfig::SPDLOG)


cpp_library(
name ${PROJECT_NAME}
depends_include ${include_dir}
depends_link_dir ${link_dirs}
depends_name ${depend_names}
)