project(force)

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

# ѭ�����ӿ�
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


set(INCLUDE_DIRS)
set(LINK_NAME)
set(LINK_DIRS)
if(USE_ROCM)
	set(INCLUDE_DIRS ${rocthrust_INCLUDE_DIRS} ${HIP_INCLUDE_DIRS})
	set(LINK_NAME "hip::device" "roc::rocthrust")
    set(LINK_DIRS ${rocthrust_LIBRARIES} ${HIP_LIBRARIES})
elseif(USE_CUDA)
	set(INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

set(inclue_files "../data_manager/include" ${include_dirs} ${TOOLS_INCLUDE_PATH} ${INCLUDE_DIRS} ${SPDLOG_INCLUDE_DIRS} ${JSONCPP_INCLUDE_DIRS})
set(link_dirs ${LIBRARY_DIR} ${TOOLS_LINK_PATH} ${LINK_DIRS}  ${SPDLOG_LIBRARY_DIRS} ${JSONCPP_LIBRARY_DIRS})
set(depend_names ${LINK_NAME} PkgConfig::JSONCPP PkgConfig::SPDLOG)

cpp_library(
name ${PROJECT_NAME}
depends_include ${inclue_files}
depends_link_dir ${link_dirs}
depends_name ${depend_names}
)
set(LIB_DIR ${CMAKE_CURRENT_LIST_DIR}/../lib)
target_link_libraries(${PROJECT_NAME} neighbor_list data_manager)