project(controller)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)

get_property(data_manager_include_dirs GLOBAL PROPERTY data_manager_lib_INCLUDE_DIRS)
get_property(data_manager_lib_path GLOBAL PROPERTY data_manager_lib_LIB_PATH)

get_property(spdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(spdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)

#cpp_library(
#name ${PROJECT_NAME}
#depends_include ${spdlog_include}
#depends_link_dir ${spdlog_link}
#depends_name "spdlog"
#)
cpp_library(
name ${PROJECT_NAME}
depends_include "${spdlog_include};${data_manager_include_dirs}"
depends_link_dir "${spdlog_link};${data_manager_lib_path}"
depends_name "spdlog;data_manager"
)