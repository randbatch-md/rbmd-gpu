project(data_manager)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)

get_property(sppdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(sppdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)
cpp_library(
name ${PROJECT_NAME}
depends_include ${sppdlog_include}
depends_link_dir ${sppdlog_link}
depends_name "spdlog"
)

