project(data_manager)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)

get_property(sppdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(sppdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)

set(include_dir ${sppdlog_include} ${rocthrust_INCLUDE_DIRS})
set(link_dir ${sppdlog_link} ${rocthrust_LIBRARIES})
set(link_name "spdlog" "roc::rocthrust")

cpp_library(
name ${PROJECT_NAME}
depends_include ${include_dir}
depends_link_dir ${link_dir}
depends_name ${link_name}
)