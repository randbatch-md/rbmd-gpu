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

## Set global properties for the library's include path and link path
#get_property(data_manager_include_dirs TARGET data_manager_lib PROPERTY INCLUDE_DIRECTORIES)
#set_property(GLOBAL PROPERTY DATA_MANAGER_INCLUDE_DIRS ${data_manager_include_dirs})
#set_property(GLOBAL PROPERTY DATA_MANAGER_LIB_PATH $<TARGET_FILE:data_manager_lib>)
