project(data_manager)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)
get_property(sppdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(sppdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)

get_property(TOOLS_INCLUDE_PATH GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
get_property(TOOLS_LINK_PATH GLOBAL PROPERTY TOOLS_LINK_PATH)
set(include_dir ${TOOLS_INCLUDE_PATH}  ${sppdlog_include}  ${rocthrust_INCLUDE_DIRS} "../../tools/Windows/Debug/jsoncpp-1.9.5/include")  # ${jsoncpp_INCLUDE_DIRS}
set(link_dir  ${TOOLS_LINK_PATH}  ${sppdlog_link}  ${rocthrust_LIBRARIES}   )  # ${jsoncpp_LIBRARIES}
set(link_name "jsoncpp_lib" "spdlog" "roc::rocthrust")

cpp_library(
name ${PROJECT_NAME}
depends_include ${include_dir}
depends_link_dir ${link_dir}
depends_name ${link_name}
)