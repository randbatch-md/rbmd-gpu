#build tools

include(../cmake/build_tool.cmake)

file(GLOB tools_path "./*.tar.gz*")

foreach(tool_path ${tools_path})
	get_filename_component(toolname ${tool_path} NAME)
	get_filename_component(tooldir ${tool_path} PATH)
	string(REGEX REPLACE "\\.tar\\.gz$" "" toolname_without_ext ${toolname})
	message("****** build ${toolname_without_ext} start ******")
	build_tool()
	message("****** build ${toolname_without_ext} end ******")
endforeach()