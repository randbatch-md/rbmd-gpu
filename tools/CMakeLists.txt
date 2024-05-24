#build tools

include(../cmake/build_tool.cmake)

file(GLOB tools_path "./*.tar.gz*")

foreach(tool_path ${tools_path})
	get_filename_component(filename ${tool_path} NAME)
	string(REGEX REPLACE "\\.tar\\.gz$" "" filename_without_ext ${filename})
	message("****** build ${filename_without_ext} start ******")
	build_tool()
	message("****** build ${filename_without_ext} end ******")
endforeach()