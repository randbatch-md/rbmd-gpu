define_property(GLOBAL PROPERTY TOOLS_INCLUDE_PATH)

macro(build_tool)
	set(TOOL_INSTALL_PATH ${tooldir}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/${toolname_without_ext})
	set_property(GLOBAL APPEND PROPERTY TOOLS_INCLUDE_PATH ${TOOL_INSTALL_PATH}/include)
	
	#avoid repeated build
	if(NOT EXISTS ${TOOL_INSTALL_PATH})
		message("****** build ${toolname_without_ext} start ******")

		message("--- tar -zxvf ${toolname} ---")
		execute_process(COMMAND ${CMAKE_COMMAND}
		-E tar -xzvf ${tool_path}
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

		message("--- cmake -S . -B build ${toolname_without_ext} ---")
		set(tool_source ${PROJECT_BINARY_DIR}/${toolname_without_ext})
		execute_process(COMMAND ${CMAKE_COMMAND} -S ${tool_source} -B ${tool_source}/build
		-DCMAKE_INSTALL_PREFIX=${TOOL_INSTALL_PATH} -DBUILD_SHARED_LIBS=false)

		message("--- cmake --build build ${toolname_without_ext} ---")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE})

		message("--- cmake --install ${toolname_without_ext} ---")
		execute_process(COMMAND ${CMAKE_COMMAND} --install ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE})

		message("****** build ${toolname_without_ext} end ******")

	endif()
endmacro()