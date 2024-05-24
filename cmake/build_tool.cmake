macro(build_tool)
	set(TOOL_INSTALL_PATH ${tool_path}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/filename_without_ext)	
	#avoid repeated build
	#if(NOT EXISTS ${TOOL_INSTALL_PATH})
		message("--- tar -zxvf ${filename} ---")
		execute_process(COMMAND ${CMAKE_COMMAND}
		-E tar -xzvf ${tool_path}
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

		message("--- cmake -S . -B build ${filename} ---")
		set(tool_source ${PROJECT_BINARY_DIR}/${filename_without_ext})
		execute_process(COMMAND ${CMAKE_COMMAND} -S ${tool_source} -B ${tool_source}/build
		-DCMAKE_INSTALL_PREFIX=${TOOL_INSTALL_PATH} -DBUILD_SHARED_LIBS=false)

		message("--- cmake --build build ${filename} ---")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE})

		message("--- cmake --install ${filename} ---")
		execute_process(COMMAND ${CMAKE_COMMAND} --install ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE})
	#endif()
endmacro()