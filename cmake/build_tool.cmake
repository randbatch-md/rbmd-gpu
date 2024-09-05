define_property(GLOBAL PROPERTY TOOLS_INCLUDE_PATH BRIEF_DOCS "tools install path" FULL_DOCS "tools install path")
define_property(GLOBAL PROPERTY TOOLS_LINK_PATH BRIEF_DOCS "tools link path" FULL_DOCS "tools link path")

define_property(GLOBAL PROPERTY SPDLOG_INCLUDE_PATH BRIEF_DOCS "spdlog install path" FULL_DOCS "spdlog install path")
define_property(GLOBAL PROPERTY SPDLOG_LINK_PATH BRIEF_DOCS "spdlog link path" FULL_DOCS "spdlog link path")

macro(build_tool name)
	set(TOOL_INSTALL_PATH ${tooldir}/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/${toolname_without_ext})
	set_property(GLOBAL APPEND PROPERTY TOOLS_INCLUDE_PATH ${TOOL_INSTALL_PATH}/include)
	set_property(GLOBAL APPEND PROPERTY TOOLS_LINK_PATH ${TOOL_INSTALL_PATH}/lib64)

string(FIND ${name} "spdlog" index)
if(NOT index EQUAL -1)
	set_property(GLOBAL  PROPERTY SPDLOG_INCLUDE_PATH ${TOOL_INSTALL_PATH}/include)
	set_property(GLOBAL  PROPERTY SPDLOG_LINK_PATH ${TOOL_INSTALL_PATH}/lib64)
endif()

	include(ProcessorCount)
	ProcessorCount(N)

	if(NOT N EQUAL 0)
		set(PARALLEL_OPTION "--parallel ${N}")
	else()
		set(PARALLEL_OPTION "")
	endif()

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

		message("--- cmake --build build ${toolname_without_ext} build type: ${CMAKE_BUILD_TYPE}---")
		execute_process(COMMAND ${CMAKE_COMMAND} --build ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE} ${PARALLEL_OPTION})

		message("--- cmake --install ${toolname_without_ext} ---")
		execute_process(COMMAND ${CMAKE_COMMAND} --install ${tool_source}/build
		--config ${CMAKE_BUILD_TYPE})

		message("****** build ${toolname_without_ext} end ******")

	endif()
endmacro()