set(RUNTIME_DIR ${CMAKE_CURRENT_LIST_DIR}/../bin)
set(LIBRARY_DIR ${CMAKE_CURRENT_LIST_DIR}/../lib)

macro(config_cpp name)
	target_include_directories(${name} PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/
	${CMAKE_CURRENT_LIST_DIR}/../
	${CMAKE_CURRENT_LIST_DIR}/include/)

	target_compile_features(${name} PRIVATE cxx_std_17)

	#default release
	if(CMAKE_BUILD_TYPE STREQUAL "")
		set(CMAKE_BUILD_TYPE Release)
	endif()

	set(CONFIG_TYPES Debug Release RelWithDebInfo MinSizeRel)
	list(APPEND CONFIG_TYPES "")
	foreach(type IN LISTS CONFIG_TYPES)
		set(conf "")
		if(type)
			string(TOUPPER _${type} conf)
			#message("conf: ${conf}")
		endif()
		set_target_properties(${name} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY${conf} ${RUNTIME_DIR}
		LIBRARY_OUTPUT_DIRECTORY${conf} ${LIBRARY_DIR}
		ARCHIVE_OUTPUT_DIRECTORY${conf} ${LIBRARY_DIR})
	endforeach()
endmacro()

#get src file and header file
macro(get_src_include)
	#find .cpp/.cxx
	aux_source_directory(${CMAKE_CURRENT_LIST_DIR} SRC)
	message("${name} SRC: " ${SRC})

	#find .h
	file(GLOB H_FILE ${CMAKE_CURRENT_LIST_DIR}/*.h)
	message("${name} H_FILE: " ${H_FILE})

	#find interface
	file(GLOB H_FILE_I ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
	message("${name} H_FILE_I: " ${H_FILE_I})
endmacro()

#default static library
function(cpp_library)
	message(STATUS "====================${name} library begin====================")
	cmake_parse_arguments(
		"lib" 
		""
		"name"
		"depends_include;depends_link_dir;depends_name"
		${ARGN}
	)

	#get src file and header file
	get_src_include()

	#add static library
	add_library(${lib_name} STATIC ${SRC} ${H_FILE} ${H_FILE_I})

	config_cpp(${lib_name})

	#depends
	target_include_directories(${lib_name} PRIVATE ${lib_depends_include})
	target_link_directories(${lib_name} PRIVATE ${lib_depends_link_dir})
	target_link_libraries(${lib_name} ${lib_depends_name})

	#install
	#to do

	message(STATUS "====================${name} library end======================")
endfunction()