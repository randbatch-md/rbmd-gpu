set(RUNTIME_DIR ${CMAKE_CURRENT_LIST_DIR}/../bin)
set(LIBRARY_DIR ${CMAKE_CURRENT_LIST_DIR}/../lib)

macro(config_cpp name)
	target_include_directories(${name} PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/
	${CMAKE_CURRENT_LIST_DIR}/../
	${CMAKE_CURRENT_LIST_DIR}/include/)

	target_compile_features(${name} PRIVATE cxx_std_17)
	#target_compile_features(${name} PRIVATE cxx_std_14)

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
	#Recursive search method
    # Recursively find all .cpp and .cxx files
    file(GLOB_RECURSE SRC ${CMAKE_CURRENT_LIST_DIR}/*.cpp ${CMAKE_CURRENT_LIST_DIR}/*.cxx)
    # message("${name} SRC: " ${SRC})
	
    # Find interface headers in the include directories
    file(GLOB_RECURSE H_FILE_I ${CMAKE_CURRENT_LIST_DIR}/include/*.h)
    # message("${name} H_FILE_I: " ${H_FILE_I})

	# Recursively find all src .h files
    file(GLOB_RECURSE H_FILE ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
    # message("${name} H_FILE: " ${H_FILE})
	
	# Recursively find all .cu files
    file(GLOB_RECURSE CU_FILE ${CMAKE_CURRENT_LIST_DIR}/src/*.cu)
    # message("${name} CU_FILE: " ${CU_FILE})
    if(NOT USE_MACE)
        set(MACE_FILES
         "${CMAKE_CURRENT_LIST_DIR}/src/mace_memory_scheduler.cpp"
         "${CMAKE_CURRENT_LIST_DIR}/src/op/rocm/mace_neighbor_list_op.hip.cu"
         "${CMAKE_CURRENT_LIST_DIR}/include/maceload.h"
         "${CMAKE_CURRENT_LIST_DIR}/src/maceload.cpp"
         "${CMAKE_CURRENT_LIST_DIR}/src/op/mace_neighbor_list_op.h"
         "${CMAKE_CURRENT_LIST_DIR}/include/neighbor_list_builder/mace_neighbor_list_builder.h"
         "${CMAKE_CURRENT_LIST_DIR}/src/neighbor_list_builder/mace_neighbor_list_builder.cpp"
         "${CMAKE_CURRENT_LIST_DIR}/include/scheduler/mace_memory_scheduler.h"
         "${CMAKE_CURRENT_LIST_DIR}/include/force_field/mace_force_field_data.h"
        )

        foreach(file ${MACE_FILES})
            list(REMOVE_ITEM SRC ${file})
            list(REMOVE_ITEM CU_FILE ${file})
            list(REMOVE_ITEM H_FILE ${file})
            list(REMOVE_ITEM H_FILE_I ${file})
        endforeach()
    endif()

endmacro()

#default static library
function(cpp_library)
	cmake_parse_arguments(
		"lib" 
		""
		"name"
		"depends_include;depends_link_dir;depends_name"
		${ARGN}
	)

	message(STATUS "====================${lib_name} library begin====================")
	#get src file and header file
	get_src_include()
	if(USE_CUDA)
		set_source_files_properties(${SRC} ${H_FILE} ${H_FILE_I} ${CU_FILE} PROPERTIES LANGUAGE CUDA)
	elseif(USE_ROCM)
		set_source_files_properties(${CU_FILE} PROPERTIES LANGUAGE CXX)
	endif()

	#add static library
	add_library(${lib_name} STATIC ${SRC} ${H_FILE} ${H_FILE_I} ${CU_FILE})

	if(USE_CUDA)
		target_compile_options(${lib_name} PRIVATE --extended-lambda -w)
	endif()

	config_cpp(${lib_name})

	#depends
	target_include_directories(${lib_name} PRIVATE ${lib_depends_include})
	target_link_directories(${lib_name} PRIVATE ${lib_depends_link_dir})
        if(USE_MACE)  
	target_link_libraries(${lib_name} PRIVATE ${lib_depends_name})
        else()
        target_link_libraries(${lib_name} ${lib_depends_name})
        endif()
	#install
	#to do

	message(STATUS "====================${lib_name} library end======================")
endfunction()
