if(USE_ROCM)
	#device rocm
	set(DEVICE_ROCM_NAME "device_rocm")
	file(GLOB_RECURSE device_rocm_srcs "./op/*.hip.cu*")
	add_library(${DEVICE_ROCM_NAME} STATIC ${device_rocm_srcs})
	#target_link_libraries(${DEVICE_ROCM_NAME} PRIVATE --emit-static-lib)
	set_target_properties(${DEVICE_ROCM_NAME} PROPERTIES LINKER_LANGUAGE CXX)
	target_include_directories(${DEVICE_ROCM_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/op/near_force/direct_truncation)

	#device
	#set(DEVICE_NAME "device")
	#file(GLOB_RECURSE device_srcs "./op/*.cpp*")
	#add_library(${DEVICE_NAME} OBJECT ${device_srcs})
	#add_dependencies(${DEVICE_NAME} ${DEVICE_ROCM_NAME})
	#target_link_libraries(${DEVICE_NAME} PRIVATE -l${DEVICE_ROCM_NAME})
	#target_include_directories(${DEVICE_NAME} PUBLIC "${CMAKE_CURRENT_LIST_DIR}/op/base")
endif()
