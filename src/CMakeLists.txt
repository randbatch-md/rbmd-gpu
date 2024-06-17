if(USE_ROCM)
	file(GLOB_RECURSE device_srcs "./op/*.cpp*")
	add_library(device OBJECT ${device_srcs})
	
	#device rocm
	file(GLOB_RECURSE device_rocm_srcs "./op/*.hip.cu*")
	hip_add_library(device_rocm STATIC device_rocm_srcs)

	#device
	target_compile_options(device PRIVATE ${HIP_HIPCC_FLAGS})
	target_link_libraries(
		device
		device_rocm
	    hip::host
		hip::device)

	message(STATUS "USE_ROCM: device link device_rocm")
endif()
