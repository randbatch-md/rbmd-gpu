if(USE_ROCM)
	file(GLOB_RECURSE device_srcs "./op/*.cpp*")
	add_library(device OBJECT ${device_srcs})
	
	#device rocm
	file(GLOB_RECURSE device_rocm_srcs "./op/*.hip.cu*")
	hip_add_library(device_rocm STATIC device_rocm_srcs)

	#set_target_properties(device_rocm PROPERTIES
	#	INTERFACE_COMPILE_OPTIONS "-fPIC"
	#	INTERFACE_LINK_LIBRARIES "hip::host")
	#target_link_libraries (device_rocm hip::host hip::device)

	#device
	target_compile_options(device_rocm PRIVATE ${HIP_HIPCC_FLAGS})
	target_compile_options(device_rocm PRIVATE _ROCM)

	target_link_libraries(
		device
		device_rocm)

	message(STATUS "USE_ROCM: device link device_rocm")
endif()
