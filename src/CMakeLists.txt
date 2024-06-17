if(USE_ROCM)
	if(NOT DEFINED ROCM_PATH)
		set(ROCM_PATH "/opt/rocm" CACHE STRING "Default ROCM installation directory.")
	endif()
	message(STATUS "ROCM_PATH: ${ROCM_PATH}")

	if(NOT DEFINED HIP_PATH)
		if(NOT DEFINED ENV{HIP_PATH})
			set(HIP_PATH "${ROCM_PATH}/hip" CACHE PATH "hip path")
		else()
			set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "env hip path")
		endif()
	endif()
	message(STATUS "HIP_PATH: ${HIP_PATH}")
	set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
	set(HIP_HIPCC_FLAGS -fno-gpu-rdc; --std=c++17)

	find_package(HIP REQUIRED)
	if(HIP_FOUND)
		message(STATUS "HIP found and the version: ${HIP_VERSION}")
		message(STATUS "HIP include directories: ${HIP_INCLUDE_DIRS}")
		message(STATUS "HIP libraries: ${HIP_LIBRARIES}")
	else()
	    message(FATAL_ERROR "Hip not found.")
	endif()

	file(GLOB_RECURSE device_srcs "./op/*.cpp*")
	add_library(device OBJECT ${device_srcs})
	
	#device rocm
	file(GLOB_RECURSE device_rocm_srcs "./op/*.hip.cu*")
	hip_add_library(device_rocm STATIC device_rocm_srcs)

	#device
	target_compile_options(device PRIVATE ${HIP_HIPCC_FLAGS})
	target_link_libraries(device PRIVATE device_rocm)
	message(STATUS "USE_ROCM: device link device_rocm")
endif()
