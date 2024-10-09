if(USE_CUDA)
	find_package(CUDA REQUIRED)
	message("CUDA verson: ${CUDA_VERSION}")
	message("CUDA include directories: ${CUDA_INCLUDE_DIRS}")
	message("CUDA libraries: ${CUDA_LIBRARIES}")
endif()


if(USE_ROCM)
	if(NOT DEFINED HIP_PATH AND "${HIP_PATH}" STREQUAL "")
		if(NOT DEFINED ENV${HIP_PATH} AND "$ENV{HIP_PATH}" STREQUAL "")
			message("HIP_PATH is not set or is empty AND set the HIP_PATH /opt/rocm/hip.")
			set(HIP_PATH "/opt/rocm/hip" CACHE STRING "default hip installation directory.")
		else()
			set(HIP_PATH "$ENV{HIP_PATH}" CACHE STRING "ENV hip installation directory.")
			message("HIP_PATH: ${HIP_PATH}")
		endif()
	else()
		message("HIP_PATH: ${HIP_PATH}")
	endif()

	#search for rocm in common locations
	#list(APPEND CMAKE_PREFIX_PATH ${HIP_PATH})
	#set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
	set(HIP_HIPCC_FLAGS -fno-gpu-rdc; --std=c++17)

	find_package(HIP REQUIRED)
	if(HIP_FOUND)
	  message("Found HIP: " ${HIP_VERSION})
	  message("HIP include directories: ${HIP_INCLUDE_DIRS}")
	  message("HIP libraries: ${HIP_LIBRARIES}")
	else()
	  message(FATAL_ERROR "Could not find HIP.")
	endif()

	include_directories(${HIP_INCLUDE_DIRS})
	link_directories(${HIP_LIBRARIES})

	message("HIP_HIPCC_EXECUTABLE:${HIP_HIPCC_EXECUTABLE}")
	#set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})
	#set(CMAKE_CXX_LINKER ${HIP_HIPCC_EXECUTABLE})

	#rocthrust
	find_package(rocthrust REQUIRED)
	if(rocthrust_FOUND)
		message("rochthrust: ${rocthrust_DIR}")
		message("rochthrust include directories: ${rocthrust_INCLUDE_DIRS}")
		message("rochthrust libraries: ${rocthrust_LIBRARIES}")
	endif()
endif()
