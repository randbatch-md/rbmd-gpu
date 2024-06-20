#if(NOT DEFINED ROCM_PATH AND "${ROCM_PATH}" STREQUAL "")
#	message("ROCM_PATH is not set or is empty.")
#	if(NOT DEFINED ENV${ROCM_PATH} AND "$ENV{ROCM_PATH}" STREQUAL "")
#		message("ENV ROCM_PATH is not set or is empty AND set the ROCM_PATH /opt/rocm.")
#		set(ROCM_PATH "/opt/rocm" CACHE STRING "default rocm installation directory.")
#	else()
#		set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE STRING "ENV rocm installation directory.")
#		message("ROCM_PATH: ${ROCM_PATH}")
#	endif()
#else()
#	message("ROCM_PATH: ${ROCM_PATH}")
#endif()

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
#list(APPEND CMAKE_PREFIX_PATH ${ROCM_PATH}/hip ${ROCM_PATH})
set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
set(HIP_CLANG_FLAGS -fno-gpu-rdc; --std=c++17)
#set(HIP_HIPCC_FLAGS "${HIP_HIPCC_FLAGS}" "-std=c++11 -fgpu-rdc" )

find_package(HIP REQUIRED)
if(HIP_FOUND)
  message("Found HIP: " ${HIP_VERSION})
else()
  message(FATAL_ERROR "Could not find HIP.")
endif()

#set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})
#set(CMAKE_CXX_LINKER ${HIP_HIPCC_EXECUTABLE})
