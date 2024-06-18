if(NOT DEFINED ROCM_PATH AND "${ROCM_PATH}" STREQUAL "")
	message("ROCM_PATH is not set or is empty.")
	if(NOT DEFINED ENV${ROCM_PATH} AND "$ENV{ROCM_PATH}" STREQUAL "")
		message("ENV ROCM_PATH is not set or is empty AND set the ROCM_PATH /opt/rocm.")
		set(ROCM_PATH "/opt/rocm" CACHE STRING "default rocm installation directory.")
	else()
		set(ROCM_PATH "$ENV{ROCM_PATH}" CACHE STRING "ENV rocm installation directory.")
		message("ROCM_PATH: ${ROCM_PATH}")
	endif()
else()
	message("ROCM_PATH: ${ROCM_PATH}")
endif()

#search for rocm in common locations
list(APPEND CMAKE_PREFIX_PATH ${ROCM_PATH}/hip ${ROCM_PATH})

find_package(HIP REQUIRED)
if(HIP_FOUND)
  message("Found HIP: " ${HIP_VERSION})
  message("HIP include directories: ${HIP_INCLUDE_DIRS}")
  message("HIP libraries: ${HIP_LIBRARIES}")
else()
  message(FATAL_ERROR "Could not find HIP.")
endif()

message("HIP_HIPCC_EXECUTABLE: ${HIP_HIPCC_EXECUTABLE}")
set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})
set(CMAKE_CXX_LINKER ${HIP_HIPCC_EXECUTABLE})
