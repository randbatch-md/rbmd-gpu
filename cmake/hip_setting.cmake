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

message("-----------THRUST_INCLUDE_DIRS: ${THRUST_INCLUDE_DIRS}")

include_directories(${HIP_INCLUDE_DIRS})
link_directories(${HIP_LIBRARIES})

set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})
set(CMAKE_CXX_LINKER ${HIP_HIPCC_EXECUTABLE})
