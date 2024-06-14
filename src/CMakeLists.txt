find_package(HIP REQUIRED)

file(GLOB_RECURSE device_srcs "./op/*.cpp*")
add_library(device OBJECT ${device_srcs})

if(USE_ROCM)
file(GLOB_RECURSE device_rocm_srcs "./op/*.hip.cu*")
hip_add_library(device_rocm STATIC device_rocm_srcs)
target_link_libraries(device device_rocm)
message(STATUS "USE_ROCM: device link device_rocm")
endif()
