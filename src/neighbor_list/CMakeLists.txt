project(neighbor_list)

get_property(cmake_dir GLOBAL PROPERTY CUSTOME_CMAKE_DIR)
include(${cmake_dir}/cpp_setting.cmake)

get_property(spdlog_include GLOBAL PROPERTY SPDLOG_INCLUDE_PATH)
get_property(spdlog_link GLOBAL PROPERTY SPDLOG_LINK_PATH)


file(GLOB_RECURSE all_include_dirs
    "${CMAKE_CURRENT_LIST_DIR}/include/*.h"
    "${CMAKE_CURRENT_LIST_DIR}/src/*.h"
)

set(include_dirs)
foreach(header_file ${all_include_dirs})
    get_filename_component(dir ${header_file} DIRECTORY)
    list(APPEND include_dirs ${dir})
endforeach()
list(REMOVE_DUPLICATES include_dirs)


set(INCLUDE_DIRS)
set(LINK_NAME)
set(LINK_DIRS)
if(USE_ROCM)
    find_package(hipcub REQUIRED)
	set(INCLUDE_DIRS ${hipcub_INCLUDE_DIRS})
	set(LINK_NAME "hip::hipcub")
    set(LINK_DIRS ${hipcub_LIBRARIES})
elseif(USE_CUDA)
	set(INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()


get_property(TOOLS_INCLUDE_PATH GLOBAL PROPERTY TOOLS_INCLUDE_PATH)
get_property(TOOLS_LINK_PATH GLOBAL PROPERTY TOOLS_LINK_PATH)
set(inclue_files ${TOOLS_INCLUDE_PATH}  "../data_manager/include"  ${include_dirs} ${spdlog_include} ${INCLUDE_DIRS})
set(link_dirs ${TOOLS_LINK_PATH}  ${spdlog_link} ${LINK_DIRS})
#set(link_name "hip::device" "spdlog" "hip::hipcub")
set(depend_names ${LINK_NAME} "data_manager")


find_package(CUDAToolkit REQUIRED)
if(CUDNN_FOUND)
    add_definitions(-DUSE_CUDNN=1)
endif()
if(CUSPARSELT_FOUND)
    add_definitions(-DUSE_CUSPARSELT=1)
endif()
# 手动定义 gflags
if(USE_MACE)

#set(TORCH_DIR "/work/home/zhaoteng_sciai/zkh/libtorch")

    set(Torch_DIR "${TORCH_DIR}/share/cmake/Torch")
    find_package(Torch REQUIRED)
#message(STATUS "Found Torch: ${Torch_FOUND}")
#message(STATUS "Torch libraries: ${TORCH_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
    list(APPEND link_dirs "${TORCH_DIR}/lib"  )
    list(APPEND depend_names ${TORCH_LIBRARIES})
    list(APPEND include_files  "${TORCH_DIR}/include" "${TORCH_DIR}/include/torch/csrc/api/include" "${TORCH_DIR}/include/torch/csrc/autograd"
      "${TORCH_DIR}/include/torch/csrc/cuda" )
endif()

cpp_library(
name ${PROJECT_NAME}
depends_include ${inclue_files}
depends_link_dir ${link_dirs}
depends_name ${depend_names}
)
