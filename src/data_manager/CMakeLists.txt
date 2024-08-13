project(data_manager)

file(GLOB_RECURSE data_manager_headers "${CMAKE_CURRENT_LIST_DIR}/*.h*" )
file(GLOB_RECURSE data_manager_sources "${CMAKE_CURRENT_LIST_DIR}/*.cpp*")

add_library(${PROJECT_NAME} SHARED ${data_manager_sources} ${data_manager_headers})
target_include_directories( ${PROJECT_NAME} PUBLIC
	"${CMAKE_CURRENT_LIST_DIR}/data"
	"${CMAKE_CURRENT_LIST_DIR}/model"
	"${CMAKE_CURRENT_LIST_DIR}/memeory_scheduler"
)


