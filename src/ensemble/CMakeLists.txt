project(ensemble)

file(GLOB_RECURSE ensemble_headers "${CMAKE_CURRENT_LIST_DIR}/*.h*" )
file(GLOB_RECURSE ensemble_sources "${CMAKE_CURRENT_LIST_DIR}/*.cpp*")

add_library(${PROJECT_NAME} SHARED ${ensemble_headers} ${ensemble_sources})
target_include_directories( ${PROJECT_NAME} PUBLIC
	"${CMAKE_CURRENT_LIST_DIR}/ensembles"
)


