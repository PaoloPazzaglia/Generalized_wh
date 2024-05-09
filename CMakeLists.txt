cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

project (wsa)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

# Include dirs.
include_directories(.)
include_directories(data)
include_directories(${CPLEX_INCLUDE_DIR})
include_directories(/opt/ibm/ILOG/CPLEX_Studio127/cplex/include/ilcplex/)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
find_package(CPLEX)

set (SOURCE_FILES
	wsa.cpp
	wa.cpp
	wa_kill.cpp
	wa_fast.cpp
	wa_kill_fast.cpp
	util.cpp
	task.cpp)


add_executable (${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} cplex-library)
target_link_libraries(${PROJECT_NAME} cplex-concert)
target_link_libraries(${PROJECT_NAME} ilocplex)
target_link_libraries(${PROJECT_NAME} Threads::Threads)
target_link_libraries(${PROJECT_NAME} m)

MESSAGE( STATUS "CPLEX_INCLUDE_DIR: " ${CPLEX_INCLUDE_DIR} )
MESSAGE( STATUS "CPLEX_SOURCE_DIR: " ${CMAKE_SOURCE_DIR} )

#DEBUG
#target_compile_definitions(${PROJECT_NAME} PRIVATE __DEBUG_MILP__)
#target_compile_definitions(${PROJECT_NAME} PRIVATE __DEBUG_ALLOCATION__)
#target_compile_definitions(${PROJECT_NAME} PRIVATE __DEBUG_OPTBROE__)

#This is extremely verbose
#target_compile_definitions(${PROJECT_NAME} PRIVATE __DEBUG_OPTBROE_L2__)
