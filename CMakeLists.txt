cmake_minimum_required(VERSION 3.1)
project(onko3D-3.0)

set (CMAKE_C_COMPILER gcc-8)
set (CMAKE_CXX_COMPILER g++-8)
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++14 -Wall -Wextra -pedantic")#



# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# CGAL
find_package(CGAL REQUIRED COMPONENTS Core io)
include( ${CGAL_USE_FILE} ) 
MESSAGE( STATUS "EXECUTABLE_OUTPUT_PATH: " ${CGAL_USE_FILE} )

# openMP & multithread
find_package(Threads REQUIRED)
find_package( OpenMP)
if(OPENMP_FOUND)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
   include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
endif()

# mecacell 
set(MECACELL_DIR src/external/MecaCell)
include_directories(${MECACELL_DIR})

# source
file(GLOB SRC
    "src/*.h"
    "src/*.hpp"
    "src/core/*.h"
    "src/core/*.hpp"
    "src/core/*.cpp"
	"src/plugin/delaunaymassspringdamper/*.h"
	"src/plugin/delaunaymassspringdamper/*.hpp"
	"src/plugin/delaunaymassspringdamper/*.cpp"
    )

# ---------- EXECUTABLES -----------
#console
add_executable(console ${SRC} src/mainconsole.cpp)
target_link_libraries(console boost_program_options)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS console)
target_link_libraries(console ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
if(THREADS_HAVE_PTHREAD_ARG)
target_compile_options(PUBLIC console "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
target_link_libraries(console "${CMAKE_THREAD_LIBS_INIT}")
endif()

# viewer
set(VIEWER_NAME viewer) # nom de l'executable
find_package(Qt5 COMPONENTS Core Quick REQUIRED) # find Qt
set(MECACELL_VIEWER_DIR ${MECACELL_DIR}/mecacellviewer)
set(CMAKE_AUTOMOC ON)
qt5_wrap_cpp(MOC_SRC ${MECACELL_VIEWER_DIR}/signalslotbase.h ${MECACELL_VIEWER_DIR}/camera.hpp)
qt5_add_resources(RESOURCES ${MECACELL_VIEWER_DIR}/resourcesLibMecacellViewer.qrc)
add_executable(${VIEWER_NAME} ${SRC} ${MOC_SRC} ${RESOURCES} src/mainviewer.cpp)
qt5_use_modules(${VIEWER_NAME} Core Quick)
target_link_libraries(${VIEWER_NAME} boost_program_options)
add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${VIEWER_NAME} )
target_link_libraries(${VIEWER_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

if(THREADS_HAVE_PTHREAD_ARG)
target_compile_options(PUBLIC ${VIEWER_NAME} "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
target_link_libraries(${VIEWER_NAME} "${CMAKE_THREAD_LIBS_INIT}")
endif()

