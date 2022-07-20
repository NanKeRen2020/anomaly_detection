cmake_minimum_required(VERSION 3.8)

project(AnomalyDetection)


message(STATUS "Configuring main project: " ${CMAKE_PROJECT_NAME})

set(PROJECT_VERSION_MAJOR  1)
set(PROJECT_VERSION_MINOR  0)
set(PROJECT_VERSION_PATCH  0)

set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
    CACHE INTERNAL "The project version.")

###############################################################################

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_BUILD_TYPE Debug)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lfftw3 -lfftw3f -g -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lfftw3 -lfftw3f  -O3 -w -ffast-math -fomit-frame-pointer -fstrength-reduce")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lfftw3 -lfftw3f  -O3 -w -ffast-math -fomit-frame-pointer -fstrength-reduce")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")


###############################################################################

cmake_minimum_required(VERSION 3.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
#set(CMAKE_VERBOSE_MAKEFILE ON)

###############################################################################

#find_package( CBLAS REQUIRED )
#find_package( LAPACKE REQUIRED )
find_package (TIFF REQUIRED)
find_package (JPEG REQUIRED)
find_package (PNG REQUIRED)
#find_package (Boost REQUIRED)
include_directories (PUBLIC ${TIFF_INCLUDE_DIR} PUBLIC ${JPEG_INCLUDE_DIR} PUBLIC ${PNG_INCLUDE_DIRS} PUBLIC ${FFTW_INCLUDES})

###############################################################################

find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3f)
link_libraries (${TIFF_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${FFTWF_LIBRARIES})

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
#	message("-- OpenMP found.")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

find_package(OpenCV REQUIRED)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11")

# If the package has been found, several variables will
# be set, you can find the full list with descriptions
# in the OpenCVConfig.cmake file.
# Print some message showing some of them
message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

if(CMAKE_VERSION VERSION_LESS "2.8.11")
  # Add OpenCV headers location to your include paths
  include_directories(${OpenCV_INCLUDE_DIRS})
endif()

find_package(Boost REQUIRED COMPONENTS system filesystem thread)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})


add_subdirectory(src/model)
add_subdirectory(src)
