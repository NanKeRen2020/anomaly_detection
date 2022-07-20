
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Utilities)
include_directories(${CMAKE_SOURCE_DIR}/model)


#message("-- Message LAPACK libs: " ${LAPACK_LIBRARIES})
#message("-- LAPACK incs: " ${LAPACK_INCLUDES})
include_directories(${LAPACK_INCLUDES})
include_directories(${CBLAS_INCLUDE_DIR})
include_directories(${FFTW_INCLUDEs})

find_package(Boost REQUIRED COMPONENTS system filesystem thread)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})


link_directories(${CMAKE_BINARY_DIR}/lib)

add_executable(main_anomaly_detection main_anomaly_detection.cpp)


target_link_libraries(main_anomaly_detection modelLibrary png ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES} ${FFTW_LIBRARIES} ${OpenCV_LIBS} -lstdc++fs)




