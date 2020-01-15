
#CMakeLists模板

cmake_minimum_required( VERSION 2.8 )
project(ShapeEightFitting)
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)
#include_directories( ${CERES_INCLUDE_DIRS})
add_executable(${PROJECT_NAME} lineCirclefit.cpp plot.cpp)
target_link_libraries(${PROJECT_NAME} ceres ${OpenCV_LIBS})
