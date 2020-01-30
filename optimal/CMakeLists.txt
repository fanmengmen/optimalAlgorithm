
#CMakeLists模板

cmake_minimum_required( VERSION 2.8 )
project(ShapeEightFitting)
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
add_subdirectory(./creatPath)
#include_directories( ${CERES_INCLUDE_DIRS})
add_executable(${PROJECT_NAME} lineCirclefit.cpp plot.cpp)
# add_executable(${PROJECT_NAME} creatPath.cpp creatPlan.cpp)
target_link_libraries(${PROJECT_NAME} ceres ${OpenCV_LIBS})
# target_link_libraries(${PROJECT_NAME} nlohmann_json::nlohmann_json)