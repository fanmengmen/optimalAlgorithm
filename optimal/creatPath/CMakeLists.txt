
#include_directories( ${CERES_INCLUDE_DIRS})
# add_executable(${PROJECT_NAME} lineCirclefit.cpp plot.cpp)
add_executable(creatPath creatPath.cpp creatPlan.cpp)
target_link_libraries(creatPath nlohmann_json::nlohmann_json)