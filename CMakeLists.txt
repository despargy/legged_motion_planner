cmake_minimum_required(VERSION 3.0.2)
project(legged_motion_planner)
# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
set ( libipopt.a )


## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  rospy
  std_msgs
)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)

find_package(PkgConfig)
pkg_check_modules(EIGEN3 eigen3 REQUIRED)


catkin_package(
  INCLUDE_DIRS include
  LIBRARIES legged_motion_planner
  CATKIN_DEPENDS roscpp roslib rospy std_msgs
#  DEPENDS system_lib
)


include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(
 include
  ${catkin_INCLUDE_DIRS}
)
include_directories("${PROJECT_BINARY_DIR}" /usr/local/include/coin-or/)


# creation of an executable
add_executable(planner src/main.cpp src/CentroidalNLP.cpp src/Agent.cpp  src/NumericalIntegrationConstraints.cpp)
target_link_libraries(planner ${catkin_LIBRARIES} ipopt)
#set_target_properties(planner PROPERTIES COMPILE_FLAGS "-std=c++11 -o3")
