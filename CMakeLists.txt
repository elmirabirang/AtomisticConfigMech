#cmake_minimum_required(VERSION 2.8)

#project(engrelease3d)

#find_package(Ceres REQUIRED)
#include_directories(${CERES_INCLUDE_DIRS})

# helloworld
#add_executable(engrelease3d test_eam_cu.cpp)
#target_link_libraries(engrelease3d ${CERES_LIBRARIES})

#cmake_minimum_required(VERSION 2.8)

project(engrelease2d)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_executable(eng_release2d energy_release_2d.cpp)
target_link_libraries(eng_release2d ${CERES_LIBRARIES})

#cmake_minimum_required(VERSION 2.8)

#project(silicon)

#find_package(Ceres REQUIRED)
#include_directories(${CERES_INCLUDE_DIRS})

# helloworld
#add_executable(silicon test_silicon.cpp)
#target_link_libraries(silicon ${CERES_LIBRARIES})
