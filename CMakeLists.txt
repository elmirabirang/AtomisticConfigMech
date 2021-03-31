cmake_minimum_required(VERSION 2.8)
project(AtomisticConfigMech)
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# 3D Energy Release (EAM)
add_executable(engrelease3d test_eam_cu.cpp)
target_link_libraries(engrelease3d ${CERES_LIBRARIES})

# 2D Energy Release (LJ)
add_executable(eng_release2d energy_release_2d.cpp)
target_link_libraries(eng_release2d ${CERES_LIBRARIES})

# Silicon (Stillinger-Weber)
add_executable(silicon test_silicon.cpp)
target_link_libraries(silicon ${CERES_LIBRARIES})
