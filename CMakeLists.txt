cmake_minimum_required(VERSION 2.8)
project(AtomisticConfigMech)
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

set(CORE_CPP_FILENAMES
    assign_region.cpp
    atom.cpp
    atom_stress_2d.cpp
    atoms_in_cell.cpp
    bond.cpp
    boundary.cpp
    cell.cpp
    cell_neighbors.cpp
    crack_atoms.cpp
    #dynamical_matrix.cpp
    first_order_tensor.cpp
    generate_cells.cpp
    generate_waves.cpp
    matrix.cpp
    modification_interaction_list_3d.cpp
    neighbors.cpp
    neighbors_cell_based.cpp
    point.cpp
    read_dump_file.cpp
    subcell.cpp
    third_order_tensor.cpp
    unrelaxed_config.cpp
    write_data_file.cpp
)

PREPEND(CORE_CPP_FILES "src/core" ${CORE_CPP_FILENAMES})

set(KERNELS_CPP_FILENAMES
    LennardJones.cpp
    Stillinger_Weber.cpp
    broken_bonds.cpp
    config_force_criterion.cpp
    config_force_criterion_lj.cpp
    energy.cpp
    force.cpp
    modified_tersoff.cpp
)

PREPEND(KERNELS_CPP_FILES "src/kernels" ${KERNELS_CPP_FILENAMES})

set(HEADER_FILENAMES
    atom.h
    bond.h
    boundary.h
    cell.h
    ceres.h
    crack_atoms.h
    dynamical_matrix.h
    energy.h
    first_order_tensor.h
    force.h
    lbfgs.h
    Lennard_Jones.h
    matrix.h
    minimizer.h
    modified_tersoff.h
    point.h
    Stillinger_Weber.h
    subcell.h
    third_order_tensor.h
    unrelaxed_config.h
)

PREPEND(HEADER_FILES "include" ${HEADER_FILENAMES})

# 3D Energy Release (EAM)
add_executable(engrelease3d src/test/test_eam_cu.cpp ${CORE_CPP_FILES} ${KERNELS_CPP_FILES} ${HEADER_FILES})
target_include_directories(engrelease3d PUBLIC "include")
target_link_libraries(engrelease3d ${CERES_LIBRARIES})

# 2D Energy Release (LJ)
add_executable(eng_release2d src/test/energy_release_2d.cpp ${CORE_CPP_FILES} ${KERNELS_CPP_FILES} ${HEADER_FILES})
target_include_directories(eng_release2d PUBLIC "include")
target_link_libraries(eng_release2d ${CERES_LIBRARIES})

# Silicon (Stillinger-Weber)
add_executable(silicon src/test/test_silicon.cpp ${CORE_CPP_FILES} ${KERNELS_CPP_FILES} ${HEADER_FILES})
target_include_directories(silicon PUBLIC "include")
target_link_libraries(silicon ${CERES_LIBRARIES})
