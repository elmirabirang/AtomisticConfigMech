/*
 * modification_interaction_list_3d.cpp
 *
 *  Created on: Feb 29, 2020
 *      Author: S.Elmira Birang.O
 *      Central Institute of Scientific Computing
 *      University of Erlangen-Nuremberg
 *      This function modifies neighbor list of atoms met configurational-force-criterion for 3D lattice
 */

#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//---
#include "atom.h"
#include "force.h"
#include "point.h"

template<int dim> void ModifyNeighborList(Atoms<dim> atoms, double alpha1, double alpha2,
                                          double r_01, double r_02, double E1, double E2, double delta,
                                          double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                          double s1, double s2, double s3,  double a, double beta1,
                                          double beta2, double r_03, double r_04, double F0, double F2,
                                          double q1, double q2, double q3, double q4, double Q1, double Q2, double crack_line) {

    Force<dim> force(atoms, 0.0, 0.0);

    LOOP_OVER_ATOMS(atoms, i,
        Point<3> config_force_pair = force.configForceEamPair(i, alpha1, alpha2,
                                                              r_01, r_02, E1, E2, delta,
                                                              r_cut, h, r_s1, r_s2, r_s3,
                                                              s1, s2, s3);

        Point<3> config_force_embedding = force.configForceEamEmbedding(i, a, beta1,
                                                                        beta2, r_03, r_04, F0, F2,
                                                                        q1, q2, q3, q4, Q1, Q2,
                                                                        h, r_cut);

        Point<3> config_force = (config_force_pair + config_force_embedding) * -1;
        double length_config_force = config_force.norm();

        if(length_config_force >= 10.367) {
            Point<dim> atom_material_position = atoms->getMaterialPosition(i);
            double atom_z_coord = atom_material_position.GetZCoord();
            bool sign = (atom_z_coord - crack_line) > 0;

            LOOP_OVER_ATOM_NEIGHBORS_CUSTOM(atoms, neighbor_list, it, i, j,
                Point<dim> neighbor_material_position = atoms->getMaterialPosition(j);
                double neighbor_z_coord = neighbor_material_position.GetZCoord();
                if(((neighbor_z_coord - crack_line) < 0 && sign) || ((neighbor_z_coord - crack_line) > 0 && !sign)) {
                    it = neighbor_list.erase(it);
                } else {
                    it++;
                }
            )
        }
    )
}
