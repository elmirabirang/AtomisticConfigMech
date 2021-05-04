/*
 * broken_bonds.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <iostream>
#include <vector>
//---
#include "force.h"
#include "atom.h"

template<int dim> bool brokenBonds(Atoms<dim> *atoms) {
    Force<dim> force(1.0, 1.0);
    int count = 0;
    bool bond_broken = false;

    LOOP_OVER_ATOMS(atoms, i,
        Point<dim> atom_position = atoms->getMaterialPosition(i);
        double atomp_x = atom_position.GetXCoord();
        double atomp_y = atom_position.GetYCoord();
        int atom_region = atoms->getRegion(i);

        if(atomp_x < 40.0 && atomp_x > 11.4 && atom_region != 1 && atom_region != 2 && atomp_y > 13.2 && atomp_y < 14.3) {
            std::vector<int> update_neighbor_list;

            LOOP_OVER_ATOM_NEIGHBORS(atoms, i, j,
                Point<dim> config_interaction_force = force.MaterialBondForce(i, j);
                // scale the pairwise configurational force
                double config_force_norm = config_interaction_force.norm();

                // scale critical configurational force
                Point<dim> critical_config_force = force.CriticalMaterialBondForce(i, j, 1.5);
                double critical_config_force_norm = critical_config_force.norm() - 0.000034;
                std::cout << critical_config_force_norm << "   " << config_force_norm << std::endl;

                if(config_force_norm < critical_config_force_norm) {
                    //cout << critical_config_force_norm << "   " << config_force_norm << endl;
                    update_neighbor_list.push_back(j);
                    // FIXME: What is really the logic here? If at least one bond is broken then return true?
                    // If that is the case, next line should be removed!
                    bond_broken = false;
                } else if(config_force_norm > critical_config_force_norm) {
                    bond_broken = true;
                }
            )

            atoms->setNeighbors(update_neighbor_list);
        }
    )

    return bond_broken;
}
