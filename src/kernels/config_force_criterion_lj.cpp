/*
 * config_force_criterion_lj.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <iostream>
#include <vector>
//---
#include "force.h"
#include "atom.h"

template<int dim> Atoms<dim> *ConfigForceCriterion(Atoms<dim> *atoms) {
	Force<dim> force(1.0, 1.0);

    LOOP_OVER_ATOMS(atoms, i,
        std::vector<int> updated_neighbor_list;
        Point<2> atomMaterialPosition = atoms->getMaterialPosition(i);
        double atomMaterialPositionX = atomMaterialPosition.GetXCoord();
        double atomMaterialPositionY = atomMaterialPosition.GetYCoord();
        double crackRegion_MinX = 14.6;
        double crackRegion_MinY = 22.86;
        double crackRegion_MaxX = 58.;
        double crackRegion_MaxY = 23.83;

        if(atomMaterialPositionX >= crackRegion_MinX && atomMaterialPositionY >= crackRegion_MinY &&
           atomMaterialPositionX <= crackRegion_MaxX && atomMaterialPositionY <= crackRegion_MaxY) {

            LOOP_OVER_ATOM_NEIGHBORS(atoms, i, j,
                Point<2> neighborMaterialPosition = atoms->getMaterialPosition(j);
                double neighborMaterialPositionX = neighborMaterialPosition.GetXCoord();
                double neighborMaterialPositionY = neighborMaterialPosition.GetYCoord();
                auto material_delta = atoms->getMaterialPosition(i) - atoms->getMaterialPosition(j);
                auto distance = material_delta.norm();
                auto initial_delta = atoms.getInitialPosition(i) - atoms->getInitialPosition(j);
                auto initial_distance = initial_delta.norm();
                double criticalConfigForce = (3.77 * 1.0 / (distance * initial_distance)) * 0.999;
                double configForce = force.interatomicConfigForceValue(i, j);

                if(abs(configForce) < criticalConfigForce) {
                    updated_neighbor_list.push_back(j);
                } else if(abs(configForce) >= criticalConfigForce &&
                          neighborMaterialPositionX >= crackRegion_MinX && neighborMaterialPositionY >= crackRegion_MinY &&
                          neighborMaterialPositionX <= crackRegion_MaxX && neighborMaterialPositionY <= crackRegion_MaxY) {

                    std::cout << "criticalConfigForce: " << criticalConfigForce <<
                                 " configForce: " << abs(configForce) <<
                                 " neighbor id : " << j <<
                                 " atom id: " << i << std::endl;
                }
	    	)

            atoms->setNeighbors(i, updated_neighbor_list);

            /*
            std::cout << "update_neighbor_list size: " << update_neighbor_list.size() <<
                         " neighbors size: " << neighbors.size() << std::endl;
            */
		}
	)

	return atoms;
}

