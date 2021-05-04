/*
 * neighbors.cpp
 * get the list of atoms from unrelax_config.cpp
 * construct the neighbors of atom
 *
 *  Created on: Dec 3, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical enigneering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "bond.h"
#include "point.h"
#include "cell.h"
#include <iostream>
#include <math.h>
#include <vector>
//---
#include "subcell.h"
#include "atom.h"

// TODO: use cell lists to build neighbor lists
template<int dim> void FindNeighbors(Atoms<dim> *atoms, double cutoff_radius) {
    LOOP_OVER_ATOMS(atoms, i,
		Point<dim> pos_i = atoms->getMaterialPosition(i);
		int region_i = atoms->getRegion(i);

        if(region_i != 3) {
            LOOP_OVER_ATOMS(atoms, j,
                int region_j = atoms->getRegion(j);

                if(i < j && region_j != 3) {
                    if((region_i == 4 && region_j == 5) || (region_i == 5 && region_j == 4)) {
                        continue;
                    }

                    Point<dim> pos_j = atoms->getMaterialPosition(j);
                    Point<dim> delta = pos_i - pos_j;
                    if(delta.norm() < cutoff_radius) {
                        atoms->appendNeighbor(i, j);
                        atoms->appendBondNeighbor(j, i);
                        //neighbors.push_back(*atomj);
                    }
                }
			)
		}
	)
}

template void FindNeighbors<2>(Atoms<2> *atoms, double cutoff_radius);
template void FindNeighbors<3>(Atoms<3> *atoms, double cutoff_radius);
