/*
 * atoms_in_cell.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */

#include "atom.h"
#include "cell.h"
#include "math.h"
#include "point.h"

template<int dim> void AtominCell(Atoms<dim> *atoms, vector<Cell<dim> *> &cells, double r_cut) {
	for(auto cell = cells.begin(); cell != cells.end(); ++cell) {
		vector<int> cell_atoms;
        LOOP_OVER_ATOMS(atoms, i,
			int cell_id = atoms->getCell(i);
            Point<dim> material_position = atoms->getMaterialPosition(i);
            double atom_x = material_position.GetXCoord();
            double atom_y = material_position.GetYCoord();
            double atom_z = material_position.GetZCoord();
            Point<dim> cell_origin = (*cell)->GetOrigin();
            double cell_x_origin = cell_origin.GetXCoord();
            double cell_y_origin = cell_origin.GetYCoord();
            double cell_z_origin = cell_origin.GetZCoord();
            Point<dim> cell_tr = (*cell)->GetTopRight();
            double cell_x_tr = cell_tr.GetXCoord();
            double cell_y_tr = cell_tr.GetYCoord();
            double cell_z_tr = cell_tr.GetZCoord();

			if(cell_id == -1) {
                switch (dim) {
                    case 2:
                        if(atom_x >= cell_x_origin && atom_x <= cell_x_tr && atom_y >= cell_y_origin && atom_y <= cell_y_tr) {
                            atoms->setCell(i, (*cell)->GetCellID());
                            cell_atoms.push_back(i);
                        }

                        break;

                    case 3:
                        if( atom_x >= cell_x_origin && atom_x <= cell_x_tr &&
                            atom_y >= cell_y_origin && atom_y <= cell_y_tr &&
                            atom_z >= cell_z_origin && atom_z <= cell_z_tr) {

                            atoms->setCell(i, (*cell)->GetCellID());
                            cell_atoms.push_back(i);
                        }

                        break;
                }
		    }

		    (*cell)->AtomsInCell(cell_atoms);
	    )
    }
}
