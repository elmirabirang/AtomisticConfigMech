/*
 * neighbors_cell_based.cpp
 *
 *  Created on: Mar 17, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Eralnegn-Nuremberg
 */

#include "atom.h"
#include "cell.h"
#include "point.h"
#include "math.h"

#include <vector>
#include <iostream>

template <int dim>
void AtomNeighbors(vector < Cell <dim> *> &cells, double r_cut)
{

	typedef typename vector < Cell <dim>* >::iterator Cel;
	typedef typename vector < Atom <dim>* >::iterator At;

	switch(dim){
	case (2):
	{

		for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
		{
			vector < Atom <dim>* > cell_atoms=(*cell)->GetCellAtoms();

			vector <Cell <dim>*>cell_neighbors=(*cell)->GetCellNeighbors();

			for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
			{
				Point <dim> atom_materialP=(*atom)->GetMaterialPosition();

				double atom_x=atom_materialP.GetXCoord();
				double atom_y=atom_materialP.GetYCoord();

				vector <Atom <dim>*> neighbors;

				int atom_id=(*atom)->GetID();

				for (Cel celln=cell_neighbors.begin();celln!=cell_neighbors.end(); ++celln)
				{

					vector < Atom <dim>* > neighbor_cell_atoms=(*celln)->GetCellAtoms();



					for (At atomn=neighbor_cell_atoms.begin(); atomn!=neighbor_cell_atoms.end(); ++atomn)
					{

						Point <dim> neighbor_atom_materialP=(*atomn)->GetMaterialPosition();

						double neighbor_atom_x=neighbor_atom_materialP.GetXCoord();
						double neighbor_atom_y=neighbor_atom_materialP.GetYCoord();

						int atomn_id=(*atomn)->GetID();

						vector <Atom <dim>*> atomn_neighbors=(*atomn)->BondNeighbor();

						double x_coord=neighbor_atom_x-atom_x;
						double y_coord=neighbor_atom_y-atom_y;

						double distance = sqrt ((x_coord*x_coord)+(y_coord*y_coord));

						if ( distance < r_cut && atom_id <atomn_id)
						{
							neighbors.push_back(*atomn);

							atomn_neighbors.push_back(*atom);
							(*atomn)->SetBondNeighbors(atomn_neighbors);

						}



					}

				}

				(*atom)->SetNeighbors(neighbors);

			}

		}

		break;

	}

	case (3):
	{
		for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
		{
			vector < Atom <dim>* > cell_atoms=(*cell)->GetCellAtoms();

			vector <Cell <dim>*>cell_neighbors=(*cell)->GetCellNeighbors();

			for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
			{
				Point <dim> atom_materialP=(*atom)->GetMaterialPosition();

				double atom_x=atom_materialP.GetXCoord();
				double atom_y=atom_materialP.GetYCoord();
				double atom_z=atom_materialP.GetZCoord();

				vector <Atom <dim>*> neighbors;
				int atom_id=(*atom)->GetID();

				for (Cel celln=cell_neighbors.begin();celln!=cell_neighbors.end(); ++celln)
				{

					vector < Atom <dim>* > neighbor_cell_atoms=(*celln)->GetCellAtoms();

					for (At atomn=neighbor_cell_atoms.begin(); atomn!=neighbor_cell_atoms.end(); ++atomn)
					{

						Point <dim> neighbor_atom_materialP=(*atomn)->GetMaterialPosition();

						double neighbor_atom_x=neighbor_atom_materialP.GetXCoord();
						double neighbor_atom_y=neighbor_atom_materialP.GetYCoord();
						double neighbor_atom_z=neighbor_atom_materialP.GetZCoord();

						int atomn_id=(*atomn)->GetID();

						vector <Atom <dim>*> atomn_neighbors=(*atomn)->BondNeighbor();

						double distance= pow ((neighbor_atom_x-atom_x),2) +
								         pow ((neighbor_atom_y-atom_y),2) +
										 pow((neighbor_atom_z-atom_z),2);

						if ( distance < pow(r_cut,2) && atom_id <atomn_id)
						{
							neighbors.push_back(*atomn);

							atomn_neighbors.push_back(*atom);
							(*atomn)->SetBondNeighbors(atomn_neighbors);


						}


					}

				}

				(*atom)->SetNeighbors(neighbors);

			}

		}

		break;

	}

	}

}





