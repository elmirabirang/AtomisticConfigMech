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



template <int dim>
void AtominCell(vector < Atom < dim>* > &atoms, vector < Cell <dim>* > &cells, double r_cut)
{

	typedef typename vector < Atom < dim>* >::iterator At;
	typedef typename vector < Cell < dim>* >::iterator Cel;


	for ( Cel cell=cells.begin(); cell!=cells.end(); ++cell)
	{
		vector < Atom <dim>*> cell_atoms;

		for(At atom=atoms.begin(); atom!=atoms.end(); ++atom )
		{


			switch (dim)
			{
			case (2):{

				int cell_id = (*atom)->GetCellID();
				if (cell_id==-1)
				{

					Point <dim> material_position=(*atom)->GetMaterialPosition();

					double atom_x=material_position.GetXCoord();
					double atom_y=material_position.GetYCoord();

					Point <dim> cell_origin=(*cell) -> GetOrigin();

					Point <dim> cell_tr=(*cell)-> GetTopRight();

					double cell_x_origin=cell_origin.GetXCoord();
					double cell_y_origin=cell_origin.GetYCoord();

					double cell_x_tr=cell_tr.GetXCoord();
					double cell_y_tr=cell_tr.GetYCoord();

					if ( atom_x >=cell_x_origin && atom_x <= cell_x_tr
						 && atom_y>= cell_y_origin && atom_y<=cell_y_tr)
					{

						int cell_id=(*cell)->GetCellID();

						(*atom)->SetCellID(cell_id);
						cell_atoms.push_back(*atom);


					}

				}


				break;
			}


			case (3):{

				int cell_id = (*atom)->GetCellID();
				if (cell_id==-1)
				{

				Point <dim> material_position= (*atom) -> GetMaterialPosition();

				double atom_x=material_position.GetXCoord();
				double atom_y=material_position.GetYCoord();
				double atom_z=material_position.GetZCoord();

				Point <dim> cell_origin=(*cell) -> GetOrigin();
				Point <dim> cell_tr=(*cell)-> GetTopRight();


				double cell_x_origin=cell_origin.GetXCoord();
				double cell_y_origin=cell_origin.GetYCoord();
				double cell_z_origin=cell_origin.GetZCoord();

				double cell_x_tr=cell_tr.GetXCoord();
				double cell_y_tr=cell_tr.GetYCoord();
				double cell_z_tr=cell_tr.GetZCoord();

				if ( atom_x >= cell_x_origin && atom_x <=cell_x_tr
					 && atom_y>=cell_y_origin && atom_y <= cell_y_tr
					 && atom_z>=cell_z_origin && atom_z <=cell_z_tr)
				{

					int cell_id=(*cell)->GetCellID();

					(*atom)->SetCellID(cell_id);
					cell_atoms.push_back(*atom);

				}
				}

				break;


			}

			}


		}
		(*cell)->AtomsInCell(cell_atoms);

	}


}




