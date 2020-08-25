/*
 * practice_parallel.cpp
 *
 *  Created on: May 4, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <iostream>
#include "cell.h"
#include "cell.cpp"
#include "point.h"
#include <vector>
#include <omp.h>
#include "generate_cells.cpp"
#include "unrelaxed_config.cpp"
#include "atoms_in_cell.cpp"
#include "neighbors_cell_based.cpp"
#include "cell_neighbors.cpp"


using namespace std;


int main()
{

	vector < Cell <2>* > cells=GenerateCells <2> (2.2, 1.1, 12, 12, 0);
	vector < Atom <2>* > atoms=UnrelaxedConfigGenerator <2> (13, 13, 0 , 1.1, 2);

	CellsNeighbors<2>(cells, 2.2);
	AtominCell(atoms, cells, 2.2);
	AtomNeighbors<2>(cells, 2.2);

	typedef typename vector < Cell <2>* >::iterator cell;
	typedef typename vector <vector < Cell <2>*> >::iterator Cells;

	vector < Cell <2>* > cells_1;
	vector < Cell <2>* > cells_2;
	vector < Cell <2>* > cells_3;
	vector < Cell <2>* > cells_4;

	for (cell C=cells.begin(); C!=cells.end(); ++C)
	{

		int cell_id=(*C)->GetCellID();

		if (cell_id==0 || cell_id==30 || cell_id==5 || cell_id== 35)
		{

			cells_1.push_back(*C);

		}

		else if (cell_id==6 || cell_id==24 || cell_id==11 || cell_id== 29)
		{

			cells_2.push_back(*C);

		}

		else if(cell_id==1 || cell_id==4 || cell_id==31 || cell_id== 34)
		{

			cells_3.push_back(*C);

		}


		else if(cell_id==7 || cell_id==25 || cell_id==10 || cell_id==28)
		{

			cells_4.push_back(*C);

		}


	}

	vector < vector < Cell <2>* > > cells_to_update;

	cells_to_update.push_back(cells_1);
	cells_to_update.push_back(cells_2);
	cells_to_update.push_back(cells_3);
	cells_to_update.push_back(cells_4);

	typedef typename vector < Atom <2>* >::iterator At;

#pragma omp parallel

	{

      #pragma omp single

		{

          #pragma omp taskgroup

			{
				for (Cells cells=cells_to_update.begin(); cells!=cells_to_update.end(); ++cells )
				{

                   #pragma omp taskgroup
					{

						for (cell c=cells->begin(); c!=cells->end(); ++c )
						{

                           #pragma omp task

							{

								vector < Atom <2>*> cell_atoms=(*c)->GetCellAtoms();

								for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
								{

									vector < Atom <2>* > atom_neighbors=(*atom)->Neighbor();
						            //vector < Atom <2>* > atom_neighbors=(*atom)->BondNeighbor();

									int atom_id=(*atom)->GetID();
									cout << "atom id: " << atom_id <<endl;

								}

							}

                           #pragma omp taskwait

						}

					}

				}

			}

		}

	}




return 0;

}








