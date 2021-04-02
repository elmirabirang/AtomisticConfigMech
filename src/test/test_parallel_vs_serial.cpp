/*
 * test_parallel_vs_serial.cpp
 *
 *  Created on: May 7, 2019
 *      Author: S.Elmira Birang.O
 *      Mcehanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include <iostream>
#include "cell.h"
#include "cell.cpp"
#include "point.h"
#include "force.h"
#include "force.cpp"
#include "point.h"
#include "bond.h"
#include "bond.cpp"
#include <vector>
#include <omp.h>
#include "generate_cells.cpp"
#include "unrelaxed_config.cpp"
#include "atoms_in_cell.cpp"
#include "neighbors_cell_based.cpp"
#include "cell_neighbors.cpp"



int main ()
{

//	vector < Atom <2>* > unrelax_atoms=UnrelaxedConfigGenerator <2> (3, 5, 0, 1.1, 2);
//	vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,2.2);
//	Force <2> force(1.0,1.0);
//
//    double force_total=0.0;
//
//	typedef typename  vector < Atom <2>* >::iterator At;
//
//	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//	{
//
//		Point <2> force_atom=force.ResultantForce(*atom);
//
//		cout << "fx:" << force_atom.GetXCoord() << endl;
//		cout << "fy:" << force_atom.GetYCoord() << endl;
//
//		double atom_force=force_atom.PointNorm();
////		cout << "force: " << atom_force<<endl;
//
//		force_total=force_total+atom_force;
//
//	}
//
//	cout << "total force: "<< force_total << endl;

	///////////////////////////////////////////////////////////////////////////////

	vector < Cell <2>* > cells= GenerateCells <2> (2.2, 1.1, 2, 4, 0);
	vector < Atom <2>* > atoms=UnrelaxedConfigGenerator <2> (3, 5, 0, 1.1, 2);

	Force <2> force(1.0,1.0);

	CellsNeighbors <2> (cells,2.2);
	AtominCell (atoms, cells, 2.2);
	AtomNeighbors <2> (cells,2.2);

	typedef typename  vector < Cell <2>* >::iterator cell;

	vector < Cell <2>* > cells_1;

	for(cell C=cells.begin(); C!=cells.end(); ++C)

	{

		cells_1.push_back(*C);

	}

	vector < vector < Cell <2>* > > cells_to_update;
	cells_to_update.push_back(cells_1);

	typedef typename vector < vector <Cell <2>*> >::iterator Cells;
	typedef typename  vector < Atom <2>* >::iterator At;

#pragma omp parallel

	{

 	  double force_cells=0;

	  double forces=0.0;

      #pragma omp single

		{


                #pragma omp taskloop shared(forces,force_cells)

				for (Cells cells=cells_to_update.begin(); cells!=cells_to_update.end(); ++cells )
				{
					double force_cells=0;
					double force_cell=0;


                        #pragma omp parallel for
						for (int i=0; i < cells->size(); ++i)
						{

//                            #pragma omp task shared(forces,force_cells,force_cell)

							Cell <2> *c=new Cell<2>;
							c=cells->begin();

							{
								double force_cell=0;

								vector < Atom <2>*> cell_atoms=(*c).GetCellAtoms();


								for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
								{

//                                    #pragma omp task shared(forces,force_cells,force_cell)

									{

										Point <2> force_atom=force.ResultantForce(*atom);
										double force_atom_val=force_atom.PointNorm();

//										cout << "force:" << force_atom_val <<endl;

//										cout << "fx:" << force_atom.GetXCoord() << endl;
//										cout << "fy:" << force_atom.GetYCoord() << endl;


										force_cell=force_cell+force_atom_val;

										int atom_id=(*atom)->GetID();
										cout << "atom id: " << atom_id <<endl;

										int thread_id=omp_get_thread_num();
										cout << "thread id: " << thread_id << endl;

									}

//                                   #pragma omp taskwait

								}

//								cout << "force cell: " << force_cell << endl;

								force_cells=force_cells+force_cell;


							}

//                           #pragma omp taskwait

						}

//						cout << "forces: " << force_cells << endl;





           #pragma omp taskwait

		}

	}
	}


	///////////////////////////////////////////////////////////////////////////////





return 0;
}




