/*
 * test_neighbor_cell_based.cpp
 *
 *  Created on: Mar 18, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */



#include "atom.h"
#include "cell.h"
#include <vector>

#include "generate_cells.cpp"
#include "unrelaxed_config.cpp"
#include "cell_neighbors.cpp"
#include "atoms_in_cell.cpp"
#include "neighbors_cell_based.cpp"
#include "cell.cpp"


int main ()
{

	vector < Atom <2>* > atoms=UnrelaxedConfigGenerator <2>(2, 4, 0, 1.1, 2);

	vector < Cell <2>* > cells=GenerateCells <2> (2.2, 1.1, 2, 4, 0) ;

	cout << "number of cells: "<< cells.size() <<endl;

	CellsNeighbors<2>(cells, 2.2);

	AtominCell <2> (atoms, cells, 2.2);

	AtomNeighbors<2>(cells, 2.2);

	typedef typename vector < Atom <2>* >::iterator At;
	typedef typename vector < Cell <2>* >::iterator Cel;

	for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
	{

		vector < Atom <2>*> cell_atoms=(*cell)->GetCellAtoms();

		for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
		{

			vector < Atom <2>* > atom_neighbors=(*atom)->Neighbor();
//			vector < Atom <2>* > atom_neighbors=(*atom)->BondNeighbor();
			int atom_id=(*atom)->GetID();

			cout << "atom id: " << atom_id <<endl;

			for(At atomn=atom_neighbors.begin(); atomn!=atom_neighbors.end(); ++atomn)
			{

				int atomn_id=(*atomn)->GetID();
				cout << "atomn id: " << atomn_id <<endl;

			}


		}



	}


//	vector < Atom <3>* > atoms=UnrelaxedConfigGenerator <3>(2, 4, 2, 1.1, 3);
//
//	vector < Cell <3>* > cells=GenerateCells <3> (2.2, 1.1, 2, 4, 2) ;
//
//	cout << "number of cells: "<< cells.size() <<endl;
//
//	CellsNeighbors<3>(cells, 2.2);
//
//	AtominCell <3> (atoms, cells, 2.2);
//
//	AtomNeighbors<3>(cells, 2.2);
//
//	typedef typename vector < Atom <3>* >::iterator At;
//	typedef typename vector < Cell <3>* >::iterator Cel;
//
//	for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
//	{
//
//
//		vector < Atom <3>*> cell_atoms=(*cell)->GetCellAtoms();
//
//		for (At atom=cell_atoms.begin(); atom!=cell_atoms.end(); ++atom)
//		{
//
//			vector < Atom <3>* > atom_neighbors=(*atom)->Neighbor();
//			int atom_id=(*atom)->GetID();
//
//			cout << "atom id: " << atom_id <<endl;
//
//			for(At atomn=atom_neighbors.begin(); atomn!=atom_neighbors.end(); ++atomn)
//			{
//
//				int atomn_id=(*atomn)->GetID();
//				cout << "atomn id: " << atomn_id <<endl;
//
//			}
//
//
//		}
//
//
//
//	}

	return 0;

}





