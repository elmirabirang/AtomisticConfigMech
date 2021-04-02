/*
 * test_cell_neighbors.cpp
 *
 *  Created on: Mar 17, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical ENGINEERING Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "cell_neighbors.cpp"
#include "cell.h"
#include "atom.h"
#include "point.h"

#include "generate_cells.cpp"
#include "cell.cpp"
#include "point.cpp"

#include <vector>
#include <iostream>

int main(){

//	vector < Cell <2>* > cells=GenerateCells <2> (2.2, 1.1, 4, 4, 0);
//
//	CellsNeighbors(cells, 2.2);

//	typedef typename vector <Cell <2> *>::iterator Cel;
//	for (Cel cell=cells.begin(); cell!=cells.end(); ++cell )
//	{
//		cout << "here" << endl;
//
//		vector < Cell <2> *> cell_neighbors=(*cell)->GetCellNeighbors();
//		cout << cell_neighbors.size()<< endl;
//
//		for (Cel celln=cell_neighbors.begin(); celln!=cell_neighbors.end(); ++celln)
//		{
//			cout << "here1" << endl;
//			cout << "neighbor cell id: " << (*celln)->GetCellID() << endl;
//
//		}
//	}


		vector < Cell <3>* > cells=GenerateCells <3> (2.2, 1.1, 4, 4, 2);

		CellsNeighbors(cells, 2.2);

		typedef typename vector <Cell <3> *>::iterator Cel;
		for (Cel cell=cells.begin(); cell!=cells.end(); ++cell )
		{

			vector < Cell <3> *> cell_neighbors=(*cell)->GetCellNeighbors();
			cout << cell_neighbors.size()<< endl;

			for (Cel celln=cell_neighbors.begin(); celln!=cell_neighbors.end(); ++celln)
			{
				cout << "neighbor cell id: " << (*celln)->GetCellID() << endl;

			}


	}

	return 0;
}



