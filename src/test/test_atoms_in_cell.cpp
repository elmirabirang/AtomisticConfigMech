/*
 * test_atoms_in_cell.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */



#include "atom.h"
#include "cell.h"

#include "unrelaxed_config.cpp"
#include "generate_cells.cpp"
#include "atoms_in_cell.cpp"
#include "cell.cpp"


int main ()
{

//	vector < Atom <2>* > atoms=UnrelaxedConfigGenerator <2> (4, 4, 0 , 1.1, 2);
//
//	vector < Cell <2>* > cells=GenerateCells <2> (2.2, 1.1, 4, 4, 0);
//	typedef typename vector < Cell <2>* >::iterator Cel;
//
//	for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
//	{
//
//		Point <2> cell_origin=(*cell)->GetOrigin();
//
//		double cell_origin_x=cell_origin.GetXCoord();
//		cout << "cell_x: "<< cell_origin_x<<endl;
//
//		double cell_origin_y=cell_origin.GetYCoord();
//		cout << "cell_y: "<< cell_origin_y<<endl;
//
//	}
//
//	cout << "cells size: "<< cells.size() << endl;
//
//	AtominCell(atoms, cells, 2.2);
//
//	typedef typename vector < Atom <2> *>::iterator At;
//	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//	{
//
//		int cell_id=(*atom)-> GetCellID();
//		cout << "cell id: " << cell_id << endl;
//
//		Point <2> material_position=(*atom)->GetMaterialPosition();
//
//		double materialPx=material_position.GetXCoord();
//		double materialPy=material_position.GetYCoord();
//
////		cout << "materialpx:"<< materialPx << "\n";
////		cout << "materialpy:"<< materialPy << endl;
//
//
//	}

	vector < Atom <3>* > atoms=UnrelaxedConfigGenerator <3> (3, 4, 2 , 1.1, 3);

	vector < Cell <3>* > cells=GenerateCells <3> (2.2, 1.1, 3, 4, 2);
	typedef typename vector < Cell <3>* >::iterator Cel;

	for (Cel cell=cells.begin(); cell!=cells.end(); ++cell)
	{

		Point <3> cell_origin=(*cell)->GetOrigin();

		double cell_origin_x=cell_origin.GetXCoord();
		cout << "cell_x: "<< cell_origin_x<<endl;

		double cell_origin_y=cell_origin.GetYCoord();
		cout << "cell_y: "<< cell_origin_y<<endl;

		double cell_origin_z=cell_origin.GetZCoord();
		cout << "cell_z: "<< cell_origin_z<<endl;

	}

	cout << "cells size: "<< cells.size() << endl;

	AtominCell(atoms, cells, 2.2);

	typedef typename vector < Atom <3> *>::iterator At;
	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	{

		int cell_id=(*atom)-> GetCellID();
		cout << "cell id: " << cell_id << endl;

		Point <3> material_position=(*atom)->GetMaterialPosition();

		double materialPx=material_position.GetXCoord();
		double materialPy=material_position.GetYCoord();
		double materialPz=material_position.GetZCoord();

//		cout << "materialpx:"<< materialPx << "\n";
//		cout << "materialpy:"<< materialPy << endl;


	}





return 0;
}




