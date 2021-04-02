/*
 * test_generate_cells.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "atom.h"
#include "cell.h"
#include "point.h"

#include "generate_cells.cpp"
#include "cell.cpp"

#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

int main ()
{

//	vector < Cell <2>* > cells_2d=GenerateCells <2> (2.2, 1.1, 5, 7, 0);
//
//	typedef typename vector < Cell <2>* >::iterator Cel;
//
//	for (Cel cell=cells_2d.begin(); cell!=cells_2d.end(); ++cell)
//
//	{
//		Point <2> cell_origin=(*cell) -> GetOrigin();
//
//		double cell_x=cell_origin.GetXCoord();
//		double cell_y=cell_origin.GetYCoord();
//
//		cout <<"cell x: " <<cell_x <<"\n";
//		cout <<"cell y: " <<cell_y <<endl;
//
//
//	}

	vector < Cell <3>* > cells_3d=GenerateCells <3> (2.2, 1.1,5,7,4);

	typedef typename vector < Cell <3>* >::iterator Cel3d;

	for (Cel3d cell=cells_3d.begin(); cell!=cells_3d.end(); ++cell)
	{

		Point <3> cell_origin_3d=(*cell) -> GetOrigin();

		double cell_x=cell_origin_3d.GetXCoord();
		double cell_y=cell_origin_3d.GetYCoord();
		double cell_z=cell_origin_3d.GetZCoord();

		cout <<"cell x: " <<cell_x <<"\n";
		cout <<"cell y: " <<cell_y <<"\n";
		cout <<"cell z: " <<cell_z <<endl;

	}

	return 0;

}



