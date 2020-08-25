/*
 * test_boundary.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlanegn-nürnberg
 */

#include "boundary.h"
#include "boundary.cpp"

#include <iostream>

using namespace std;
int main()
{

	Boundary <3> boundary;

	boundary.SetBoundaryRegion(0,1,2,1,2,3,0);

	double x0=boundary.GetBoundaryX0();
	double y0=boundary.GetBoundaryY0();
	double z0=boundary.GetBoundaryZ0();

	double x1=boundary.GetBoundaryX1();
	double y1=boundary.GetBoundaryY1();
	double z1=boundary.GetBoundaryZ1();

	cout << "x0: " <<x0 << " " <<
			"y0: " <<y0 << " " <<
			"z0: " <<z0 << endl;

	cout << "x1: " <<x1 << " " <<
			"y1: " <<y1 << " " <<
			"z1: " <<z1 << endl;

	Boundary <2> boundary2d;
	boundary2d.SetTriangBoundaryRegion(0,0,2,1,0,2,0);


	double x0_2d=boundary2d.GetBoundaryX0();
	double y0_2d=boundary2d.GetBoundaryY0();

	double x1_2d=boundary2d.GetBoundaryX1();
	double y1_2d=boundary2d.GetBoundaryY1();

	double x2_2d=boundary2d.GetBoundaryX2();
	double y2_2d=boundary2d.GetBoundaryY2();


	cout << "x0: " <<x0_2d << " " <<
			"y0: " <<y0_2d << endl;

	cout << "x1: " <<x1_2d << " " <<
			"y1: " <<y1_2d << endl;

	cout << "x2: " <<x2_2d << " " <<
			"y2: " <<y2_2d << endl;


return 0;

}




