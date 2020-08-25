/*
 * test_point.cpp
 *
 *  Created on: Nov 30, 2018
 *      Author: S.Elmira Birang.O
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "point.h"
#include "point.cpp"

#include <iostream>

using namespace std;
int main()
{

   Point <2> P1(1,0);
   Point <2> P0(4,1);
   Point <3> P2(1,0,1);

//   P.SetXCoord(0);
//   P.SetYCoord(10);

   double x1=P1.GetXCoord();
   double y1=P1.GetYCoord();

   double x2=P2.GetXCoord();
   double y2=P2.GetYCoord();
   double z2=P2.GetZCoord();

   Point <2> point;
   point=P1-P0;
   double x=point.GetXCoord();
   double y=point.GetYCoord();

   cout << "x: " << x <<", "<<"y: "<< y <<endl;

   cout << "x1: " << x1 << endl;
   cout << "y1: " << y1 << endl;

   cout << "x2: " << x2 << endl;
   cout << "y2: " << y2 << endl;
   cout << "z2: " << z2 << endl;

   Point <3> po(4,3,0);
   int dimension=po.GetDim();
   cout << "dimension: " << dimension << endl;

   double pointnorm=po.PointNorm();
   cout << "norm: " << pointnorm << endl;


	Point <2> point1(-1,0);
	double fac=2;
	Point <2> pointfac=point1*fac;
	double xp=pointfac.GetXCoord();
	double yp=pointfac.GetYCoord();

	cout << "x: " << xp <<" " << "y: " << yp<< endl;

   return 0;

}

