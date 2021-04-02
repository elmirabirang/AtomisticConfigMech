/*
 * test_unrelaxed_config.cpp
 *
 *  Created on: Nov 30, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "unrelaxed_config.cpp"
#include "point.h"
#include "bond.h"
//#include "point.cpp"
#include "atom.h"
//#include "atom.cpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>

using namespace std;

int main()
{

//the parameter dim is inserted, non-type template should get its value in compilation time.

vector < Atom <2>* > atms;

atms=UnrelaxedConfigGenerator <2> (10, 11, 0, 2.5,2);


//vector < Atom <3>* > atoms_3d ;
//
//atoms_3d=UnrelaxedConfigGenerator <3>(5, 5, 2, 1.1, 3);


//for (int i=0; i<atms.size(); ++i)
//{
//    Point <2> p;
//    Point <2> p1;
//    int id;
//
//	p=atms[i].GetMaterialPosition();
//	p1=atms[i].GetSpatialPosition();
//	id=atms[i].GetID();
//
//	vector < Atom <2> > neighbors;
//	neighbors=atms[i].Neighbor();
//
//	cout << "siz: " << neighbors.size() << endl;
//
//	typedef typename vector < Atom <2> >::iterator Nghbr;
//
//	for (Nghbr nghbr=neighbors.begin(); nghbr!=neighbors.end(); ++nghbr)
//
//	{
//		cout << "atom " << atms[i].GetID()
//			 << " is neighbors with: "
//			 << (*nghbr).GetID() <<endl;
//
//	}
//
//
//
//	double x=p.GetXCoord();
//	double y=p.GetYCoord();
//
//	double x1=p1.GetXCoord();
//	double y1=p1.GetYCoord();
//
//	cout << "material position: (" << x << ", " << y << ")"<<endl;
//	cout << "spatial position: (" << x1 << ", " << y1 << ")"<<endl;
//	cout << "id: " << id << endl;
//
//}
return 0;

}



