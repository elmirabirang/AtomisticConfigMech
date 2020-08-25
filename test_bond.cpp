/*
 * test_bond.cpp
 *
 *  Created on: Dec 2, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "bond.h"
#include "atom.h"
#include "atom.cpp"
#include "unrelaxed_config.cpp"
#include "point.h"
#include "point.cpp"
#include <vector>
#include <iostream>

using namespace std;
int main()
{

	vector < Atom <2> > atms;

	atms=UnrelaxedConfigGenerator <2> (2, 3, 0, 1.1,2.5);

	Point <2> distance;

	Bond <2> bond;
	distance= bond.MaterialBondVec(atms[0],atms[1]);
	double x=distance.GetXCoord();
	double y=distance.GetYCoord();

	cout<< "x: " << x << " "<< "y: " << y <<endl;


	Point <2> spatial_distance;
	spatial_distance=bond.SpatialBondVec(atms[0],atms[1]);
	double x1=spatial_distance.GetXCoord();
	double y1=spatial_distance.GetYCoord();

	cout << "x1: " << x1 << " "<< "y1: " << y1 <<endl;

	double dis;

	dis=bond.MaterialBondDistance(atms[0],atms[1]);

	cout << "material distance: " << dis << endl;

	double spatialdis;

	spatialdis=bond.SpatialBondDistance(atms[0],atms[1]);

	cout << "spatial dis: " << spatialdis  << endl;

	Point <2> materialnormal;
	materialnormal=bond.MaterialBondNormal(atms[0],atms[1]);
	double xpoint=materialnormal.GetXCoord();
	double ypoint=materialnormal.GetYCoord();
	cout << "x_normal: " << xpoint << " " << "y_normal: "<< ypoint << endl;

	Point <2> spatialnormal;
	spatialnormal=bond.SpatialBondNormal(atms[0],atms[1]);
	double xspatial=spatialnormal.GetXCoord();
	double yspatial=spatialnormal.GetYCoord();
	cout << "x_normal: " << xspatial<< " " << "y_normal" << yspatial<<endl;

	double material_stretch=bond.MaterialBondStretch(atms[0],atms[1]);
	cout <<"material stretch:" << material_stretch<<endl;

	double spatial_stretch=bond.SpatialBondStretch(atms[0],atms[1]);
	cout << "spatial stretch: " << spatial_stretch<< endl;
	
	Point <2> material_result;
	
	material_result=bond.ResultantBondVec_Material(atms[1]);
	
	cout << "xm_resultant: " << material_result.GetXCoord() 
		<< " " 
		<< "ym_resultant: " << material_result.GetYCoord() << endl;
	
	Point <2> spatial_result;
	
	spatial_result=bond.ResultantBondVec_Spatial(atms[1]);
	
	cout << "xs_resultant: " << spatial_result.GetXCoord() 
		<< " " 
		<< "ys_resultant: " << spatial_result.GetYCoord() << endl;
			
return 0;

}




