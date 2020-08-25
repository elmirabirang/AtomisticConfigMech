/*
 * test_force.cpp
 *
 *  Created on: Dec 4, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Phd Candidate
 *      University of Erlangen-Nürnberg
 */


#include "force.h"
#include "force.cpp"
#include "unrelaxed_config.cpp"

#include <vector>
#include <iostream>

using namespace std;

int main(){

	Force <2> force(1.0,1.0);

	vector < Atom <2>* > atms;
	atms=UnrelaxedConfigGenerator <2> (2, 3, 0, 1.1,2);

	vector < Atom <2>* > atoms = FindNeighbors(atms, 10);

//	Point <2> force_ij;
//
//	force_ij=force.BondForce(*atoms[0],*atoms[1]);
//
//	double fx=force_ij.GetXCoord();
//	double fy=force_ij.GetYCoord();
//
//	cout << "fx: " <<fx<< " "
//
//		 <<"fy: " <<fy<<endl;
//
//	Point <2> force_ij_mat;
//
//	force_ij_mat=force.MaterialBondForce(*atoms[0],*atoms[1]);
//	double fx_mat=force_ij_mat.GetXCoord();
//	double fy_mat=force_ij_mat.GetYCoord();
//
//	cout << "fx: " <<fx_mat<< " "
//
//		 <<"fy: " <<fy_mat<<endl;
//
//
//	Point <2> force_ij_spa;
//
//	force_ij_spa=force.SpatialBondForce(*atoms[0],*atoms[1]);
//	double fx_spa=force_ij_spa.GetXCoord();
//	double fy_spa=force_ij_spa.GetYCoord();
//
//	cout << "fx: " <<fx_spa<< " "
//
//		 <<"fy: " <<fy_spa<<endl;

	for (int i=0; i<atoms.size() ; ++i)
	{


		force.ResultantForce(atoms[i]);
//
//		double resultant_force_x=resultant_force.GetXCoord();
//		double resultant_force_y=resultant_force.GetYCoord();

		Point <2> atom_force=atoms[i] -> GetForce();

		double force_x = atom_force.GetXCoord();
		double force_y = atom_force.GetYCoord();
	//	double force_z = atom_force.GetZCoord();

		cout << "fx: " << force_x << " "<<
			 "fy: " <<force_y <<endl;

	}



//		Point <2> resultant_force;
//		resultant_force=force.ResultantForce(atoms[4]);
//
//		double resultant_force_x=resultant_force.GetXCoord();
//		double resultant_force_y=resultant_force.GetYCoord();
//
//		Point <2> atom_force=atoms[4] -> GetForce();
//
//		double force_x = atom_force.GetXCoord();
//		double force_y = atom_force.GetYCoord();
//	//	double force_z = atom_force.GetZCoord();
//
//		cout << "fx: " << resultant_force_x << " "<<
//			 "fy: " <<resultant_force_y << endl;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	Point <3> resultant_force_3d;
//
//	Force <3> force_3d(1.0,1.0);
//
//	vector < Atom <3>* > atms_3d;
//	atms_3d=UnrelaxedConfigGenerator <3> (2, 3, 2, 1.1,3);
//
//	vector < Atom <3>* > atoms_3d=FindNeighbors(atms_3d, 10);
//
//	resultant_force_3d=force_3d.ResultantForce(atoms_3d[4]);
//
//	double resultant_force_x3d=resultant_force_3d.GetXCoord();
//	double resultant_force_y3d=resultant_force_3d.GetYCoord();
//	double resultant_force_z3d=resultant_force_3d.GetZCoord();
//
//
//	cout << "fx: " << resultant_force_x3d << " "
//		 <<"fy: "   <<resultant_force_y3d << " "
//		 << "fz: " << resultant_force_z3d <<endl;


	return 0;

}

