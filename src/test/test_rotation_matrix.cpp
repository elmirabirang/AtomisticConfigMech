/*
 * test_rotation_matrix.cpp
 *
 *  Created on: Aug 2, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      University of Erlangen-Nuremberg
 */

#include "atom.h"
#include "point.h"
#include "matrix.h"
#include "bond.h"
#include "third_order_tensor.h"

#include "bond.cpp"
#include "atom.cpp"
#include "point.cpp"
#include "third_order_tensor.cpp"

#include <iostream>

int main()
{
	Point <3> atom_alpha_material_position(0,0,0);
	Point <3> atom_beta_material_position(0,1.,0);

	Point <3> atom_alpha_spatial_position(0,0,0);
	Point <3> atom_beta_spatial_position(0,1.,0);

	Atom <3> atom_alpha;
	atom_alpha.SetID(0);

	atom_alpha.SetMaterialPosition(atom_alpha_material_position);
	atom_alpha.SetSpatialPosition(atom_alpha_spatial_position);

	Atom <3> atom_beta;
	atom_beta.SetID(1);

	atom_beta.SetMaterialPosition(atom_beta_material_position);
	atom_beta.SetSpatialPosition(atom_beta_spatial_position);

	Bond <3> bond;

	SecondTensor <double> rotation_mat(3,3,0.);

	rotation_mat=bond.RotationMatrix(atom_alpha, atom_beta);

	ThirdTensor <double> R(3,3,3,0.);

//	R=bond.DerivRotMat_RotAxis(atom_alpha, atom_beta);
	R=bond.DerivRotMat_RotAngle(atom_alpha, atom_beta);


    return 0;

}



