/*
 * test_ghost_atoms.cpp
 *
 *  Created on: Sep 16, 2020
 *      Author: S.Elmira Birang O.
 *      Mechanical Engineering PhD candidate
 *      University of Erlangen-Nürnberg
 */
#include "unrelaxed_config.cpp"
#include "atom.h"
#include "point.h"
#include "energy.h"
#include "bond.h"
#include "bond.cpp"
#include "energy.cpp"
#include "force.h"
#include "force.cpp"
#include "matrixx.h"
#include "matrix.cpp"

#include <vector>
#include <iostream>

int main()
{

	vector < Atom <3>* > atoms;
	atoms=UnrelaxedConfigGenerator <3> (1, 1, 1, 5.431, 3, "Si");

	vector <double> box_bounds;

	box_bounds.push_back(0);
	box_bounds.push_back(2*5.4310);

	box_bounds.push_back(0);
	box_bounds.push_back(1*5.4310);

	box_bounds.push_back(0);
	box_bounds.push_back(1*5.4310);

	vector <int> periodic_bc;
	periodic_bc.push_back(0);

	double GridSize_x=1*5.431;

	typedef typename vector < Atom <3>* >::const_iterator At;

    Force <3> force;

    Energy <3> energy;

	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	{

       (*atom)->SetAtomRegion(0);

	}

	vector < Atom <3>* > modified_atoms=FindNeighbors(atoms, 3.8);

	for (At atom=modified_atoms.begin(); atom!=modified_atoms.end(); ++atom)
	{

		 force.ResultantSWThreeBodyForce(*atom,2.0951, 2.0951,
                                         2.1683, 1.2, 21.,-0.333333,
                                         1.8, 1.8);

	}

	double total_energy=0.;

	for(At atom=modified_atoms.begin(); atom!=modified_atoms.end(); ++atom)
	{

			total_energy+=energy.totalStillingerWeberEnergy(*atom, 2.0951, 2.0951,
															2.1683, 1.2, 21.,-0.333333,
															1.8, 1.8,
															7.049556277, 0.6022245584,
															4.0, 0.0);

	}

	cout << "energy: " << total_energy << endl;

	string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/dump.periodic_Si";
	ofstream file;
	file.open(path.c_str());

	file <<"ITEM: TIMESTEP" <<endl;
	file << "0" <<endl;
	file <<"ITEM: NUMBER OF ATOMS" << endl;
	file << modified_atoms.size() << endl;
	file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
	file << "0" << " " << "1.6" << endl;
	file << "0" << " " << "1.6" << endl;
	file << "0" << " " << "1.6" << endl;
	file <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

	for(At atom = modified_atoms.begin(); atom < modified_atoms.end(); ++atom)
	{

		Point <3> material_position;

		material_position=(*atom)->GetMaterialPosition();

		double material_position_x=material_position.GetXCoord();
		double material_position_y=material_position.GetYCoord();
		double material_position_z=material_position.GetZCoord();

		int atom_id=(*atom)->GetID();
		int atom_region=(*atom)->GetAtomRegion();

		Point <3> resultant_force2body(0.,0.,0.);

 		resultant_force2body=force.ResultantSWTwoBodyForce(*atom,2.0951,
				                                            2.1683,
				                                            7.049556277, 0.6022245584,
				                                            4.0, 0.0, 1.8);
 		Point <3> resultatnt_force3body=(*atom)->GetForce();
//
		file << atom_id << " "
			  << atom_region <<" "
			  << material_position_x <<" "
			  << material_position_y <<" "
			  << material_position_z <<
			  " " << -resultant_force2body.GetXCoord()+resultatnt_force3body.GetXCoord()
			  << " "<<-resultant_force2body.GetYCoord()+resultatnt_force3body.GetYCoord()
			  << " "<< -resultant_force2body.GetZCoord()+resultatnt_force3body.GetZCoord()
			  << endl;

//			file << atom_id+1 << " "
//				  << atom_id+1 <<" "
//				  << "1"
//				  << " "
//				  << "1.0"
//				  << " "
//				  << material_position_x <<" "
//				  << material_position_y <<" "
//				  << material_position_z
//				  << endl;

		(*atom)->SetAtomRegion(0);

	}


}




