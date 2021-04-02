/*
 * test_energy.cpp
 *
 *  Created on: Dec 2, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */

#include "energy.cpp"
#include "energy.h"
#include "unrelaxed_config.cpp"

#include <vector>



using namespace std;
int main ()
{

	vector < Atom <2>* > atms;
	atms=UnrelaxedConfigGenerator <2> (2, 3, 0, 1.1,2);

	vector < Atom <2>* > atoms=FindNeighbors(atms,10);

	Energy <2> toteng(1.0,1.0);
	double energy=toteng.TotPotentialEnergy(atoms);

	cout << "energy is: " << energy << endl;

	vector < Atom <3>* > atms_3d;
	atms_3d=UnrelaxedConfigGenerator <3> (5, 7, 1, 1.1, 3);

	vector < Atom <3>* > atoms_3d=FindNeighbors(atms_3d, 10);

	Energy <3> toteng3d(1.0,1.0);

	double energy_3d=toteng3d.TotPotentialEnergy(atoms_3d);

	cout << "energy of 3d model is: " << energy_3d << endl;




return 0;

}



