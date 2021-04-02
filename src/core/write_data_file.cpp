/*
 * write_data_file.cpp
 *
 *  Created on: Jul 8, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nuremberg
 */

#include <vector>
#include "atom.h"
#include "point.h"
#include "force.h"
#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

using namespace std;

template <int dim>
void writeDataFile (vector <Atom <dim>* > atoms, int load_step, string path)
{

	Force<dim> force;

    string out_string;
    stringstream ss;
    ss << load_step;
    out_string=ss.str();

    string file_name=path+out_string;
	ofstream Tfile;
	Tfile.open(file_name.c_str());

	int id=0;

	Tfile << "LAMMPS Description" << endl;
	Tfile << "\n" << endl;
	Tfile << "     "<< atoms.size() << "  " << "atoms"<<endl;
	// In advance number of bonds is not known.

	Tfile << "     "<< "7489" << "  "<< "bonds"<<endl;

	Tfile << "     "<< "0" << "  " << "angles"<<endl;
	Tfile << "     "<< "0" << "  " << "dihedrals"<<endl;
	Tfile << "     "<< "0" << "  " << "impropers"<<endl;
	Tfile << "\n" << endl;
	Tfile << "     "<< "1" << " "<< "atom" << " " << "types" <<endl;
	Tfile << "     "<< "1" << " "<< "bond" << " " << "types" <<endl;
	Tfile << "\n" << endl;
	Tfile << "  "<< "-1" << " "<< "40" << " " << "xlo" << " "<< "xhi" << endl;
	Tfile << "  "<< "-1" << " "<< "40" << " " << "ylo" << " "<< "yhi" << endl;
	Tfile << "  "<< "-1" << " "<< "40" << " " << "zlo" << " "<< "zhi" << endl;
	Tfile << endl;
	Tfile << "Atoms " << endl;
	Tfile << endl;

	typedef typename vector <Atom <dim>*> :: iterator At;

	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)

	{

		int atom_id=(*atom)-> GetID();
		Point <dim> spatial_position = (*atom)-> GetSpatialPosition();

		double spatialp_x=spatial_position.GetXCoord();
		double spatialp_y=spatial_position.GetYCoord();
		double spatialp_z=spatial_position.GetZCoord();

		// in data file can not have an atom from type 0
		//"Add this command when lattice structure has several regions."(*atom)->GetAtomRegion()+1

		  Tfile << id << " " << "1"<< " " << "1" << " " << "1.0" << " "<< spatialp_x
		  << " " << setprecision(5)<< spatialp_y <<" " << spatialp_z
		  <<endl;
		   id+=1;

	}

	Tfile << endl;
	Tfile << "Bonds" << endl;
	Tfile << endl;
	int count_num_neighbors=0;

	for (At at=atoms.begin(); at!=atoms.end(); ++at)
	{

		vector < Atom <dim>* > neighbors= (*at)-> Neighbor();

		for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
		{

			count_num_neighbors+=1;
			Tfile << count_num_neighbors << " " << "1" << " " << (*at)->GetID() << " "  << (*neighbor)->GetID() << endl;

		}

	}
	Tfile.close();

}

template void writeDataFile<2>(vector <Atom <2>* > atoms, int load_step, string path);
template void writeDataFile<3>(vector <Atom <3>* > atoms, int load_step, string path);
