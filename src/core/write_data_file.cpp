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

template<int dim> void writeDataFile(Atoms<dim> *atoms, int load_step, string path) {
	Force<dim> force(atoms, 0.0, 0.0);
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
	Tfile << "     "<< atoms->getNumberOfAtoms() << "  " << "atoms"<<endl;
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

    LOOP_OVER_ATOMS(atoms, i,
		Point<dim> spatial_position = atoms->getSpatialPosition(i);
		double spatialp_x=spatial_position.GetXCoord();
		double spatialp_y=spatial_position.GetYCoord();
		double spatialp_z=spatial_position.GetZCoord();

		// in data file can not have an atom from type 0
		//"Add this command when lattice structure has several regions."(*atom)->GetAtomRegion()+1

		  Tfile << i << " " << "1"<< " " << "1" << " " << "1.0" << " "<< spatialp_x
		  << " " << setprecision(5)<< spatialp_y <<" " << spatialp_z
		  <<endl;
	)

	Tfile << endl;
	Tfile << "Bonds" << endl;
	Tfile << endl;
	int count_num_neighbors=0;

    LOOP_OVER_ATOMS(atoms, i,
        LOOP_OVER_ATOM_HALF_NEIGHBORS(atoms, i, j,
			count_num_neighbors+=1;
			Tfile << count_num_neighbors << " " << "1" << " " << i << " "  << j << endl;
        )
    )

	Tfile.close();

}

template void writeDataFile<2>(Atoms<2> *atoms, int load_step, string path);
template void writeDataFile<3>(Atoms<3> *atoms, int load_step, string path);
