/*
 * test_dumpfile.cpp
 *
 *  Created on: Aug 29, 2020
 *      Author: elmira
 */

#include "read_dump_file.cpp"
#include <vector>
#include <iostream>
#include "atom.h"
#include "atom.cpp"
#include "point.h"
#include "point.cpp"


int main()
{

	std::vector < Atom <3>* > unrelax_atoms=ReadDumpFile();
   	ofstream MFile;
   	string MPath="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/dump.test";
   	string Mfile_name=MPath;
   	MFile.open(Mfile_name.c_str());

   	typedef typename vector < Atom <3>* >::iterator At;
   	if ( MFile.is_open())
   	{
   		cout << "here" << "\n";

   		MFile <<"ITEM: TIMESTEP" <<endl;
   		MFile <<  0 <<endl;
   		MFile <<"ITEM: NUMBER OF ATOMS" << endl;
   		MFile << 2673 << endl;
   		MFile <<"ITEM: BOX BOUNDS ss pp ss"<<endl;
   		MFile << "0" << " " << "1.6" << endl;
   		MFile << "0" << " " << "1.6" << endl;
   		MFile << "0" << " " << "1.6" << endl;
   		MFile <<"ITEM: ATOMS id type x y z" <<endl;

   		for(At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
   		{
		Point <3> material_position=(*atom)->GetSpatialPosition();


		double materialp_x=material_position.GetXCoord();
		double materialp_y=material_position.GetYCoord();
		double materialp_z=material_position.GetZCoord();

		int atom_id=(*atom)->GetID();

		MFile << atom_id << " " << (*atom)->GetAtomRegion() <<" " << materialp_x<<" "
			  << materialp_y <<" " << materialp_z<<endl;

   	   }
   	}

   	MFile.close();
    return 0;

}






