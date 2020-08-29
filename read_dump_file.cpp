/*
 * read_dump_file.cpp
 *
 *  Created on: Aug 29, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "atom.h"
#include "point.h"

//Read data file

std::vector < Atom <3>* > ReadDumpFile()
{

	std::ifstream infile;
	infile.open("test_dump.txt");
	std::string line;
	int lines_read=0;
	int line_number=9;

	vector < Atom <3>* > atoms;

	while (std::getline(infile,line))
	{

	   lines_read++;

	   if (lines_read > line_number)
	   {

		   int type=0,id=0;
		   double x=0.,y=0.,z=0.;

		   std::istringstream iss(line);
		   iss>> id >> type >> x >> y >> z;

		   Point <3> *spatial_position=new Point <3>;
		   Atom <3> *atom=new Atom <3>;

		   spatial_position->SetXCoord(x);
		   spatial_position->SetYCoord(y);
		   spatial_position->SetZCoord(z);

		   atom->SetID(id);
		   atom->SetSpatialPosition(*spatial_position);
		   atom->SetAtomRegion(type);

		   atoms.push_back(atom);

		}

    }

	return atoms;

}








