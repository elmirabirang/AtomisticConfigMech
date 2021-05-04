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
//---
#include "atom.h"
#include "point.h"

//Read data file

Atoms<3> *ReadDumpFile() {
    std::ifstream infile;
    infile.open("test_dump.txt");
    std::string line;
    int lines_read = 0;
    int line_number = 9;
    Atoms<3> *atoms = new Atoms<3>(1024);

    while (std::getline(infile,line)) {
        lines_read++;

        if (lines_read > line_number) {
            int type, id;
            double x, y, z;
            std::istringstream iss(line);
            iss >> id >> type >> x >> y >> z;

            int i = atoms->addAtom();
            atoms->setSpatialPosition(i, Point<3>(x, y, z));
            atoms->setRegion(i, type);
        }
    }

    return atoms;
}
