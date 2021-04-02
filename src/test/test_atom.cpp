/*
 * tets_atom.cpp
 *
 *  Created on: Dec 2, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical engineering PhD candidate
 *      Chiar of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "atom.h"
#include "atom.cpp"

#include "unrelaxed_config.cpp"
#include "point.h"
#include "point.cpp"

#include <vector>
#include <iostream>
#include <fstream>


using namespace std;
int main ()
{
	vector < Point <2> > atms;

	atms=UnrelaxedConfigGenerator <2> (2, 1, 0, 1.1);

	for






return 0;
}



