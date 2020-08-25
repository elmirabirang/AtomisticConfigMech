/*
 * test_crack_atoms.cpp
 *
 *  Created on: Jun 11, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "crack_atoms.cpp"
#include "atom.h"
#include "point.h"
#include "bond.h"

#include <iostream>
#include <vector>

#include "unrelaxed_config.cpp"
#include "bond.cpp"


using namespace std;

int main()
{

    vector < Atom <2>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <2> (3, 5, 0, 1.1, 2);
    FindNeighbors(unrelax_atoms,2.0);

    typedef typename vector <Atom <2>* >::iterator At;

//    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
//    {
//
//    	int atom_id=(*atom) -> GetID();
//    	Point <2> atom_material_position=(*atom)-> GetMaterialPosition();
//    	double atom_x=atom_material_position.GetXCoord();
//    	double atom_y=atom_material_position.GetYCoord();
//
//        cout << "atom id: " << atom_id << "\n";
//        cout << "atom position x: " << atom_x << "\n";
//        cout << "atom position y: " << atom_y << "\n";
//
//    }

    Atom <2>* at=unrelax_atoms[6];

    CrackAtom <2>* crack_atom;

    crack_atom->setCrackRegion(at);

//    vector < Atom <2>* > bottom_atoms=crack_atom->getCrackBottomAtoms();
//
//    for (At atom=bottom_atoms.begin(); atom!=bottom_atoms.end(); ++atom)
//    {
//
//    	int atom_id=(*atom)->GetID();
//    	cout<< "atom id: " << atom_id <<endl;
//
//    }


return 0;

}





