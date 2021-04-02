/*
 * test_LennardJones.cpp
 *
 *  Created on: Oct 22, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "bond.h"
#include "LennardJones.cpp"
#include "unrelaxed_config.cpp"
#include "bond.cpp"

using namespace std;
int main()
{

	LennardJones <2> LJ;

    vector < Atom <2>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <2> (10, 10, 0, 1.1, 2, "model");

    typedef typename vector < Atom <2>* >::iterator At;

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {

        (*atom)->SetAtomRegion(0);

    }

    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,1.5);

    double total_energy=0;

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {

    	vector < Atom <2>* > neighbors=(*atom)->Neighbor();

   	    Atom <2> *atm;
   	    atm=*atom;

    	for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
    	{

       	    Atom <2> *neigh;
       	    neigh=*neighbor;

    		double interatomic_energy=LJ.LennardJones_Configurational(*atm, *neigh, 1.0, 1.0);

    		total_energy=total_energy+interatomic_energy;

    	}

    }

    cout << "Total Energy:" << total_energy << endl;

    return 0;

}
