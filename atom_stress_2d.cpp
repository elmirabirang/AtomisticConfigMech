/*
 * atom_stress_2d.cpp
 *
 *  Created on: Jan 22, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing (ZISC)
 *      University of ERlangen-Nuremberg
 *      Calculation of Atomic Stress for 2D Crystalline Lattice modelled by Lennard-Jones Potential with coefficients sigma=1.0, epsilon=1.0.
 *      Source Paper: Atomistic Simulation of J-Integral in 2D Graphene Nanosystems, Y.Jin 2005
 *      DOI: https://www.ncbi.nlm.nih.gov/pubmed/16430147
 *      Source Web-Page: https://lammps.sandia.gov/doc/compute_stress_atom.html
 *
 */


#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

using namespace std;



template <int dim>
void AtomicStress(vector <Atom <dim>* > atoms)
{
	typedef typename vector < Atom <2>* >::const_iterator At;
	Bond <dim> bond;
	Force <dim> force(1.0,1.0);

    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	vector <Atom <dim>*> neighbors= (*atom)->Neighbor();

    	double stress_11=0.0;
    	double stress_12=0.0;
    	double stress_21=0.0;
    	double stress_22=0.0;
    	double vonMisesStress=0.;
    	double vonMisesStressSquared=0.;

		Atom <2> *atm;
		atm=*atom;

    	for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
    	{

			Atom <2> *neigh;
			neigh=*neighbor;

    		Point <dim> InteratomicSpatialVec=bond.SpatialBondVec(*atm,*neigh);

    		double InteratomicSpatialVec_X=InteratomicSpatialVec.GetXCoord();
    		double InteratomicSpatialVec_Y=InteratomicSpatialVec.GetYCoord();

    		Point <dim> InteratomicForce=force.BondForce(*atm,*neigh);

    		double InteratomicForce_X=InteratomicForce.GetXCoord();
    		double InteratomicForce_Y=InteratomicForce.GetYCoord();

    		stress_11=stress_11+InteratomicSpatialVec_X*InteratomicForce_X;
    		stress_12=stress_12+InteratomicSpatialVec_X*InteratomicForce_Y;
    		stress_21=stress_21+InteratomicSpatialVec_Y*InteratomicForce_X;
    		stress_22=stress_22+InteratomicSpatialVec_Y*InteratomicForce_Y;

    	}

    	vector <double> AtomicStress;

    	AtomicStress.push_back(stress_11);
    	AtomicStress.push_back(stress_12);
    	AtomicStress.push_back(stress_21);
    	AtomicStress.push_back(stress_22);

    	vonMisesStressSquared=pow(stress_11,2)+pow(stress_22,2)+6*pow(stress_12,2)-2*(stress_11*stress_22);
    	vonMisesStress=sqrt(vonMisesStressSquared);

    	AtomicStress.push_back(vonMisesStress);

    	(*atom)->setAtomicStress(AtomicStress);

    }

}


