/*
 * test_SW.cpp
 *
 *  Created on: Aug 3, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#include "energy.h"
#include "energy.cpp"

#include "force.h"
#include "force.cpp"

#include "point.h"
#include "atom.h"
#include "bond.h"
#include "bond.cpp"

#include "matrixx.h"
#include "Stillinger_Weber.h"

#include "third_order_tensors.h"
#include "write_data_file.cpp"
#include "unrelaxed_config.cpp"
#include "Stillinger_Weber.cpp"
#include "matrix.cpp"

#include <math.h>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <list>
#include <algorithm>
#include <iomanip> 

int main()
{

	#define epsilon 2.1683
	#define sigma 2.0951
	#define a 1.80
	#define eta 21.0
	#define mu 1.20
	#define costheta0 -0.333333
	#define A 7.049556277
	#define B 0.6022245584
	#define p 4.0
	#define q 0.0
	#define tol 0.0

	Point <3> atom_alpha_matP( 1.35775,1.35775,1.35775 );
	Point <3> atom_beta_matP(0.,0.,0.);
	Point <3> atom_gamma_matP(2.7155,0.,2.9);

	Point <3> atom_alpha_spaP(1.35775,1.35775,1.35775);
	Point <3> atom_beta_spaP(0.,0.,0.);
	Point <3> atom_gamma_spaP(2.7155,0.,2.9);

	Atom <3> *atom_alpha=new Atom <3>;
	Atom <3> *atom_beta=new Atom <3>;
	Atom <3> *atom_gamma=new Atom <3>;

	atom_alpha->SetID(0);
	atom_beta->SetID(1);
	atom_gamma->SetID(2);

	atom_alpha->SetMaterialPosition(atom_alpha_matP);
	atom_beta->SetMaterialPosition(atom_beta_matP);
	atom_gamma->SetMaterialPosition(atom_gamma_matP);

	atom_alpha->SetSpatialPosition(atom_alpha_spaP);
	atom_beta->SetSpatialPosition(atom_beta_spaP);
	atom_gamma->SetSpatialPosition(atom_gamma_spaP);

	vector < Atom <3>* > unrelax_atoms;
	unrelax_atoms.push_back(atom_alpha);
	unrelax_atoms.push_back(atom_beta);
	unrelax_atoms.push_back(atom_gamma);

	vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,2.6);

	Bond <3> bond;

	Energy <3> ENG;

	Force <3> FORC;

	StillingerWeber<3> SW;

	typedef typename vector < Atom <3>* >::const_iterator At;

	double total_energy=0;

	for(At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	{

		vector < Atom <3>* > neighbors=(*atom)->Neighbor();

		Atom <3> * atm;
		atm=*atom;

		for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
		{

			Atom <3> * neigh;
			neigh=*neighbor;

			double distance=bond.SpatialBondDistance(*atm,*neigh);

			int neigh_beta_id=neigh->GetID();

			for (At neighbor_gamma=neighbors.begin(); neighbor_gamma!=neighbors.end(); ++neighbor_gamma)
			{

				Atom <3> * neigh_gamma;
				neigh_gamma=*neighbor_gamma;

				int neigh_gamma_id=neigh_gamma->GetID();

				if( neigh_gamma_id > neigh_beta_id)
				{

				    double triplet_deform_energy=SW.StillingerWeberThreeBody_Deformational
					                             (*atm, *neigh, *neigh_gamma, epsilon, eta, costheta0);

				    Point <3> triplet_deform_force(0.,0.,0.);

				    triplet_deform_force=(*atm).GetForce();

				    double triplet_config_energy=SW.StillingerWeberThreeBody_Configurational
					                             (*atm, *neigh, *neigh_gamma, epsilon, eta, costheta0);

				    Point <3> triplet_config_force(0.,0.,0.);

				    triplet_config_force=(*atm).GetConfigForce();

				}

			}

			double pair_deform_energy=SW.StillingerWeberTwoBody_Deformational(*atm, *neigh, sigma, epsilon, A, B, p, q, a);

			double pair_config_energy=SW.StillingerWeberTwoBody_Configurational(*atm, *neigh, sigma, epsilon, A, B, p, q, a);

		}

	}

	for(At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	{

//		FORC.ResultantSWThreeBodyForce(*atom,2.0951, 2.0951, 2.1683, 1.2, 21.,-0.333333,1.8, 1.8);
//
//		Point <3> resultant_force2body(0.,0.,0.);
//
//		resultant_force2body=FORC.ResultantSWTwoBodyForce(*atom,2.0951, 2.1683,
//				                                           7.049556277, 0.6022245584,
//				                                           4.0, 0.0, 1.8);

		Point <3> atom_force=(*atom)->GetForce();

		double atom_force_x=atom_force.GetXCoord();
		double atom_force_y=atom_force.GetYCoord();
		double atom_force_z=atom_force.GetZCoord();

		cout << (*atom)->GetID() << " "
				 << " atom_force_x: " <<
				 atom_force_x <<
				 " atom_force_y: "<<
				 atom_force_y <<
				 " atom_force_z: " <<
				 atom_force_z << endl;

	}

    return 0;
}
	







