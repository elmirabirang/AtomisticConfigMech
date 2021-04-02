/*
 * config_force_criterion.cpp
 *
 *  Created on: Jun 25, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "force.h"
#include "bond.h"
#include "atom.h"

#include <iostream>
#include <vector>

using namespace std;

template <int dim>
bool ConfigForceCriterion (vector < Atom <dim>* > atoms)
{

typedef typename vector < Atom <dim>* >::iterator At;

Bond <dim> bond;
Force <dim> force(1.0,1.0);
//double critical_config_force=4.8;
//double crack_line_y=13.9;

double critical_config_force=10.32;
double crack_line_z=33.5;

double x_max=67.;
double x_min=1.;

double y_max=16.;
double y_min=-1.0;

double z_max=35.;
double z_min=32.;

double alpha1=2.97758;
double alpha2=1.54927;

double r_01=0.83591;
double r_02=4.46867;

double E1=2.01458e2;
double E2=6.59288e-3;

double delta=0.86225e-2;
double h= 0.50037;
double r_cut=5.50679;

double r_s1=2.24;
double r_s2=1.8;
double r_s3=1.2;

double s1=4;
double s2=40;
double s3=1150;


double a= 3.80362;

double beta1= 0.17394;
double beta2= 5.3566e2;

double r_03= -2.19885;
double r_04= -261.984;

double F0= -2.28235;
double F2= 1.35535;

double q1= -1.27775 ;
double q2=-0.86074;
double q3= 1.78804;
double q4= 2.97571;

double Q1= 0.4 ;
double Q2= 0.3 ;



bool bond_broken=false;

for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)

{

	Point <dim> atom_position= (*atom)-> GetMaterialPosition();

	double atomp_x=atom_position.GetXCoord();
	double atomp_y=atom_position.GetYCoord();
	double atomp_z=atom_position.GetZCoord();

	int atom_region=(*atom)->GetAtomRegion();


	if ( atomp_x < x_max && atomp_x > x_min && atomp_y > y_min && atomp_y < y_max && atomp_z < z_max && atomp_z > z_min)
	{

	    //Point <dim> config_force =force.ConfigurationalForce(*atom);



		Point <dim> config_force_pair = force.configForceEamPair(*atom, alpha1, alpha2,
                                                               r_01, r_02, E1, E2, delta,
                                                               r_cut, h, r_s1, r_s2, r_s3,
                                                               s1, s2, s3);

		Point <dim> config_force_embedding = force.configForceEamEmbedding(*atom, a, beta1,
                                                                         beta2, r_03, r_04, F0, F2,
                                                                         q1, q2, q3, q4, Q1, Q2,
                                                                         h, r_cut);

		Point <dim>  config_force=(config_force_pair+config_force_embedding);

		double config_force_value=config_force.PointNorm();

	    if (config_force_value > critical_config_force)

	    {

	    	Atom <dim> *atm;
	    	atm=*atom;
	    	bond_broken=true;

	    	vector < Atom <dim>* > neighbors=(*atom)->Neighbor();
	    	vector < Atom <dim>* > bond_neighbors=(*atom)->BondNeighbor();

	    	vector <Atom <dim>*> update_neighbor_list;
	    	vector <Atom <dim>*> update_bond_neighbor_list;

	    	   for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	    	   {

			   	   Atom <dim> *neigh;
			   	   neigh=*neighbor;

		   		   Point <dim> neighbor_position= (*neighbor)-> GetMaterialPosition();
		   		   double neighborp_x=neighbor_position.GetXCoord();
		   		   double neighborp_y=neighbor_position.GetYCoord();
		   		   double neighborp_z=neighbor_position.GetZCoord();

		   		   if (neighborp_z > crack_line_z && atomp_z > crack_line_z)

		   		   {

		   					//remove (neighbors.begin(), neighbors.end(),*neighbor);
		   					update_neighbor_list.push_back(*neighbor);

		   		   }

		   		   else if (neighborp_z < crack_line_z && atomp_z < crack_line_z)

		   		   {

		   			       //remove (neighbors.begin(), neighbors.end(),*neighbor);
		   			       update_neighbor_list.push_back(*neighbor);

		   		   }

	    	   }

	    	   for (At bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	    	   {

			   	   Atom <dim> *bond_neigh;
			   	   bond_neigh=*bond_neighbor;

		   		   Point <dim> bond_neighbor_position= (*bond_neighbor)-> GetMaterialPosition();

		   		   double neighborp_x=bond_neighbor_position.GetXCoord();
		   		   double neighborp_y=bond_neighbor_position.GetYCoord();
		   		   double neighborp_z=bond_neighbor_position.GetZCoord();

		   		   if (neighborp_z > crack_line_z && atomp_z > crack_line_z)

		   		   {

		   					//remove (bond_neighbors.begin(), bond_neighbors.end(),*bond_neighbor);
		   					update_bond_neighbor_list.push_back(*bond_neighbor);

		   		   }

		   		   else if (neighborp_z < crack_line_z && atomp_z < crack_line_z)

		   		   {

		   			       //remove (bond_neighbors.begin(), bond_neighbors.end(),*bond_neighbor);
		   			       update_bond_neighbor_list.push_back(*bond_neighbor);

		   		   }

	    	   }

	    	   (*atom)->SetNeighbors(update_neighbor_list);
	    	   (*atom)->SetBondNeighbors(update_bond_neighbor_list);

	    }

	    else
	    {

	    	bond_broken=false;
	    }

	}

}

return (bond_broken);

}

template bool ConfigForceCriterion<3>(vector < Atom <3>* > atoms);
