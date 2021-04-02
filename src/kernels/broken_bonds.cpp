/*
 * broken_bonds.cpp
 *
 *  Created on: Jun 18, 2019
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
bool brokenBonds (vector < Atom <dim>* > atoms)
{

typedef typename vector < Atom <dim>* >::iterator At;

Bond <dim> bond;
Force <dim> force(1.0,1.0);
int count =0;
bool bond_broken=false;



for (At atom=atoms.begin(); atom!=atoms.end(); ++atom )

  {

	int atom_region=(*atom)->GetAtomRegion();
	Point <dim> atom_position= (*atom)-> GetMaterialPosition();
	double atomp_x=atom_position.GetXCoord();
	double atomp_y=atom_position.GetYCoord();

	int atom_id=(*atom)->GetID();

	if ( atomp_x <40.0 && atomp_x>11.4 && atom_region!=1 && atom_region!=2 && atomp_y>13.2 && atomp_y<14.3 )
	{

	Atom <dim> *atm;
	atm=*atom;

	vector < Atom <dim>* > neighbors=(*atom)->Neighbor();
	vector < Atom <dim>* > bond_neighbors=(*atom)->BondNeighbor();

	vector <Atom <dim>*> update_neighbor_list;
	vector <Atom <dim>*> update_bond_neighbor_list;

	   for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	   {

		   	    Atom <dim> *neigh;
		   	    neigh=*neighbor;

				Point <dim> config_interaction_force;
				config_interaction_force = force.MaterialBondForce(*atm,*neigh);
//				scale the pairwise configurational force
				double config_force_norm=config_interaction_force.PointNorm();


				Point <dim> critical_config_force;
//				scale critical configurational force
				critical_config_force=force.CriticalMaterialBondForce(*atm,*neigh,1.5);
				double critical_config_force_norm=critical_config_force.PointNorm()-0.000034;


				cout << critical_config_force_norm << "   " << config_force_norm << endl;


		   	 if ( config_force_norm < critical_config_force_norm)

				{
					//cout << critical_config_force_norm << "   " << config_force_norm << endl;

					update_neighbor_list.push_back(*neighbor);
					bond_broken=false;


				}

		   	 else if(config_force_norm > critical_config_force_norm) {
//
////		   		    cout << "atom id: " <<  atom_id <<endl;
//
		   		   bond_broken=true;
//					vector <Atom <dim>*> neighbor_bond_neighbors=(*neighbor)-> BondNeighbor();
//
//					remove (neighbor_bond_neighbors.begin(), neighbor_bond_neighbors.end(),*atom);
//					(*neighbor) -> SetBondNeighbors(neighbor_bond_neighbors);
//
		   	 }


	   }



	   for (At bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	   {

		   	    Atom <dim> *neigh;
		   	    neigh=*bond_neighbor;

				Point <dim> config_interaction_force;
				config_interaction_force = force.MaterialBondForce(*atm,*neigh);
//				scale the pairwise configurational force
				double config_force_norm=config_interaction_force.PointNorm();


				Point <dim> critical_config_force;
//				scale critical configurational force
				critical_config_force=force.CriticalMaterialBondForce(*atm,*neigh,1.5);
				double critical_config_force_norm=critical_config_force.PointNorm()-0.000034;


				cout << critical_config_force_norm << "   " << config_force_norm << endl;


		   	 if ( config_force_norm < critical_config_force_norm)

				{

					update_bond_neighbor_list.push_back(*bond_neighbor);
					bond_broken=false;

				}

		   	 else if(config_force_norm > critical_config_force_norm)
		   	 {

		   	        bond_broken=true;
//
//					vector <Atom <dim>*> neighbor_neighbors=(*bond_neighbor)-> Neighbor();
//
//					remove (neighbor_neighbors.begin(), neighbor_neighbors.end(),*atom);
//					(*bond_neighbor) -> SetNeighbors(neighbor_neighbors);
//
		   	 }


	   }



	   if (update_neighbor_list.size()!= 0 && update_neighbor_list.size()!=neighbors.size() )

	   {

		   (*atom)->SetNeighbors(update_neighbor_list);

		   cout <<update_neighbor_list.size() << "   " << neighbors.size() <<"\n";

		   cout << "neighbor out" << endl;

	   }

	   if (update_bond_neighbor_list.size()!= 0 && update_bond_neighbor_list.size()!=bond_neighbors.size() )

	   {

		   (*atom)->SetBondNeighbors(update_bond_neighbor_list);

		   cout <<update_bond_neighbor_list.size() << "   " << bond_neighbors.size() <<"\n";

		   cout << "bond neighbor out" << endl;

	   }

	}

	}
return (bond_broken);


}



