/*
 * config_force_criterion_lj.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
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
vector < Atom <dim>* > ConfigForceCriterion (vector < Atom <dim>* > atoms)
{
	typedef typename vector < Atom <dim>* >::iterator At;

	Bond <dim> bond;
	Force <dim> force(1.0,1.0);

		for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)

		{
	    	Atom <2> *atm;
	    	atm=*atom;

	    	vector <Atom <2>*> update_neighbor_list;
	    	vector <Atom <2>*> update_bond_neighbor_list;

	    	Point <2> atomMaterialPosition= atm->GetMaterialPosition();

	    	double atomMaterialPositionX=atomMaterialPosition.GetXCoord();
	    	double atomMaterialPositionY=atomMaterialPosition.GetYCoord();

	    	double crackRegion_MinX=14.6;
	    	double crackRegion_MinY=22.86;

	    	double crackRegion_MaxX=58.;
	    	double crackRegion_MaxY=23.83;

			if ( atomMaterialPositionX >= crackRegion_MinX  && atomMaterialPositionY >= crackRegion_MinY
				 && atomMaterialPositionX <= crackRegion_MaxX  && atomMaterialPositionY <= crackRegion_MaxY  )
			{
	    	vector < Atom <2>* > neighbors=atm->Neighbor();
	    	vector < Atom <2>* > bond_neighbors=atm->BondNeighbor();

	    	vector <Atom <2>*> update_neighbor_list;
	    	vector <Atom <2>*> update_bond_neighbor_list;



	    	for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	    	{

	        	Point <2> neighborMaterialPosition= (*neighbor) -> GetMaterialPosition();

	        	double neighborMaterialPositionX=neighborMaterialPosition.GetXCoord();
	        	double neighborMaterialPositionY=neighborMaterialPosition.GetYCoord();

			   	   Atom <2> *neigh;
			   	   neigh=*neighbor;

			       vector <Atom <2>*> update_neighbor_bond_neighbor_list;

			   	   vector < Atom <2>* > neighbor_neighbors=neigh->BondNeighbor();
			   	   int neigh_size=neighbor_neighbors.size();

		    	   double distance=bond.MaterialBondDistance(*atm, *neigh);

		    	   double initial_distance=bond.InitialBondDistance(*atm, *neigh);
		    	   double criticalConfigForce=(3.77*1.0/(distance*initial_distance))*0.999;

		    	   double configForce=0.0;

		    	   configForce=force.interatomicConfigForceValue(*atm,*neigh);

		    	   if ( abs(configForce) < criticalConfigForce )

		    	   {

		    		   update_neighbor_list.push_back(*neighbor);
	//	    		   update_neighbor_bond_neighbor_list.push_back(*atom);



		    	   }

		    	   else if (abs(configForce) >= criticalConfigForce && neighborMaterialPositionX >= crackRegion_MinX  && neighborMaterialPositionY >= crackRegion_MinY
			  				 && neighborMaterialPositionX<= crackRegion_MaxX  && neighborMaterialPositionY <= crackRegion_MaxY)
		    	   {

		    		   cout << "criticalConfigForce: " << criticalConfigForce << " configForce: "<< abs(configForce) << " neighbor id : "<< neigh->GetID()<< " atom id: "<< atm->GetID() << endl;

		    		   neighbor_neighbors.erase(find(neighbor_neighbors.begin(), neighbor_neighbors.end(),(*atom)));


		    	   }

		    	   (*neighbor)->SetBondNeighbors(neighbor_neighbors);

			       	if (neighbor_neighbors.size()!= neigh_size)
			       	{

			       	  	cout << "neighbor_neighbors size: " << neighbor_neighbors.size() << " neighbor size: " <<  neigh_size << endl;

			       	}

	    	}


	    	for (At bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	    	{

	        	Point <2> bond_neighborMaterialPosition= (*bond_neighbor) -> GetMaterialPosition();

	        	double bond_neighborMaterialPositionX=bond_neighborMaterialPosition.GetXCoord();
	        	double bond_neighborMaterialPositionY=bond_neighborMaterialPosition.GetYCoord();


			   	   Atom <2> *bond_neigh;
			   	   bond_neigh=*bond_neighbor;


			   	   vector <Atom <2>*> update_bond_neighbor_neighbor_list;

			   	   vector < Atom <2>* > bond_neighbor_neighbors=bond_neigh->Neighbor();

			   	   int bond_neigh_size=bond_neighbor_neighbors.size();

		    	   double distance=bond.MaterialBondDistance(*atm,*bond_neigh);
		    	   double initial_distance=bond.InitialBondDistance(*atm,*bond_neigh);

		    	   double criticalConfigForce=(3.77*1.0/(distance*initial_distance))*0.999;

		    	   double configForce=0.0;

		    	   configForce=force.interatomicConfigForceValue(*atm,*bond_neigh);

		    	   if ( abs(configForce) < criticalConfigForce )

		    	   {

		    		   update_bond_neighbor_list.push_back(*bond_neighbor);
	//	    		   update_bond_neighbor_neighbor_list.push_back(*atom);

		    	   }

		    	   else if ( abs(configForce) >= criticalConfigForce && bond_neighborMaterialPositionX >= crackRegion_MinX  && bond_neighborMaterialPositionY >= crackRegion_MinY
		  				 && bond_neighborMaterialPositionX<= crackRegion_MaxX  && bond_neighborMaterialPositionY <= crackRegion_MaxY )
		    	   {

		    		   cout << "criticalConfigForce: " << criticalConfigForce << " configForce: "<< abs(configForce) << " bond_neighbor id : "<< bond_neigh->GetID()<< " atom id: "<< atm->GetID() << endl;

		    		   bond_neighbor_neighbors.erase(find(bond_neighbor_neighbors.begin(), bond_neighbor_neighbors.end(),(*atom)));

		    	   }

		      	   (*bond_neighbor)->SetNeighbors(bond_neighbor_neighbors);

		       	if (bond_neighbor_neighbors.size()!= bond_neigh_size)
		       	{

		       	  	cout << "bond_neighbor_neighbors size: " << bond_neighbor_neighbors.size() << " bond_neighbor size: " <<  bond_neigh_size << endl;

		       	}

	    	}

	    	if (update_neighbor_list.size()!= neighbors.size()
	    		|| update_bond_neighbor_list.size()!=bond_neighbors.size())
	    	{

	    		(*atom)->SetNeighbors(update_neighbor_list);
	    	  	(*atom)->SetBondNeighbors(update_bond_neighbor_list);

	    	  	cout << "update_neighbor_list size: " << update_neighbor_list.size() << " neighbors size: " << neighbors.size() << " update_bond_neighbor_list size: " <<
	    	  			update_bond_neighbor_list.size() << " bond_neighbors size: " << bond_neighbors.size() << endl;

	    	}

			}

		}
	return (atoms);

}

