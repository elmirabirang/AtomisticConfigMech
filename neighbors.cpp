/*
 * neighbors.cpp
 * get the list of atoms from unrelax_config.cpp
 * construct the neighbors of atom
 *
 *  Created on: Dec 3, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical enigneering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "bond.h"
#include "point.h"
#include "cell.h"
#include "subcell.h"
#include "atom.h"
#include <math.h>

#include <iostream>
#include <vector>

using namespace std;


template <int dim>
vector < Atom <dim>* > FindNeighbors(vector < Atom <dim>* > &unrelax_atoms, double cut_radius)
{

	Bond <dim> bond;

	typedef typename vector < Atom <dim>* >:: iterator Atomi;
	typedef typename vector < Atom <dim>* >:: iterator Atomj;

	//list of neighbors to calculate pair-wise energy.
	//list of neighbors to calculate ResultantbOndVec_Spatial & ResultantBondVec_Material
if (dim==2){

	for (Atomi atomi=unrelax_atoms.begin(); atomi!=unrelax_atoms.end(); ++atomi)
	{
		vector < Atom <dim>* > neighbors;
		int atomi_id=(*atomi)->GetID();

		int region=(*atomi)->GetAtomRegion();

		Point <dim> materialP=(*atomi)->GetMaterialPosition();

		double x_i=materialP.GetXCoord();
		double y_i=materialP.GetYCoord();

		switch (region){

		case 4:{


			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				int atomj_id=(*atomj)->GetID();
				int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=5 && region_j!=3)
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();
			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;

				     double distance = sqrt((x_coord*x_coord)+(y_coord*y_coord));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;
		}

		case 5:{
			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				int atomj_id=(*atomj)->GetID();
				int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=4 && region_j!=3 )
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();
			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;

				     double distance = sqrt(pow(x_coord,2)+pow(y_coord,2));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;

		}



		case 0:
		case 1:
		case 2:{

			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				 int atomj_id=(*atomj)->GetID();
				 int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=3)
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();
			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;

				     double distance = sqrt((x_coord*x_coord)+(y_coord*y_coord));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;

		}

		}

		(*atomi)->SetNeighbors(neighbors);

	}

	return(unrelax_atoms);
}



else{

	for (Atomi atomi=unrelax_atoms.begin(); atomi!=unrelax_atoms.end(); ++atomi)
	{
		vector < Atom <dim>* > neighbors;
		int atomi_id=(*atomi)->GetID();

		int region=(*atomi)->GetAtomRegion();

		Point <dim> materialP=(*atomi)->GetMaterialPosition();

		double x_i=materialP.GetXCoord();
		double y_i=materialP.GetYCoord();
		double z_i=materialP.GetZCoord();

		switch (region){

		case 4:{


			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				int atomj_id=(*atomj)->GetID();
				int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=5 && region_j!=3 )
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();

			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();
			    	 double z_j=neighborP.GetZCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;
				     double z_coord=z_j-z_i;

				     double distance = sqrt(pow(x_coord,2)+pow(y_coord,2)+pow(z_coord,2));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;
		}

		case 5:{

			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				int atomj_id=(*atomj)->GetID();
				int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=4 && region_j!=3 )
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();
			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();
			    	 double z_j=neighborP.GetZCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;
				     double z_coord=z_j-z_i;

				     double distance = sqrt((x_coord*x_coord)+(y_coord*y_coord)+(z_coord*z_coord));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;

		}

		case 0:
		case 1:
		case 2:
		{

			for (Atomj atomj=unrelax_atoms.begin(); atomj!=unrelax_atoms.end(); ++atomj)

			{
				 int atomj_id=(*atomj)->GetID();
				 int region_j=(*atomj)->GetAtomRegion();

			     if (atomi_id <  atomj_id && region_j!=3)
			     {

			    	 Point <dim> neighborP=(*atomj)->GetMaterialPosition();
			    	 double x_j=neighborP.GetXCoord();
			    	 double y_j=neighborP.GetYCoord();
			    	 double z_j=neighborP.GetZCoord();

				     double x_coord=x_j-x_i;
				     double y_coord=y_j-y_i;
				     double z_coord=z_j-z_i;

				     double distance = sqrt((x_coord*x_coord)+(y_coord*y_coord)+(z_coord*z_coord));

			    	 if (distance < cut_radius)

			    	 {
			    		 vector < Atom <dim>* > neighbors_of_atomj=(*atomj)->BondNeighbor();

			    		 neighbors_of_atomj.push_back(*atomi);
			    		 (*atomj)->SetBondNeighbors(neighbors_of_atomj);

			    		 neighbors.push_back(*atomj);

			    	 }
			     }
			}
			break;

		}

		}

		(*atomi)->SetNeighbors(neighbors);

	}

	return(unrelax_atoms);
}

}
