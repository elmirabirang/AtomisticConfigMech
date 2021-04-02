/*
 * modification_interaction_list_3d.cpp
 *
 *  Created on: Feb 29, 2020
 *      Author: S.Elmira Birang.O
 *      Central Institute of Scientific Computing
 *      University of Erlangen-Nuremberg
 *      This function modifies neighbor list of atoms met configurational-force-criterion for 3D lattice
 */




#include <vector>
#include "atom.h"
#include "point.h"
#include "force.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

using namespace std;

template <int dim>
void ModifyNeighborList(vector < Atom<dim>* > atoms, double alpha1, double alpha2,
                        double r_01, double r_02, double E1, double E2, double delta,
                        double r_cut, double h, double r_s1, double r_s2, double r_s3,
                        double s1, double s2, double s3,  double a, double beta1,
						double beta2, double r_03, double r_04, double F0, double F2,
						double q1, double q2, double q3, double q4, double Q1, double Q2, double crack_line)
{

	    typedef typename vector < Atom <3>* >::const_iterator At;
	    Force <dim> force;

   		for (At at=atoms.begin(); at!=atoms.end(); ++at)
    		{

			    vector < Atom <dim>* > neighbors=(*at)->Neighbor();
			    vector < Atom <dim>* > bond_neighbors=(*at)->BondNeighbor();

    			Point <3> config_force_pair = force.configForceEamPair(*at, alpha1, alpha2,
                                                                       r_01, r_02, E1, E2, delta,
                                                                       r_cut, h, r_s1, r_s2, r_s3,
                                                                       s1, s2, s3);

    			Point <3> config_force_embedding = force.configForceEamEmbedding(*at, a, beta1,
                                                                                 beta2, r_03, r_04, F0, F2,
                                                                                 q1, q2, q3, q4, Q1, Q2,
                                                                                 h, r_cut);

    			Point <3>  config_force=(config_force_pair+config_force_embedding)*-1;
    			double length_config_force=config_force.PointNorm();

    			if (length_config_force>=10.367)
    			{

    			  Point <dim> atom_material_position=(*at)->GetMaterialPosition();
    			  double atom_z_coord=atom_material_position.GetZCoord();

    			  bool sign=false;

    			  if ((atom_z_coord-crack_line) > 0)
    			  {
    				  sign=true;

    			  }

    			  for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
    			  {

    				  Point <dim> neighbor_material_psoition=(*neighbor)->GetMaterialPosition();

    				  double neighbor_z_coord=neighbor_material_psoition.GetZCoord();

    				  if (((neighbor_z_coord-crack_line) < 0 && sign==true)|| ((neighbor_z_coord-crack_line) > 0 && sign==false))
    				  {

    					  vector < Atom <dim>* > neighbor_bond_neighbor=(*neighbor)->BondNeighbor();

    					  neighbor_bond_neighbor.erase(find(neighbor_bond_neighbor.begin(), neighbor_bond_neighbor.end(),(*at)));
    					  (*neighbor)->SetBondNeighbors(neighbor_bond_neighbor);

    					  neighbors.erase(find(neighbors.begin(), neighbors.end(),(*neighbor)));

    				  }



    			  }

    			  (*at)->SetNeighbors(neighbors);


    			  for (At bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
    			  {

    				  Point <dim> bond_neighbor_material_psoition=(*bond_neighbor)->GetMaterialPosition();

    				  double bond_neighbor_z_coord=bond_neighbor_material_psoition.GetZCoord();

    				  if (((bond_neighbor_z_coord-crack_line) < 0 && sign==true) || ((bond_neighbor_z_coord-crack_line) > 0 && sign==false))
    				  {

    					  vector < Atom <dim>* > bond_neighbor_neighbor=(*bond_neighbor)->Neighbor();

    					  bond_neighbor_neighbor.erase(find(bond_neighbor_neighbor.begin(), bond_neighbor_neighbor.end(),(*at)));
    					  (*bond_neighbor)->SetNeighbors(bond_neighbor_neighbor);

    					  bond_neighbors.erase(find(bond_neighbors.begin(), bond_neighbors.end(),(*bond_neighbor)));

    				  }


    			  }

    			  (*at)->SetNeighbors(neighbors);
    			  (*at)->SetBondNeighbors(bond_neighbors);

    		  }

    		}

}




