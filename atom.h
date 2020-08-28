/*
 * atom.h
 *
 *  Created on: Nov 16, 2018
 *      Author: S.Elmira Birang.O
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#pragma once

#include "point.h"
#include <vector>

using namespace std;

template <int dim>
class Atom {

 public:

   Atom();
   ~Atom();
   Atom(Point <dim>, Point <dim>, int ID);

   void SetMaterialPosition(Point <dim> reference_position);
   Point <dim> GetMaterialPosition();

   void SetSpatialPosition(Point <dim> current_position);
   Point <dim> GetSpatialPosition();

   void SetInitialPosition(Point <dim> unrelaxed_position);
   Point <dim> GetInitialPosition();

   void SetID(int id);
   int GetID();

   //This list of neighbors to calculate total potential energy of the system.
   vector < Atom <dim>* > Neighbor();
   //The list of neighbors to calculate ResultantBondVec and ResultantBondVec_Spatial
   vector < Atom <dim>* > BondNeighbor();
   
   void SetNeighbors(vector < Atom <dim>* > neighbors);
   void SetBondNeighbors(vector < Atom <dim>* > &bond_neighbors );

   void SetAtomRegion(int bID);
   int GetAtomRegion();

   void SetForce(Point <dim> pair_force);
   Point <dim> GetForce();

   void  SetCellID(int cell_id);
   int GetCellID();

   Point <dim> GetConfigForce();
   void SetConfigForce(Point <dim> atom_config_force);

   Point <dim> GetDeformForce();
   void SetDeformForce(Point <dim> atom_deform_force);

   void setAtomZone(int zone);
   int getAtomZone();

   vector <double> getAtomicStress();
   void setAtomicStress(vector <double> atomicStress);


 private:

   Point <dim> material_position;
   Point <dim> spatial_position;
   Point <dim> initial_position;

   int ID;
   int Cell_ID;
   double cutRadius;
   int boundary_id;

   vector < Atom <dim>* > atom_neighbor;
   vector < Atom <dim>* > bond_neighbor;

   Point <dim> force;
   Point <dim> config_force;
   Point <dim> deform_force;

   int atomZone;

   vector <double> AtomicStress;

};

