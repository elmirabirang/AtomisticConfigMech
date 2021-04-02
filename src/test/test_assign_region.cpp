/*
 * test_assign_region.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nürnberg
 */

#include "assign_region.cpp"

#include "boundary.h"
#include "atom.h"
#include "point.h"

#include "boundary.cpp"
#include "atom.cpp"
#include "point.cpp"

#include <vector>
#include <iostream>


using namespace std;
int main (){


   Boundary <2> boundary2d;

   boundary2d.SetTriangBoundaryRegion(1, 5, 5, 1, 3, 7,11);

   vector < Boundary <2> > boundaries;
   boundaries.push_back(boundary2d);

   Atom <2> atm;
   Point <2> mposition;

   mposition.SetXCoord(3);
   mposition.SetYCoord(8);

   atm.SetMaterialPosition(mposition);

   Atom <2> *neigh;
   neigh=&atm;



   AssignRegion(neigh,boundaries);

   int id=neigh->GetAtomRegion();

   cout << "id: " << id << endl;

   return 0;

}














































