/*
 * assign_region.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "boundary.h"
#include "point.h"

#include <vector>
#include <iostream>
#include<cmath>

using namespace std;

template <int dim>
float triangleArea( double xA, double yA, double xB, double yB, double xC, double yC) {

   return abs((xA*(yB-yC) + xB*(yC-yA)+xC*(yA-yB))/2.0);

}


template <int dim>
void AssignRegion(Atom <dim>* &atom, vector < Boundary <dim> > boundaries)
{
    Point <dim> material_position =atom->GetMaterialPosition();
//    int dimension=material_position.GetDim();
    // if atom does not belong to any boundary, its boundary id is zero,
    //otherwise it will be modified thorugh following loops.
    atom->SetAtomRegion(0);

	if (dim==2)
	{
	for (int i=0; i<boundaries.size(); ++i)
	{

		int bid=boundaries[i].GetBid();
		if (bid <10 && i==0)
		{

			double x_zone_min=boundaries[0].GetBoundaryX0();;
			double y_zone_min=boundaries[0].GetBoundaryY0();;

			double x_zone_max=boundaries[0].GetBoundaryX1();
			double y_zone_max=boundaries[0].GetBoundaryY1();

			Point <dim> material_position;
			material_position=atom->GetMaterialPosition();

			double xa=material_position.GetXCoord();
			double ya=material_position.GetYCoord();

			if ( xa> x_zone_min && xa < x_zone_max && ya>y_zone_min && ya<y_zone_max )
			{
				atom->setAtomZone(1);
			}
			else
			{
				atom->setAtomZone(0);
			}

		}

		else if (bid <10 && i!=0)
		{


			double x0=boundaries[i].GetBoundaryX0();
			double y0=boundaries[i].GetBoundaryY0();

			double x1=boundaries[i].GetBoundaryX1();
			double y1=boundaries[i].GetBoundaryY1();

			Point <dim> material_position;
			material_position=atom->GetMaterialPosition();

			double xa=material_position.GetXCoord();
			double ya=material_position.GetYCoord();

			if (xa>=x0 && xa<=x1 && ya>=y0 && ya<=y1)
			{
				atom->SetAtomRegion(bid);
//				cout << "id: " << bid << endl;
				break;
			}

		}

		else{

			double xA=boundaries[i].GetBoundaryX0();
			double yA=boundaries[i].GetBoundaryY0();


			double xB=boundaries[i].GetBoundaryX1();
			double yB=boundaries[i].GetBoundaryY1();


			double xC=boundaries[i].GetBoundaryX2();
			double yC=boundaries[i].GetBoundaryY2();


			Point <dim> material_position;
			material_position=atom->GetMaterialPosition();

			double xP=material_position.GetXCoord();
			double yP=material_position.GetYCoord();


			float area = triangleArea <dim> (xA, yA, xB, yB, xC, yC);          //area of triangle ABC

			float area1 = triangleArea <dim> (xP, yP, xB, yB, xC, yC);         //area of PBC
			float area2 = triangleArea <dim> (xA, yA, xP, yP, xC, yC);         //area of APC
			float area3 = triangleArea <dim> (xA, yA, xB, yB, xP, yP);        //area of ABP

			if (roundf(area)==roundf(area1+area2+area3))
			{
				atom->SetAtomRegion(bid);
				break;
			}




		}


	}
	}

	else{

		for (int i=0; i<boundaries.size(); ++i)
		{
			int bid=boundaries[i].GetBid();

			if (bid<10)
			{

			double x0=boundaries[i].GetBoundaryX0();
			double y0=boundaries[i].GetBoundaryY0();
			double z0=boundaries[i].GetBoundaryZ0();

			double x1=boundaries[i].GetBoundaryX1();
			double y1=boundaries[i].GetBoundaryY1();
			double z1=boundaries[i].GetBoundaryZ1();

			Point <dim> material_position;
			material_position=atom->GetMaterialPosition();

			double xa=material_position.GetXCoord();
			double ya=material_position.GetYCoord();
			double za=material_position.GetZCoord();

			if (xa>=x0 && xa<=x1 && ya>=y0 && ya<=y1 && za>=z0 && za<=z1)
			{
				atom->SetAtomRegion(bid);
				break;

			}
			}

//
//			else{
//
//
//
//
//			}


		}
	}


}
