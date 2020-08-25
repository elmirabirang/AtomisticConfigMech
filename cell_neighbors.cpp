/*
 * cell_neighbors.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */


#include "cell.h"
#include "math.h"
#include <vector>
#include <iostream>

using namespace std;

template <int dim>
void CellsNeighbors(vector < Cell<dim>* > &cells, double r_cut)
{
	typedef typename vector < Cell <dim>* >::iterator Cel;

	switch (dim)
	{

	case (2):
	{

		for (Cel celli=cells.begin(); celli!=cells.end(); ++celli)
		{
			vector < Cell <dim>* > celli_neighbors;
			int celli_id=(*celli)->GetCellID();

			Point <dim> celli_origin=(*celli)->GetOrigin();

			double celli_x=celli_origin.GetXCoord();
			double celli_y=celli_origin.GetYCoord();

			for (Cel cellj=cells.begin(); cellj!=cells.end(); ++cellj)
			{

				    int cellj_id=(*cellj)->GetCellID();

					Point <dim> cellj_origin=(*cellj)->GetOrigin();

					double cellj_x=cellj_origin.GetXCoord();
					double cellj_y=cellj_origin.GetYCoord();

					double distance = pow ( (cellj_x-celli_x),2 )+ pow ( (cellj_y-celli_y) , 2);

					if (distance <= 2*pow(r_cut,2) )
					{

						celli_neighbors.push_back(*cellj);

					}



				}

			(*celli) -> SetCellNeighbors(celli_neighbors);


			}


			break;

		}







	case (3):{


		for (Cel celli=cells.begin(); celli!=cells.end(); ++celli)
		{
			vector < Cell <dim>* > celli_neighbors;
			int celli_id=(*celli)->GetCellID();

			Point <dim> celli_origin=(*celli)->GetOrigin();

			double celli_x=celli_origin.GetXCoord();
			double celli_y=celli_origin.GetYCoord();
			double celli_z=celli_origin.GetZCoord();

			for (Cel cellj=cells.begin(); cellj!=cells.end(); ++cellj)
			{

				    int cellj_id=(*cellj)->GetCellID();

					Point <dim> cellj_origin=(*cellj)->GetOrigin();

					double cellj_x=cellj_origin.GetXCoord();
					double cellj_y=cellj_origin.GetYCoord();
					double cellj_z=cellj_origin.GetZCoord();

					double distance = pow((cellj_x-celli_x),2 )+ pow((cellj_y-celli_y),2) +pow((cellj_z-celli_z),2);

					if (distance <= 3*pow(r_cut,2) )
					{
						celli_neighbors.push_back(*cellj);

					}


				}

			(*celli) -> SetCellNeighbors(celli_neighbors);

			}

	        break;

		}


	}
}






