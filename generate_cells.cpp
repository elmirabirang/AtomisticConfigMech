/*
 * generate_cells.cpp
 *
 *  Created on: Mar 15, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nuremberg
 */

#include "cell.h"
#include "point.h"
#include <vector>
#include <iostream>

using namespace std;


template <int dim>
vector <  Cell <2>* > GenerateCells(double r_cut, double lattice_constant, int number_x, int number_y, int number_z)
{

//r cut should be a multiple of lattice_constant

	vector <  Cell <2>* > cells;

	switch (dim)
	{

	case (2):
	{

		int number_cell_x=number_x*lattice_constant/r_cut;
		int number_cell_y=number_y*lattice_constant/r_cut;

		double side_length_x=0;
		int cell_id=0;
		int cell_index=0;


		for (int i=0; i<number_cell_x; ++i)
		{
			double side_length_y=0;

			for (int j=0; j<number_cell_y; ++j)
			{

			    Cell <dim> *cell=new Cell <dim>;
			    Point <dim> origin;
			    Point <dim> top_right;

			    cell_index=i+j;

			    origin.SetXCoord(side_length_x);
			    origin.SetYCoord(side_length_y);

			    top_right.SetXCoord(side_length_x+r_cut);
			    top_right.SetYCoord(side_length_y+r_cut);

			    cell->SetOrigin(origin);
			    cell->SetTopRight(top_right);
			    cell->SetCellID(cell_id);
			    cell->SetCellIndex(cell_index);
			    cells.push_back(cell);

			    side_length_y+=r_cut;
			    cell_id+=1;

			}

			side_length_x+=r_cut;
//			cell_id+=1;

		}

		break;
	}


	case (3):
	{

		int number_cell_x=number_x*lattice_constant/r_cut;
		int number_cell_y=number_y*lattice_constant/r_cut;
		int number_cell_z=number_z*lattice_constant/r_cut;

		double side_length_z=0;
		int cell_id=0;
		int cell_index=0;


		for (int i=0; i<number_cell_z; ++i )
		{
			double side_length_x=0;

		  for (int j=0; j<number_cell_x; ++j)
		  {
			  double side_length_y=0;

			for (int k=0; k<number_cell_y; ++k)
			{

			    Cell <dim> *cell=new Cell <dim>;
			    Point <dim> origin;
			    Point <dim> top_right;

			    cell_index=i+j+k;

			    origin.SetXCoord(side_length_x);
			    origin.SetYCoord(side_length_y);
			    origin.SetZCoord(side_length_z);

			    top_right.SetXCoord(side_length_x+r_cut);
			    top_right.SetYCoord(side_length_y+r_cut);
			    top_right.SetZCoord(side_length_z+r_cut);

			    cell->SetOrigin(origin);
			    cell->SetTopRight(top_right);
			    cell->SetCellID(cell_id);
			    cell->SetCellIndex(cell_index);
			    cells.push_back(cell);

			    side_length_y+=r_cut;
			    cell_id+=1;

			}

			side_length_x+=r_cut;
//			cell_id+=1;

		  }

		  side_length_z+=r_cut;
//		  cell_id+=1;

	  }

		break;

	}


	}


return (cells);


}


