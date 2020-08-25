/*
 * test_cell.cpp
 *
 *  Created on: Jan 30, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include "point.h"
#include <vector>
#include "cell.h"
#include "cell.cpp"

#include "subcell.h"
#include "subcell.cpp"
#include "atom.h"
#include <iostream>

int main()
{
	Point <2> corner_lb(0,0);
	Point <2> corner_rt(1,1);


	Cell <2> cell0;
	cell0.SetCellCorner_LB(corner_lb);
	cell0.SetCellCorner_RT(corner_rt);
	cell0.SetCellID(0);

	Point <2> lb_corner=cell0.GetCellCorner_LB();
	Point <2> rt_corner=cell0.GetCellCorner_RT();
	int id_0=cell0.GetCellID();




	Point <2> subcell_corner_lb(0,0);
	Point <2> subcell_corner_rt(0.5,0.5);

	SubCell <2> cell1;
	cell1.SetsubCellCorner_RT(subcell_corner_rt);
	cell1.SetsubCellCorner_LB(subcell_corner_lb);
	cell1.SetsubCellID(1);

	Point <2> lb_subcell_corner=cell1.GetsubCellCorner_LB();
	Point <2> rt_subcell_corner=cell1.GetsubCellCorner_RT();
	int id_1=cell1.GetsubCellID();

	cout << "id: " << id_0 << endl;
	cout << "cell_corner_rt_x: "<< rt_corner.GetXCoord()<< " cell_corner_rt_y: " << rt_corner.GetYCoord() << endl;
	cout << "cell_corner_lb_x: "<< lb_corner.GetXCoord()<< " cell_corner_lb_y: " << lb_corner.GetYCoord() << endl;

	cout << "id: " << id_1 << endl;
	cout << "subcell_corner_rt_x: "<< rt_subcell_corner.GetXCoord()<< " subcell_corner_rt_y: " << rt_subcell_corner.GetYCoord() << endl;
	cout << "subcell_corner_lb_x: "<< lb_subcell_corner.GetXCoord()<< " subcell_corner_lb_y: " << lb_subcell_corner.GetYCoord() << endl;

}




