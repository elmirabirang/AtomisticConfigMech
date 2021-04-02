/*
 * test_assign_atom_to_cell.cpp
 *
 *  Created on: Feb 1, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */


#include "cell.h"
#include "subcell.h"
#include "point.h"
#include "atom.h"
#include <vector>
#include <iostream>


#include "assign_atom_to_cell.cpp"
#include "unrelaxed_config.cpp"

using namespace std;

int main()
{

	vector < Atom <2> > atoms;

	atoms=UnrelaxedConfigGenerator <2> (5, 9, 0, 1.1,2.5);

	vector < Point <2> > positions;

	for (int i=0; i< atoms.size(); ++i)
	{

		Point <2> atom_position=atoms[i].GetMaterialPosition();

		positions.push_back(atom_position);

	}

	vector < Cell <2> > cells;
	vector < SubCell <2> > subcells;


	Cell <2> cell0;
	cell0.SetCellID(0);

	Point <2> cell_rt(4.41,7.63);
	Point <2> cell_lb(0,0);
	cell0.SetCellCorner_LB(cell_lb);
	cell0.SetCellCorner_RT(cell_rt);

	SubCell <2> subcell0;
	subcell0.SetsubCellID(0);

	Point <2> subcell0_rt(2.22,3.83);
	Point <2> subcell0_lb(0,0);
	subcell0.SetsubCellCorner_LB(subcell0_lb);
	subcell0.SetsubCellCorner_RT(subcell0_rt);

	SubCell <2> subcell1;
	subcell1.SetsubCellID(1);

	Point <2> subcell1_rt(4.42,3.83);
	Point <2> subcell1_lb(2.23,0);
	subcell1.SetsubCellCorner_LB(subcell1_lb);
	subcell1.SetsubCellCorner_RT(subcell1_rt);

	SubCell <2> subcell2;
	subcell2.SetsubCellID(2);

	Point <2> subcell2_lb(0,3.84);
	Point <2> subcell2_rt(2.22,7.63);
	subcell2.SetsubCellCorner_RT(subcell2_rt);
	subcell2.SetsubCellCorner_LB(subcell2_lb);

	SubCell <2> subcell3;
	subcell3.SetsubCellID(3);

	Point <2> subcell3_lb(2.25,3.84);
	Point <2> subcell3_rt(4.42,7.63);
	subcell3.SetsubCellCorner_RT(subcell3_rt);
	subcell3.SetsubCellCorner_LB(subcell3_lb);

	subcells.push_back(subcell0);
	subcells.push_back(subcell1);
	subcells.push_back(subcell2);
	subcells.push_back(subcell3);

	cell0.SetCellsubCells(subcells);

	cells.push_back(cell0);

	vector < Cell <2> > Model_cells=AssignAtomToCell <2> (atoms,positions, cells);

	typedef typename vector < Cell <2> >::iterator cell;
	vector < Point <2>> points_in_cell;
	points_in_cell=Model_cells[0].GetCellPoints();

	typedef typename vector <Point <2> >::iterator point;
	cout << "size: " <<points_in_cell.size() << endl;

	for(point P=points_in_cell.begin(); P!=points_in_cell.end(); ++P)
	{
		double xp=P->GetXCoord();
		double yp=P->GetYCoord();

		cout << "( " << xp << ", " << yp << " )" <<endl;
	}

	vector <SubCell <2>> cell_subcells=Model_cells[0].GetCellsubCells();

	typedef typename vector < SubCell <2> >::iterator subcell;
	typedef typename vector <Point <2> >::iterator point_subcell;

	for (subcell subC=cell_subcells.begin(); subC!= cell_subcells.end(); ++subC)
	{
		vector <Point <2>> points_in_subcell=subC->GetsubCellPoints();

		cout << "subcell id: " << subC->GetsubCellID()<<endl;

		for (point_subcell Psubcell=points_in_subcell.begin(); Psubcell!=points_in_subcell.end(); ++Psubcell)
		{
			double Px=Psubcell->GetXCoord();
			double Py=Psubcell->GetYCoord();

			cout << "( " << Px << ", " << Py << " )" << endl;

		}


	}




	vector < Atom <2> > atoms_incell;

	atoms_incell=Model_cells[0].GetCellAtoms();
	cout << "atoms_incell: " << atoms_incell.size() << endl;

	typedef typename vector < Atom <2> > :: iterator atom;

	for (atom A=atoms_incell.begin(); A!=atoms_incell.end(); ++A)
	{
		Point <2> material_position=A->GetMaterialPosition();

		double A_X=material_position.GetXCoord();
		double A_Y=material_position.GetYCoord();

		cout << "( " << A_X << ", " << A_Y <<" )" << endl;

	}

	typedef typename vector < SubCell <2> >::iterator sub_cells;
	typedef typename vector <Atom <2> >::iterator atom_subcell;

	for (sub_cells subC=cell_subcells.begin(); subC!= cell_subcells.end(); ++subC)
	{
		vector < Atom <2> > atoms_insubcell=subC->GetsubCellAtoms();

		cout << "subcell id: " << subC->GetsubCellID()<<endl;

		for (atom_subcell Asubcell=atoms_insubcell.begin(); Asubcell!=atoms_insubcell.end(); ++Asubcell)
		{
			Point <2> material_position_subcell=Asubcell->GetMaterialPosition();
			double Px=material_position_subcell.GetXCoord();
			double Py=material_position_subcell.GetYCoord();

			cout << "( " << Px << ", " << Py << " )" << endl;

		}


	}





return 0;
}



