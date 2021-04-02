/*
 * subcell.h
 *
 *  Created on: Jan 31, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#pragma once

#include "point.h"
#include "atom.h"
#include <vector>

using namespace std;

template <int dim>
class SubCell{

public :

	SubCell();
	~SubCell();

	Point <dim> GetsubCellCorner_RT();
	void SetsubCellCorner_RT(Point <dim> subc_corner_rt);

	Point <dim> GetsubCellCorner_LB();
	void SetsubCellCorner_LB(Point <dim> subc_corner_lb);

	Point <dim> GetsubCellCorner_Inner();
	void SetsubCellCorner_Inner(Point <dim> subc_corner_inner);

	int GetsubCellID();
	void SetsubCellID(int subc_id);

	vector<int> GetsubCellNeighbors();
	void SetsubCellNeighbors(vector <int> subc_neighbors);

	vector < Point <dim> > GetsubCellPoints();
	void PointsInsubCell(vector < Point<dim> > points_of_subcell);

	vector < Atom <dim> > GetsubCellAtoms();
	void AtomsInSubCell(vector < Atom <dim> > subcell_atoms);


private:

	Point <dim> subcell_corner_rt;
	Point <dim> subcell_corner_lb;
	Point <dim> subcell_corner_inner;

	int subcell_id;

	vector <int> subcell_neighbors;

	vector < Point <dim> > points_in_subcell;

	vector < Atom <dim> > atoms_in_subcell;

};


