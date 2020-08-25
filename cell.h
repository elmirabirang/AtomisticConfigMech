/*
 * cell.h
 *
 *  Created on: Jan 30, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-nuremberg
 */

#pragma once

#include <vector>
#include "point.h"
#include "subcell.h"
#include "atom.h"

using namespace std;

template <int dim>
class Cell{

public:

	Cell();
	~Cell();

	int GetCellID();
	void SetCellID(int id);

	int GetCellIndex();
	void SetCellIndex(int index);

	vector < Cell <dim>* > GetCellNeighbors();
	void SetCellNeighbors(vector < Cell <dim>* > neighbors);

	vector < Atom <dim>* > GetCellAtoms();
	void AtomsInCell (vector < Atom <dim>* > cell_atoms);

	void SetOrigin(Point <dim> Celorigin );
	Point <dim> GetOrigin();

	void SetTopRight(Point <dim> CelTR );
	Point <dim> GetTopRight();




private:

	Point <dim> cell_origin;
	Point <dim> cell_TR;
	int cell_id;
	int cell_index;
	vector < Atom <dim>* > atoms_in_cell ;
	vector < Cell <dim>* > cell_neighbors;

};
