/*
 * Cell.cpp
 *
 *  Created on: Jan 30, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <vector>
//---
#include "atom.h"
#include "cell.h"
#include "point.h"

using namespace std;

template<int dim> Cell<dim>::Cell() {}
template<int dim> Cell<dim>::~Cell() {}
template<int dim> int Cell <dim>:: GetCellID() { return (cell_id); }
template<int dim> void Cell<dim>::SetCellID(int id) { this-> cell_id=id; }
template<int dim> vector < Cell <dim>* > Cell<dim>::GetCellNeighbors() { return (cell_neighbors); }
template<int dim> void Cell<dim>::SetCellNeighbors(vector < Cell <dim>* > neighbors) { this->cell_neighbors=neighbors; }
template<int dim> vector <int> Cell<dim>::GetCellAtoms() { return (atoms_in_cell); }
template<int dim> void Cell <dim>:: AtomsInCell(vector<int> cell_atoms) { this -> atoms_in_cell=cell_atoms; }
template<int dim> void Cell <dim>:: SetOrigin(Point <dim> Celorigin) { this -> cell_origin=Celorigin; }
template<int dim> Point <dim> Cell <dim>:: GetOrigin() { return (cell_origin); }
template<int dim> void Cell <dim>:: SetTopRight(Point <dim> CelTR) { this -> cell_TR=CelTR; }
template<int dim> Point <dim> Cell <dim>:: GetTopRight() { return (cell_TR); }
template<int dim> void Cell<dim>::SetCellIndex(int index) { this->cell_index=index; }
template<int dim> int Cell<dim>::GetCellIndex() { return (cell_index); }
