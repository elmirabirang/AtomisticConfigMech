/*
 * subcell.cpp
 *
 *  Created on: Jan 31, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */

#include "subcell.h"

#include "point.h"
#include "point.cpp"

#include "atom.h"
#include "atom.cpp"

using namespace std;

template <int dim>
SubCell<dim>::SubCell()
{}

template <int dim>
SubCell<dim>::~SubCell()
{}

template <int dim>
Point <dim> SubCell<dim>::GetsubCellCorner_RT()
{
	return (subcell_corner_rt);
}


template <int dim>
void SubCell<dim>::SetsubCellCorner_RT(Point <dim> subc_corner_rt)
{
   this-> subcell_corner_rt=subc_corner_rt;
}


template <int dim>
Point <dim> SubCell<dim>::GetsubCellCorner_LB()
{
	return (subcell_corner_lb);
}


template <int dim>
void SubCell<dim>::SetsubCellCorner_LB(Point <dim> subc_corner_lb)
{
   this-> subcell_corner_lb=subc_corner_lb;
}


template <int dim>
Point <dim> SubCell<dim>::GetsubCellCorner_Inner()
{
	return (subcell_corner_inner);
}


template <int dim>
void SubCell<dim>::SetsubCellCorner_Inner(Point <dim> subc_corner_inner)
{
   this-> subcell_corner_inner=subc_corner_inner;
}



template <int dim>
int SubCell<dim>::GetsubCellID()
{
	return(subcell_id);
}


template <int dim>
void SubCell<dim>::SetsubCellID(int subc_id)
{
	this-> subcell_id=subc_id;
}


template <int dim>
vector <int> SubCell<dim>::GetsubCellNeighbors()
{
	return (subcell_neighbors);

}


template <int dim>
void SubCell<dim>::SetsubCellNeighbors(vector <int> subc_neighbors)
{
	this -> subcell_neighbors=subc_neighbors;

}

template <int dim>
vector < Point <dim> > SubCell<dim>::GetsubCellPoints()
{
	return (points_in_subcell);
}


template <int dim>
void SubCell<dim>::PointsInsubCell(vector < Point<dim> > points_of_subcell)
{
	this ->points_in_subcell=points_of_subcell;
}

template <int dim>
vector < Atom <dim> > SubCell<dim>:: GetsubCellAtoms()
{
	return (atoms_in_subcell);
}

template <int dim>
void SubCell <dim>:: AtomsInSubCell(vector < Atom <dim> > subcell_atoms)
{
	this -> atoms_in_subcell=subcell_atoms;
}









