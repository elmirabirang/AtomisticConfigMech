/*
 * atom.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: S.Elmira Birang.O
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#include "atom.h"
#include <vector>
using namespace std;


template <int dim>
Atom<dim>::Atom()
{}

template <int dim>
Atom<dim>::~Atom()
{}

template <int dim>
Atom<dim>::Atom(Point<dim> initial_position, Point <dim> current_position, int id)
{
  this->material_position=initial_position;
  this->spatial_position=current_position;
  this->ID=id;
}

template <int dim>
void Atom<dim>::SetMaterialPosition(Point <dim> reference_position)
{
  this->material_position=reference_position;
}

template <int dim>
Point<dim> Atom <dim>::GetMaterialPosition()
{
  return(material_position);
}

template <int dim>
void Atom<dim>::SetSpatialPosition(Point<dim> current_position)
{
  this -> spatial_position=current_position;
}

template <int dim>
Point <dim> Atom <dim>::GetSpatialPosition()
{
  return(spatial_position);
}

template <int dim>
void Atom <dim>::SetInitialPosition(Point <dim> unrelaxed_position)
{
	this -> initial_position = unrelaxed_position;
}

template <int dim>
Point <dim> Atom <dim>::GetInitialPosition()
{
	return (initial_position);
}

template <int dim>
void Atom <dim>::SetID(int id)
{
  this -> ID=id;
}

template <int dim>
int Atom <dim>:: GetID()
{
  return(ID);
}

template <int dim>
vector < Atom <dim>* > Atom<dim>::Neighbor()
{
	return(atom_neighbor);

}

template <int dim>
void Atom <dim>:: SetNeighbors(vector < Atom <dim>* > neighbors)
{
	this -> atom_neighbor=neighbors;

}

template <int dim> 
vector < Atom <dim>* > Atom <dim> :: BondNeighbor()
{
 return(bond_neighbor);
}

template <int dim> 
void Atom <dim> :: SetBondNeighbors( vector < Atom <dim>* > &bond_neighbors)
{
   this-> bond_neighbor=bond_neighbors;	
}


template <int dim>
void Atom <dim> :: SetAtomRegion(int bID)
{

      this -> boundary_id=bID;

}


template <int dim>
int Atom <dim> ::GetAtomRegion()
{

     return(boundary_id);

}

template <int dim>
void Atom <dim>:: SetForce(Point <dim> initial_force)
{
	this -> force = initial_force;

}

template <int dim>
Point <dim> Atom <dim>:: GetForce()
{
	return(force);

}

template <int dim>
int Atom <dim>::GetCellID()
{
	return (Cell_ID);
}

template <int dim>
void Atom <dim>::SetCellID(int cell_id)
{
	this -> Cell_ID= cell_id;
}

template <int dim>
void Atom <dim>::SetConfigForce(Point <dim> atom_config_force)
{
   this->config_force=atom_config_force;
}

template <int dim>
Point <dim> Atom <dim>::GetConfigForce()
{
   return(config_force);
}

template <int dim>
void Atom <dim>::SetDeformForce(Point <dim> atom_deform_force)
{
   this->deform_force=atom_deform_force;
}

template <int dim>
Point <dim> Atom <dim>::GetDeformForce()
{
   return(deform_force);
}

template <int dim>
void Atom <dim>::setAtomZone(int zone)
{
	this -> atomZone=zone;
}

template <int dim>
int Atom<dim>::getAtomZone()
{
	return(atomZone);
}

template <int dim>
void Atom<dim>::setAtomicStress(vector <double> atomicStress)
{
	this -> AtomicStress=atomicStress;
}

template <int dim>
vector <double> Atom<dim>::getAtomicStress()
{
	return(AtomicStress);
}








