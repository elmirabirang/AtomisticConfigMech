/*
 * crack_atoms.cpp
 *
 *  Created on: Jun 11, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */
// This class presents a type that describes the atoms which define the crack line in brittle crack propagation.

#include <cmath>
#include <iostream>
#include <vector>
//--
#include "atom.h"
#include "crack_atoms.h"
#include "energy.h"
#include "point.h"

using namespace std;

template <int dim> CrackAtom <dim>::CrackAtom() {}
template <int dim> CrackAtom <dim>::~CrackAtom() {}
template <int dim> vector < Atom <dim>* > CrackAtom <dim>::getCrackBottomAtoms() {
	return (crack_bottom_atoms);
}

template <int dim> void CrackAtom <dim>::setCrackBottomAtoms(vector < Atom <dim>* > bottom_atoms) {
	this -> crack_bottom_atoms = bottom_atoms;
}

template <int dim> vector <Atom <dim>* > CrackAtom <dim>::getCrackTopAtoms() {
   return(crack_top_atoms);
}

//template <int dim>
//void CrackAtom <dim>::setCrackTopAtoms(vector <Atom <dim>* > top_atoms)
//{
//
//	this -> crack_top_atoms = top_atoms;
//
//}

//crackAtom <2>* crack_atom(Atom <2>* atom)
template <int dim>
void CrackAtom <dim>::setCrackRegion(Atom <dim>* atom)
{

//	brittle_crack_atom=this;


	vector < Atom <dim>* > crack_atom_neighbors;
	vector < Atom <dim>* > crack_atom_bond_neighbors;
	crack_atom_neighbors=atom->Neighbor();
	crack_atom_bond_neighbors=atom->BondNeighbor();

	vector < Atom <dim>* > bottom_atoms;
	vector < Atom <dim>* > top_atoms;

	Point <dim> material_position_crack_atom = atom->GetMaterialPosition();

	double crack_atom_y=material_position_crack_atom.GetYCoord();
	cout << "atom y: " << crack_atom_y<< endl;
	cout << "atom x: " << material_position_crack_atom.GetXCoord() << endl;
	cout << "size: " << crack_atom_neighbors.size() <<endl;

	typedef typename vector < Atom <dim>* >::iterator Neigh;

	for (Neigh neighbor=crack_atom_neighbors.begin() ; neighbor!=crack_atom_neighbors.end() ; ++neighbor)
	{

		Point <dim> material_position_neighbor=(*neighbor)-> GetMaterialPosition();
		double neighbor_y= material_position_neighbor.GetYCoord();


		if (neighbor_y <= crack_atom_y)
		{

			bottom_atoms.push_back(*neighbor);
			cout << "here00" <<endl;

		}

		else
		{

			top_atoms.push_back(*neighbor);
			cout << "here0" <<endl;


		}

	}

	for (Neigh neighbor=crack_atom_bond_neighbors.begin() ; neighbor!=crack_atom_bond_neighbors.end() ; ++neighbor)
	{

		Point <dim> material_position_neighbor=(*neighbor)-> GetMaterialPosition();
		double neighbor_y= material_position_neighbor.GetYCoord();


		if (neighbor_y <= crack_atom_y)
		{

			bottom_atoms.push_back(*neighbor);
			cout << "here11" <<endl;

		}

		else
		{

			top_atoms.push_back(*neighbor);
			cout << "here1" <<endl;


		}

	}
//TODO: set the atom at bottom region and top region. change energy crack atom


}

template <int dim>
double CrackAtom <dim>::energyCrackAtom(CrackAtom <dim>* crack_atom, double Rcut , double sigma, double epsilon)
{

	Energy <dim> energy (sigma, epsilon);
	double debond_energy=0;

	this -> brittle_crack_atom=crack_atom;
	double cut_radius=Rcut;

	vector < Atom <dim>* > top_atoms=crack_atom->getCrackTopAtoms();
	vector < Atom <dim>* > bottom_atoms=crack_atom->getCrackBottomAtoms();

	typedef typename vector < Atom <dim>* >::iterator At;

	for (At bottom_atom=bottom_atoms.begin(); bottom_atom!= bottom_atoms.end(); ++bottom_atom)
	{

		Point <dim> bottom_atom_material_position=bottom_atom.GetMaterialPosition();

		double bottom_atom_x=bottom_atom_material_position.GetXCoord();
		double bottom_atom_y=bottom_atom_material_position.GetYCoord();

		for (At top_atom=top_atoms.begin(); top_atom!=top_atoms.end() ; ++top_atom)
		{

			Point <dim> top_atom_material_position=top_atom.GetMaterialPosition();

			double top_atom_x=top_atom_material_position.GetXCoord();
			double top_atom_y=top_atom_material_position.GetYCoord();

			double distance = pow((top_atom_x-bottom_atom_x),2)+pow((top_atom_y-bottom_atom_y),2);

			if (distance < (cut_radius*cut_radius) )
			{

				double interatomic_energy = energy.InteratomicEnergy (top_atom, bottom_atom);
				debond_energy+=interatomic_energy;

			}

		}

	}

	return (debond_energy);

}














































