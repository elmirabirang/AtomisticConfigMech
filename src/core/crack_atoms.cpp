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

template<int dim> CrackAtom <dim>::CrackAtom() {}
template<int dim> CrackAtom <dim>::~CrackAtom() {}
template<int dim> vector<int> CrackAtom<dim>::getCrackBottomAtoms() { return crack_bottom_atoms; }
template<int dim> void CrackAtom <dim>::setCrackBottomAtoms(vector<int> bottom_atoms) { this->crack_bottom_atoms = bottom_atoms; }
template<int dim> vector<int> CrackAtom<dim>::getCrackTopAtoms() { return crack_top_atoms; }

//template <int dim>
//void CrackAtom <dim>::setCrackTopAtoms(vector <Atom <dim>* > top_atoms)
//{
//
//	this -> crack_top_atoms = top_atoms;
//
//}

//crackAtom <2>* crack_atom(Atom <2>* atom)
template<int dim> void CrackAtom<dim>::setCrackRegion(int i) {
//	brittle_crack_atom=this;
	Point<dim> material_position_crack_atom = this->atoms->getMaterialPosition(i);
	double crack_atom_y = material_position_crack_atom.GetYCoord();
	vector<int> bottom_atoms;
	vector<int> top_atoms;

	cout << "atom y: " << crack_atom_y<< endl;
	cout << "atom x: " << material_position_crack_atom.GetXCoord() << endl;
	//cout << "size: " << crack_atom_neighbors.size() <<endl;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		Point<dim> material_position_neighbor = this->atoms->getMaterialPosition(j);
		double neighbor_y = material_position_neighbor.GetYCoord();

		if(neighbor_y <= crack_atom_y) {
			bottom_atoms.push_back(j);
			cout << "here00" <<endl;
		} else {
			top_atoms.push_back(j);
			cout << "here0" <<endl;
		}
	)

    //TODO: set the atom at bottom region and top region. change energy crack atom
}

template<int dim> double CrackAtom<dim>::energyCrackAtom(CrackAtom<dim> *crack_atom, double Rcut, double sigma, double epsilon) {
	Energy<dim> energy(sigma, epsilon);
	this->brittle_crack_atom = crack_atom;
	double cut_radius = Rcut;
	double debond_energy = 0.0;

    LOOP_OVER_ATOMS(this->bottom_atoms, i,
		Point<dim> bottom_atom_material_position = this->bottom_atoms->getMaterialPosition(i);
		double bottom_atom_x = bottom_atom_material_position.GetXCoord();
		double bottom_atom_y = bottom_atom_material_position.GetYCoord();

        LOOP_OVER_ATOMS(this->top_atoms, j,
			Point<dim> top_atom_material_position = this->top_atoms->getMaterialPosition(j);
			double top_atom_x = top_atom_material_position.GetXCoord();
			double top_atom_y = top_atom_material_position.GetYCoord();
			double distance = pow((top_atom_x - bottom_atom_x), 2) + pow((top_atom_y - bottom_atom_y), 2);
			if(distance < (cut_radius * cut_radius)) {
                // FIXME: Get original positions for atoms on global set or just pass dist/rsq to this function
				double interatomic_energy = energy.InteratomicEnergy(i, j);
				debond_energy += interatomic_energy;
			}
		)
	)

	return debond_energy;
}
