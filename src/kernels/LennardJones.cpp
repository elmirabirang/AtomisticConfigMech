/*
 * LennardJones.cpp
 *
 *  Created on: Oct 20, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "bond.h"
#include "point.h"
#include "Lennard_Jones.h"
#include "math.h"
#include "matrix.h"
#include <vector>
#include <string>

template <int dim>
LennardJones <dim>::LennardJones()
{}

template <int dim>
LennardJones <dim>::~LennardJones()
{}

template <int dim>
double LennardJones <dim>::LennardJones_Classic(Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta, double epsilon, double sigma)
{

	Bond <dim> bond;

	this -> atomi = atom_Alpha;
	this -> atomj = atom_Beta;

	double spatial_distance=bond.SpatialBondDistance(atomi, atomj);

	double term1=sigma/spatial_distance;
	double term2=pow(term1, 6);
	double term3=term2*term2;

	double energy= 4*epsilon*(term3-term2);

	Point <dim> spatial_bond_vec(0.,0.,0.);

	spatial_bond_vec=bond.SpatialBondNormal(atomi, atomj);

	Point <dim> interatomic_force=spatial_bond_vec*(4*epsilon/spatial_distance)*(-12*term3+6*term2);

	atomi.SetForce(atomi.GetForce()+interatomic_force);
	atomj.SetForce(atomj.GetForce()-interatomic_force);

	return energy;

}

template <int dim>
double LennardJones <dim>::LennardJones_Configurational(Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta, double epsilon, double sigma)
{

	Bond <dim> bond;

	this -> atomi = atom_Alpha;
	this -> atomj = atom_Beta;

	double material_distance=bond.MaterialBondDistance(atomi, atomj);
	double material_stretch=bond.MaterialBondStretch(atomi, atomj);

	double sigma0=sigma/material_distance;
	double epsilon0=epsilon/material_distance;

	double term1=sigma0*material_stretch;

	double term2=pow(term1, 6);
	double term3=term2*term2;

	double energy= 4*epsilon0*material_stretch*(term3-term2);

	Point <dim> material_bond_vec=bond.MaterialBondNormal(atomi, atomj);

	Point <dim> configurational_interatomic_force=material_bond_vec*4*epsilon0*(13*term3-7*term2);

	atomi.SetForce(atomi.GetForce()+configurational_interatomic_force);
	atomj.SetForce(atomj.GetForce()-configurational_interatomic_force);

	return energy;

}

template <int dim>
double LennardJones <dim>::LennardJones_Deformational(Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta, double epsilon, double sigma)
{

	Bond <dim> bond;

	this -> atomi = atom_Alpha;
	this -> atomj = atom_Beta;

	double material_distance=bond.MaterialBondDistance(atomi, atomj);
	double spatial_stretch=bond.SpatialBondStretch(atomi, atomj);

	double sigma0=sigma/material_distance;
	double epsilon0=epsilon/material_distance;

	double term1=sigma0/spatial_stretch;
	double term2=pow(term1, 6);
	double term3=term2*term2;

	double energy= 4*epsilon*(term3-term2);

	Point <dim> spatial_bond_vec(0.,0.,0.);

	spatial_bond_vec=bond.SpatialBondNormal(atomi, atomj);

	Point <dim> interatomic_force=spatial_bond_vec*(4*epsilon/spatial_stretch)*(-12*term3+6*term2);

	atomi.SetForce(atomi.GetForce()+interatomic_force);
	atomj.SetForce(atomj.GetForce()-interatomic_force);

	return energy;

}



