/*
 * bond.h
 *
 *  Created on: Nov 17, 2018
 *      Author: S.Elmira Birang.O
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 *
 */

#pragma once

#include <vector>
//---
#include "atom.h"
#include "point.h"
#include "matrix.h"
#include "third_order_tensor.h"

using namespace std;
template <int dim>
class Bond{

public:
	Bond();
	~Bond();

	Point <dim> MaterialBondVec(Atom <dim> atoma, Atom <dim> atomb);
	Point <dim> SpatialBondVec(Atom<dim> atoma, Atom<dim> atomb);
	Point <dim> InitialBondVec(Atom<dim> atoma, Atom<dim> atomb);

	double MaterialBondDistance(Atom<dim> atoma, Atom<dim> atomb);
	double SpatialBondDistance(Atom<dim> atoma, Atom<dim> atomb);
	double InitialBondDistance(Atom <dim> atoma, Atom <dim> atomb);

	Point <dim> MaterialBondNormal(Atom<dim> atoma,Atom<dim> atomb);
	Point <dim> SpatialBondNormal(Atom<dim> atoma, Atom<dim> atomb);

	double MaterialBondStretch(Atom<dim> atoma,Atom<dim> atomb);
	double SpatialBondStretch(Atom<dim> atoma, Atom<dim> atomb);
	double InitialBondStretch (Atom <dim> atoma, Atom <dim> atomb);

	Point <dim> ResultantBondVec_Material(Atom <dim> atoma);
	Point <dim> ResultantBondVec_Spatial(Atom <dim> atoma);

	Point <dim> RotationAxis(Atom <dim> atoma, Atom <dim> atomb);
	double RotationAngle(Atom <dim> atoma, Atom <dim> atomb);
	SecondTensor <double> RotationMatrix(Atom <dim> atoma, Atom <dim> atomb);
	ThirdTensor <double> DerivRotMat_RotAxis(Atom <dim> atoma, Atom <dim> atomb);
	ThirdTensor <double> DerivRotMat_RotAngle(Atom <dim> atoma, Atom <dim> atomb);




private:

	Atom<dim> atomi;
	Atom<dim> atomj;
};
