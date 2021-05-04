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

#include "atom.h"
#include "point.h"
#include "matrix.h"
#include "third_order_tensor.h"

template<int dim> class Bond {
public:
    Bond(Atoms<dim> *atoms_) : atoms(atoms_) {};
    ~Bond() {};

    Point<dim> MaterialBondVec(int i, int j);
    Point<dim> SpatialBondVec(int i, int j);
    Point<dim> InitialBondVec(int i, int j);

    double MaterialBondDistance(int i, int j);
    double SpatialBondDistance(int i, int j);
    double InitialBondDistance(int i, int j);

    Point<dim> MaterialBondNormal(int i, int j);
    Point<dim> SpatialBondNormal(int i, int j);

    double MaterialBondStretch(int i, int j);
    double SpatialBondStretch(int i, int j);
    double InitialBondStretch(int i, int j);

    Point<dim> ResultantBondVec_Material(int i);
    Point<dim> ResultantBondVec_Spatial(int i);

    Point<dim> RotationAxis(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec);
    double RotationAngle(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec);
    SecondTensor<double> RotationMatrix(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec);
    ThirdTensor<double> DerivRotMat_RotAxis(int i, int j);
    ThirdTensor<double> DerivRotMat_RotAngle(int i, int j);

private:
    Atoms<dim> *atoms;
};
