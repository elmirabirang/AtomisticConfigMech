/*
 * dynamical_matrix.cpp
 *
 *  Created on: Jan 17, 2021
 *      Author: S.Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include "dynamical_matrix.h"
#include "matrix.h"
#include "atom.h"
#include "bond.h"
#include "math.h"

using namespace std;

template<int dim> void DynamicalMatrix<dim>::DynamicalMatrixLJ(int i, double sigma, double epsilon) {
    SecondTensor <double> dynamical_matrix(dim, dim, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto spatial_bond_vec = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        auto spatial_bond_distance = spatial_bond_vec.norm();
        double term1=4*epsilon;
        double term2=pow((sigma/spatial_bond_distance),12);
        double term3=pow((sigma/spatial_bond_distance),6);
        double term4=1/spatial_bond_distance;
        double term5=1/pow(spatial_bond_distance,2);
        double first_deriv_LJ=term1*((-12*term4)*term2+(6*term4)*term3);
        double second_deriv_LJ=term1*((156*term5)*term2-(56*term5)*term3);
        double spatial_bond_vecX=spatial_bond_vec.GetXCoord();
        double spatial_bond_vecY=spatial_bond_vec.GetYCoord();
        double spatial_bond_vecZ=spatial_bond_vec.GetZCoord();

        SecondTensor <double> identity(3,3,0.0);
        identity(0,0)=1.0;
        identity(1,1)=1.0;
        identity(2,2)=1.0;

        SecondTensor <double> MultMatBondVectors(3,3,0.0);
        MultMatBondVectors(0,0)=spatial_bond_vecX*spatial_bond_vecX;
        MultMatBondVectors(0,1)=spatial_bond_vecX*spatial_bond_vecY;
        MultMatBondVectors(0,2)=spatial_bond_vecX*spatial_bond_vecZ;
        MultMatBondVectors(1,0)=spatial_bond_vecY*spatial_bond_vecX;
        MultMatBondVectors(1,1)=spatial_bond_vecY*spatial_bond_vecY;
        MultMatBondVectors(1,2)=spatial_bond_vecY*spatial_bond_vecZ;
        MultMatBondVectors(2,0)=spatial_bond_vecZ*spatial_bond_vecX;
        MultMatBondVectors(2,1)=spatial_bond_vecZ*spatial_bond_vecY;
        MultMatBondVectors(2,2)=spatial_bond_vecZ*spatial_bond_vecZ;

        dynamical_matrix = first_deriv_LJ * term4 * identity + (second_deriv_LJ - (first_deriv_LJ * term4)) *
                           term5 * MultMatBondVectors;
	)

	double frequency = (0.3333333) * (dynamical_matrix(0, 0) + dynamical_matrix(1, 1) + dynamical_matrix(2, 2));
    this->atoms->setDynamicalMatrix(i, dynamical_matrix);
    this->atoms->setFrequency(i, sqrt(frequency));
}

