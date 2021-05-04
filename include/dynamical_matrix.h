/*
 * dynamical_matrix.h
 *
 *  Created on: Sep 17, 2020
 *      Author: S.Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "matrix.h"

template<int dim> class DynamicalMatrix {
public:
    DynamicalMatrix();
    ~DynamicalMatrix();
    void DynamicalMatrixLJ(int i, double sigma, double epsilon);
    SecondTensor<double> DynamicalMatrixEAM(int i, int j);

private:
    Atoms<dim> *atoms;
};




