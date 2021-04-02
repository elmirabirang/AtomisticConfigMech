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
#include <vector>

using namespace std;

template <int dim>
class DynamicalMatrix
{

    public:
	DynamicalMatrix();
	~DynamicalMatrix();

	void DynamicalMatrixLJ(Atom <dim>* atom, double sigma, double epsilon);

	SecondTensor<double> DynamicalMatrixEAM(Atom <dim> &atoma, Atom <dim> &atomb);

    private:

	Atom <dim> atomi;
	Atom <dim> atomj;

};




