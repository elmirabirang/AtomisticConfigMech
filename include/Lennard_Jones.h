/*
 * Lennard_Jones.h
 *
 *  Created on: Oct 20, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#pragma once

#include <vector>
#include <string>
//---
#include "atom.h"
#include "bond.h"
#include "point.h"
#include "math.h"
#include "matrix.h"

using namespace std;

template<int dim> class LennardJones{
public:
    LennardJones();
    ~LennardJones();
    double LennardJones_Classic(int i, int j, double epsilon, double sigma);
    double LennardJones_Configurational(int i, int j, double epsilon, double sigma);
    double LennardJones_Deformational(int i, int j, double epsilon, double sigma);

private:
    Atoms<dim> *atoms;
};

