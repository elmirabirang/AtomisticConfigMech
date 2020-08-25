/*
 * test_optimizer.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "point.h"
#include "energy.h"
#include "force.h"
#include "optimization.h"
#include "optimizer.cpp"

#include <stdlib.h>
#include <stdio.h>

using namespace alglib;
int main{

//update the material position inside load step loop

    vector < Atom <2> > atoms;
    atoms=UnrelaxedConfigGenerator <2>(2,1,0,1.1,2.5);

    vector < Atom <2> > GetAtoms()
		{
    	return(atoms);
		}

    Energy <2> GetEnergy()
    {
    	return (Energy <2> energy(1.0,1.0));
    }

    Force <2> GetForce()
	{
        return(Force <2> force(1.0,1.0));
	}

    typedef typename vector < Atom <2> >::iterator At;
    vector <double> material_position;

    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {
        Point <2> mpositon;
        mpositon=*atoms.GetMaterialPosition();
        material_position.push_back(mposition.GetXCoord());
        material_position.push_back(mposition.GetYCoord());

    }

    real_1d_array x = "material_position";
    double epsg = 0.0000000001;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0.1;
    ae_int_t maxits = 0;
    minlbfgsstate state;
    minlbfgsreport rep;

    minlbfgscreate(1, x, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    minlbfgssetstpmax(state, stpmax);
    alglib::minlbfgsoptimize(state, AtomisticOptimizer);
    minlbfgsresults(state, x, rep);

    printf("%s\n", x.tostring(2).c_str());

	return 0;

}
