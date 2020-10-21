/*
 * Lennard_Jones.h
 *
 *  Created on: Oct 20, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include "atom.h"
#include "bond.h"
#include "point.h"
#include "math.h"
#include "matrixx.h"
#include <vector>
#include <string>

using namespace std;

template <int dim>
class LennardJones{

	public:

		LennardJones();
		~LennardJones();

		double LennardJones_Classic(Atom <dim> &atom_Alpha, Atom <dim> &atom_beta, double epsilon, double sigma);

		double LennardJones_Configurational(Atom <dim> &atom_Alpha, Atom <dim> &atom_beta, double epsilon, double sigma);

		double LennardJones_Deformational(Atom <dim> &atom_Alpha, Atom <dim> &atom_beta, double epsilon, double sigma);

	private:

		Atom <dim> atomi;
		Atom <dim> atomj;

};

