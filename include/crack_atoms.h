/*
 * crack_atoms.h
 *
 *  Created on: Jun 11, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nuremberg
 */

#include "atom.h"

using namespace std;

template <int dim>
class CrackAtom {

public:

	CrackAtom ();
	~CrackAtom();

	void setCrackRegion(Atom <dim>* atom);

	double energyCrackAtom( CrackAtom <dim>* crack_atom, double Rcut, double sigma, double epsilon);

	void setCrackBottomAtoms(vector < Atom <dim>* > bottom_atoms);
	vector < Atom <dim>* > getCrackBottomAtoms();

	void setCrackTopAtoms(vector < Atom <dim>* > top_atoms);
	vector < Atom <dim>* > getCrackTopAtoms();


private:

	CrackAtom <dim>* brittle_crack_atom;
	vector < Atom <dim>* > crack_bottom_atoms;
	vector < Atom <dim>* > crack_top_atoms;
	double cut_radius;

};
