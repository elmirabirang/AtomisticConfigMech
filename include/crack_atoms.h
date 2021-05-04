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

template<int dim> class CrackAtom {
public:
	CrackAtom ();
	~CrackAtom();
	void setCrackRegion(int i);
	double energyCrackAtom(CrackAtom<dim> *crack_atom, double Rcut, double sigma, double epsilon);
	void setCrackBottomAtoms(vector<int> bottom_atoms);
	vector<int> getCrackBottomAtoms();
	void setCrackTopAtoms(vector<int> top_atoms);
	vector<int> getCrackTopAtoms();

private:
	CrackAtom<dim> *brittle_crack_atom;
	vector<int> crack_bottom_atoms;
	vector<int> crack_top_atoms;
	double cut_radius;
};
