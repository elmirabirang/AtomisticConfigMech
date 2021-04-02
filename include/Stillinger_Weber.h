/*
 * Stillinger_Weber.h
 *
 *  Created on: Oct 19, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#pragma once

#include <vector>
#include <string>
//---
#include "atom.h"
#include "bond.h"
#include "point.h"

using namespace std;

template <int dim>
class StillingerWeber{

public:

	StillingerWeber();
	~StillingerWeber();

	double StillingerWeberTwoBody_Deformational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
			                                    double sigma_AlphaBeta, double epsilon_AlphaBeta,
						                        double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
									            double J_AlphaBeta, double a_AlphaBeta);

	double StillingerWeberThreeBody_Deformational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta, Atom <dim> &atom_gamma,
			                                      double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma,
												  double cosine_teta0);

	double StillingerWeberTwoBody_Configurational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
			                                      double sigma_AlphaBeta, double epsilon_AlphaBeta,
						                          double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
									              double J_AlphaBeta, double a_AlphaBeta);

	double StillingerWeberThreeBody_Configurational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
			                                        Atom <dim> &atom_gamma,double epsilon_AlphaBetaGamma,
													double eta_AlphaBetaGamma, double cosine_teta0);


private:

	Atom <dim> atomi;
	Atom <dim> atomj;
	Atom <dim> atomk;

};

