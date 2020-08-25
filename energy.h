/*
 * energy.h
 *
 *  Created on: Nov 17, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 *
 */

#pragma once

#include "bond.h"
#include "atom.h"
#include <vector>

using namespace std;

template <int dim>
class Energy{
public:

	Energy(double sigma, double epsilon);
	Energy();
	~Energy();

	Bond <dim> bond;

	double InteratomicEnergy(Atom <dim> &atoma, Atom <dim> &atomb);
	double CutoffEnergy(Atom <dim> atoma);
//	double ClusterEnergy(Atom atoma, double clusterR);
	double Material_Interatomic_Engergy(Atom <dim> atoma, Atom <dim> atomb);
	double Spatial_Interatomic_Energy(Atom <dim> atoma, Atom <dim> atomb);
	double TotPotentialEnergy(vector < Atom<dim>* > atoms_list);
	double atomTotalPotentialEnergy(Atom <dim> *atoma);

	double ConfigInteratomicEnergy(Atom <dim> &atoma, Atom <dim> &atomb);
	double ConfigTotPotentialEnergy(vector < Atom<dim>* > atoms_list);

	double electronDensity(Atom <dim> &atoma, double atomic_number);
//	double embeddingEnergy(vector < Atom <dim>* > atoms_list, double fitting_parameter, double atomic_number);
	double tersoffEnergy(vector < Atom <dim>* > atoms_list, double lambda, double A);

	double morseFunction ( vector < Atom <dim>* > atoms_list, double alpha1, double alpha2,
			               double r_01, double r_02, double E1, double E2, double delta,
						   double r_cut, double h, double r_s1, double r_s2, double r_s3,
						   double s1, double s2, double s3);

	double embeddingEnergy( vector <Atom <dim> * > atoms_list, double a, double beta1,
			               double beta2, double r_03, double r_04, double F0, double F2,
						   double q1, double q2, double q3, double q4, double Q1, double Q2, double h, double r_cut);

	double stillingerWeberTwoBody (Atom <dim> atom_alpha, Atom <dim> atom_beta,
			                       double sigma_alphabeta, double epsilon_alphabeta,
                                   double A_alphabeta, double B_alphabeta, double p_alphabeta,
								   double q_alphabeta, double a_alphabeta);

	double stillingerWeberThreeBody (Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta, Atom <dim> &atom_Gamma, double sigma_AlphaBeta,
                                     double sigma_AlphaGamma,double epsilon_AlphaBetaGamma, double gamma, double lambda_AlphaBetaGamma
			                        ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	double totalStillingerWeberEnergy(Atom <dim> *atoma,
			                             double sigma_AlphaBeta,
					                     double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
										 double gamma, double lambda_AlphaBetaGamma,
										 double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma,
										 double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
										 double q_AlphaBeta);

	double deformationalSWtwoBody(Atom <dim> atom_alpha, Atom <dim> atom_beta,
                                  double sigma_alphabeta, double epsilon_alphabeta,
                                  double A_alphabeta, double B_alphabeta, double I_alphabeta,
			                      double J_alphabeta, double a_alphabeta);

	double deformationalSWthreeBody(Atom <dim> atom_Alpha, Atom <dim> atom_Beta, Atom <dim> atom_Gamma,
			                        double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
									double mu, double eta_AlphaBetaGamma
						            ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	double configurationalSWtwoBody(Atom <dim> atom_alpha, Atom <dim> atom_beta,
                                    double sigma_alphabeta, double epsilon_alphabeta,
                                    double A_alphabeta, double B_alphabeta, double I_alphabeta,
			                        double J_alphabeta, double a_alphabeta);

	double configurationalSWthreeBody(Atom <dim> atom_Alpha, Atom <dim> atom_Beta, Atom <dim> atom_Gamma,
			                          double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
									  double mu, double eta_AlphaBetaGamma
						              ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

private:

	Atom <dim> atomi;
	Atom <dim> atomj;
	Atom <dim> atomk;

	vector < Atom<dim>* > atoms;

	double sigma;
	double epsilon;

	double cluster_radius;

	double atomic_nember;
	double fitting_parameter;

	double lambda;
	double A;

};
