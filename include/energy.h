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

#include "atom.h"
#include "bond.h"

template<int dim> class Energy {
public:
	Energy(Atoms<dim> *atoms);
	Energy(Atoms<dim> *atoms, double sigma, double epsilon);
	~Energy();

	double InteratomicEnergy(int i, int j);
	double CutoffEnergy(int i);
//	double ClusterEnergy(Atom atoma, double clusterR);
	double Material_Interatomic_Engergy(int i, int j);
	double Spatial_Interatomic_Energy(int i, int j);
	double TotPotentialEnergy();
	double atomTotalPotentialEnergy(int i);
	double ConfigInteratomicEnergy(int i, int j);
	double ConfigTotPotentialEnergy();

	double electronDensity(int i, double atomic_number);
//	double embeddingEnergy(vector < Atom <dim>* > atoms_list, double fitting_parameter, double atomic_number);
	double tersoffEnergy(double lambda, double A);

	double morseFunction(double alpha1, double alpha2,
			             double r_01, double r_02, double E1, double E2, double delta,
						 double r_cut, double h, double r_s1, double r_s2, double r_s3,
						 double s1, double s2, double s3);

	double embeddingEnergy(double a, double beta1,
			               double beta2, double r_03, double r_04, double F0, double F2,
						   double q1, double q2, double q3, double q4, double Q1, double Q2, double h, double r_cut);

	double stillingerWeberTwoBody(int i, int j,
			                      double sigma_alphabeta, double epsilon_alphabeta,
                                  double A_alphabeta, double B_alphabeta, double p_alphabeta,
								  double q_alphabeta, double a_alphabeta);

	double stillingerWeberThreeBody(int i, int j, int k, double sigma_AlphaBeta,
                                     double sigma_AlphaGamma, double epsilon_AlphaBetaGamma, double gamma,
                                     double lambda_AlphaBetaGamma, double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	double totalStillingerWeberEnergy(int i,
			                          double sigma_AlphaBeta,
					                  double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
									  double gamma, double lambda_AlphaBetaGamma,
									  double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma,
									  double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
								  	  double q_AlphaBeta);

	double deformationalSWtwoBody(int i, int j,
                                  double sigma_alphabeta, double epsilon_alphabeta,
                                  double A_alphabeta, double B_alphabeta, double I_alphabeta,
			                      double J_alphabeta, double a_alphabeta);

	double deformationalSWthreeBody(int i, int j, int k,
			                        double sigma_AlphaBeta, double sigma_AlphaGamma, double epsilon_AlphaBetaGamma,
									double mu, double eta_AlphaBetaGamma, double cosine_teta0,
                                    double a_AlphaBeta, double a_AlphaGamma);

	double configurationalSWtwoBody(int i, int j,
                                    double sigma_alphabeta, double epsilon_alphabeta,
                                    double A_alphabeta, double B_alphabeta, double I_alphabeta,
			                        double J_alphabeta, double a_alphabeta);

	double configurationalSWthreeBody(int i, int j, int k,
			                          double sigma_AlphaBeta, double sigma_AlphaGamma, double epsilon_AlphaBetaGamma,
									  double mu, double eta_AlphaBetaGamma,
                                      double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

    double LennrdJonesPotential(double sigma, double epsilon, double distance_ab);
    double phaseAverageLJ(double sigma, double epsilon, double boltzman_const, double plancks_const, int num_quad_points);
    double internalEnergyLJ(double temperature, double boltzman_const, double planck_const, double sigma, double epsilon);

private:
	Atoms<dim> *atoms;
    Bond<dim> *bond;
	double sigma;
	double epsilon;
	double cluster_radius;
	double atomic_nember;
	double fitting_parameter;
	double lambda;
	double A;
};
