/*
 * force.h
 *
 *  Created on: Nov 16, 2018
 *      Author: Elmira Birang
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#pragma once
#include "atom.h"
#include "bond.h"
#include "point.h"
#include <vector>
#include <string>

using namespace std;

template <int dim>
class Force{

public:

	Force();
	~Force();
	Force(double sigma, double epsilon);

	Point <dim> BondForce(Atom <dim> atoma, Atom <dim>atomb);
	Point <dim> MaterialBondForce(Atom <dim> atoma, Atom <dim> atomb);
	Point <dim> CriticalMaterialBondForce(Atom <dim> atoma, Atom <dim> atomb, double cutR);
	Point <dim> SpatialBondForce(Atom <dim> atoma, Atom <dim> atomb);
	Point <dim> ConfigurationalForce(Atom<dim> *atom);
	Point <dim> ResultantForce(Atom <dim>* atoma);
	double DerivMaterialEnergy(Atom <dim> atoma, Atom <dim> atomb);
	Point <dim> ResultantConfigForce(Atom<dim> *atom);

	Point <dim> derivativeMorseFunction( Atom <dim>* atom, double alpha1, double alpha2,
			                             double r_01, double r_02, double E1, double E2, double delta,
						                 double r_cut, double h, double r_s1, double r_s2, double r_s3,
						                 double s1, double s2, double s3);

	Point <dim> embeddingForce(Atom <dim> * atom, double a, double beta1,
                               double beta2, double r_03, double r_04, double F0, double F2,
			                   double q1, double q2, double q3, double q4, double Q1, double Q2,
							   double h, double r_cut);

	double electronDensity (Atom<dim>* atoma, double a, double beta1,
	                        double beta2, double r_03, double r_04,
	                        double h, double r_cut);

	Point <dim> configForceEamPair(Atom <dim>* atom, double alpha1_mod, double alpha2_mod,
                                   double r_01_mod, double r_02_mod, double E1_mod, double E2_mod, double delta_mod,
                                   double r_cut_mod, double h_mod, double r_s1_mod, double r_s2_mod, double r_s3_mod,
                                   double s1_mod, double s2_mod, double s3_mod);

	Point <dim> configForceEamEmbedding(Atom <dim> * atom, double a_mod, double beta1_mod,
                                        double beta2_mod, double r_03_mod, double r_04_mod, double F0_mod, double F2_mod,
                                        double q1_mod, double q2_mod, double q3_mod, double q4_mod, double Q1_mod, double Q2_mod,
			                            double h_mod, double r_cut_mod);

	Point <dim> StillingerWeberThreeBodyForceH2H3(Atom <dim> &atom_Beta, Atom <dim> &atom_Alpha, Atom <dim> &atom_Gamma, double sigma_AlphaBeta,
                                                  double sigma_AlphaGamma,double epsilon_AlphaBetaGamma, double gamma, double lambda_AlphaBetaGamma
                                                 ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> StillingerWeberThreeBodyForceH1(Atom <dim> *atom_Beta, Atom <dim> *atom_Alpha, Atom <dim> *atom_Gamma, double sigma_AlphaBeta,
                                                  double sigma_AlphaGamma,double epsilon_AlphaBetaGamma, double gamma, double lambda_AlphaBetaGamma
                                                 ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> StillingerWeberTwoBodyForce(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
                                            double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                            double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
			                                double q_AlphaBeta, double a_AlphaBeta);

	Point <dim> ResultantSWTwoBodyForce(Atom <dim> *atom, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                        double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                        double q_AlphaBeta, double a_AlphaBeta);

	Point <dim> ResultantSWThreeBodyForce(Atom <dim> *atoma, double sigma_AlphaBeta,
                                          double sigma_AlphaGamma,double epsilon_AlphaBetaGamma, double gamma, double lambda_AlphaBetaGamma
                                         ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> deformationalBondSWForce(Atom <dim> &atom_alpha,Atom <dim> &atom_beta,
			                             double sigma_AlphaBeta, double epsilon_AlphaBeta,
			                             double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
						                 double q_AlphaBeta, double a_AlphaBeta);

	Point <dim> deformationalcpSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
			                                double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
										    double mu, double eta_AlphaBetaGamma,
                                            double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> deformationalprSW3BodyForce(Atom <dim> &atom_alpha, Atom <dim> &atom_beta, Atom <dim> &atom_gamma,
			                                double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
										    double mu, double eta_AlphaBetaGamma,
                                            double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> configurationalBondSWForce(Atom <dim> &atom_alpha,Atom <dim> &atom_beta,
			                               double sigma_AlphaBeta, double epsilon_AlphaBeta,
			                               double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
						                   double q_AlphaBeta, double a_AlphaBeta);

	Point <dim> configurationalcpSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
			                                  double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
										      double mu, double eta_AlphaBetaGamma,
                                              double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> configurationalprSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
			                                  double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
										      double mu, double eta_AlphaBetaGamma,
                                              double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> ResultantConfigSWThreeBodyForce( Atom <dim> *atoma, double sigma_AlphaBeta,
                                                 double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
                                                 double Gamma, double lambda_AlphaBetaGamma
                                                ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma);

	Point <dim> ResultantConfigSWTwoBodyForce(Atom <dim> *atoma, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                              double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                              double q_AlphaBeta, double a_AlphaBeta);




	double interatomicConfigForceValue(Atom <dim> atoma, Atom <dim> atomb);

	double EnergyRelease(Atom<dim> *atoma, double delta_x);

	Bond <dim> bond;



private:

	Atom <dim> atomi;
	Atom <dim> atomj;
	Atom <dim> atomk;

	double sigma;
	double epsilon;
	double cutR;
	double atopmic_number;
	double fitting_parameter;
	double lambda;
	double A;
	double delta_x;

	vector< Atom<dim> > atoms;

};
