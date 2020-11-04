/*
 * Stillinger_Weber.cpp
 *
 *  Created on: Oct 19, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#include "atom.h"
#include "bond.h"
#include "point.h"
#include "Stillinger_Weber.h"
#include "math.h"
#include "matrixx.h"
#include <vector>
#include <string>
#include <iostream>

template <int dim>
StillingerWeber <dim>::StillingerWeber()
{}

template <int dim>
StillingerWeber <dim>::~StillingerWeber()
{}

template <int dim>
double StillingerWeber <dim>::StillingerWeberTwoBody_Deformational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
		                                                           double sigma_AlphaBeta, double epsilon_AlphaBeta,
					                                               double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
								                                   double J_AlphaBeta, double a_AlphaBeta)
{

	Bond <dim> bond;

	this -> atomi=atom_alpha;
	this -> atomj=atom_beta;

	double energy=0;

	double spatial_stretch=bond.SpatialBondStretch(atomi,atomj);
	double material_distance=bond.MaterialBondDistance(atomi, atomj);

	double sigma_0=sigma_AlphaBeta/material_distance;
	double epsilon_0=epsilon_AlphaBeta/material_distance;

	double term1=1/(spatial_stretch-a_AlphaBeta*sigma_0);
	double exp_term1=exp(sigma_0*term1);

	double term2=1/spatial_stretch;

	double term3=pow(sigma_0*term2, I_AlphaBeta);
	double term4=pow(sigma_0*term2, J_AlphaBeta);

	energy=-epsilon_0*A_AlphaBeta*(B_AlphaBeta*term3-term4)*exp_term1;

	double coefficient0=I_AlphaBeta*B_AlphaBeta*term2;
	double coefficient1=J_AlphaBeta*term2;

	Point <dim> bond_normal_vec;
	bond_normal_vec=bond.SpatialBondNormal(atom_alpha, atom_beta);

	Point <dim> deformational_twobody_force(0.,0.,0.);

	deformational_twobody_force=bond_normal_vec*energy*
			                    (((coefficient0*term3-coefficient1*term4)/(B_AlphaBeta*term3-term4))
			                    +(term1*term1*sigma_0));

	Point <dim> atom_alpha_position=atom_alpha.GetSpatialPosition();
	Point <dim> atom_beta_position=atom_beta.GetSpatialPosition();

	atom_alpha.SetForce(atom_alpha.GetForce()+deformational_twobody_force);
	atom_beta.SetForce(atom_beta.GetForce()+deformational_twobody_force*-1);


	return energy;

}

template <int dim>
double StillingerWeber <dim>:: StillingerWeberThreeBody_Deformational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta, Atom <dim> &atom_gamma,
		                                                              double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0)
{

	Bond <dim> bond;

	this -> atomi=atom_alpha;
	this -> atomj=atom_beta;
	this -> atomk=atom_gamma;
	double mu=1.2;
	double a=1.8;
	double sigma=2.0951;

	double spatial_distance_AlphaBeta=bond.SpatialBondDistance(atom_alpha, atom_beta);
	double spatial_distance_AlphaGamma=bond.SpatialBondDistance(atom_alpha, atom_gamma);

	double terma=exp((mu*sigma)/(spatial_distance_AlphaBeta-a*sigma)
			         +(mu*sigma)/(spatial_distance_AlphaGamma-a*sigma));

	Point <dim> material_normal_AlphaBeta=bond.MaterialBondNormal(atomi,atomj);
	double material_normal_AlphaBeta_x=material_normal_AlphaBeta.GetXCoord();
	double material_normal_AlphaBeta_y=material_normal_AlphaBeta.GetYCoord();
	double material_normal_AlphaBeta_z=material_normal_AlphaBeta.GetZCoord();

	Point <dim> material_normal_AlphaGamma=bond.MaterialBondNormal(atomi,atomk);
	double material_normal_AlphaGamma_x=material_normal_AlphaGamma.GetXCoord();
	double material_normal_AlphaGamma_y=material_normal_AlphaGamma.GetYCoord();
	double material_normal_AlphaGamma_z=material_normal_AlphaGamma.GetZCoord();

	double cos_Theta_AlphaBetaGamma=material_normal_AlphaBeta_x*material_normal_AlphaGamma_x+
			                        material_normal_AlphaBeta_y*material_normal_AlphaGamma_y+
									material_normal_AlphaBeta_z*material_normal_AlphaGamma_z;

	Point <dim> spatial_normal_AlphaBeta=bond.SpatialBondNormal(atomi,atomj);
	double spatial_normal_AlphaBeta_x=spatial_normal_AlphaBeta.GetXCoord();
	double spatial_normal_AlphaBeta_y=spatial_normal_AlphaBeta.GetYCoord();
	double spatial_normal_AlphaBeta_z=spatial_normal_AlphaBeta.GetZCoord();

	Point <dim> spatial_normal_AlphaGamma=bond.SpatialBondNormal(atomi,atomk);
	double spatial_normal_AlphaGamma_x=spatial_normal_AlphaGamma.GetXCoord();
	double spatial_normal_AlphaGamma_y=spatial_normal_AlphaGamma.GetYCoord();
	double spatial_normal_AlphaGamma_z=spatial_normal_AlphaGamma.GetZCoord();

	double cos_theta_AlphaBetaGamma=spatial_normal_AlphaBeta_x*spatial_normal_AlphaGamma_x+
			                        spatial_normal_AlphaBeta_y*spatial_normal_AlphaGamma_y+
									spatial_normal_AlphaBeta_z*spatial_normal_AlphaGamma_z;

	double inv_cos_Theta_AlphaBetaGamma=1/cos_Theta_AlphaBetaGamma;

	double omega_AlphaBetaGamma=cos_theta_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double omega0_AlphaBetaGamma=cosine_teta0*inv_cos_Theta_AlphaBetaGamma;

	double epsilon0_AlphaBetaGamma=epsilon_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double eta0_AlphaBetaGamma=eta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma;

	double diff_omegas=omega_AlphaBetaGamma-omega0_AlphaBetaGamma;

	double force_coefficient=2*epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma*diff_omegas;

	double energy=epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma
			      *diff_omegas*diff_omegas;

	SecondTensor <double> project_tensor_AlphaBeta (3,3,0);
	SecondTensor <double> project_tensor_AlphaGamma (3,3,0);

	Point <dim> bond_normal_vec_AlphaBeta=bond.SpatialBondNormal(atomi,atomj);
	double spatial_bond_distance_AlphaBeta=bond.SpatialBondDistance(atomi,atomj);

	double normal_AlphaBeta_x=bond_normal_vec_AlphaBeta.GetXCoord();
	double normal_AlphaBeta_y=bond_normal_vec_AlphaBeta.GetYCoord();
	double normal_AlphaBeta_z=bond_normal_vec_AlphaBeta.GetZCoord();
	double inv_distance_AlphaBeta=1/spatial_bond_distance_AlphaBeta;

	Point <dim> bond_normal_vec_AlphaGamma=bond.SpatialBondNormal(atomi,atomk);
	double spatial_bond_distance_AlphaGamma=bond.SpatialBondDistance(atomi,atomk);

	double normal_AlphaGamma_x=bond_normal_vec_AlphaGamma.GetXCoord();
	double normal_AlphaGamma_y=bond_normal_vec_AlphaGamma.GetYCoord();
	double normal_AlphaGamma_z=bond_normal_vec_AlphaGamma.GetZCoord();
	double inv_distance_AlphaGamma=1/spatial_bond_distance_AlphaGamma;

	Point <dim> term1(0.,0.,0.);
	Point <dim> term2(0.,0.,0.);

	term1.SetXCoord((inv_distance_AlphaGamma
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			         *normal_AlphaBeta_x);

	term1.SetYCoord((inv_distance_AlphaGamma
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			         *normal_AlphaBeta_y);

	term1.SetZCoord((inv_distance_AlphaGamma
			        -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			        *normal_AlphaBeta_z);

	term2.SetXCoord((inv_distance_AlphaBeta
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			         *normal_AlphaGamma_x);

	term2.SetYCoord((inv_distance_AlphaBeta
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			         *normal_AlphaGamma_y);

	term2.SetZCoord((inv_distance_AlphaBeta
			        -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			        *normal_AlphaGamma_z);

//	//Alpha Beta
//	project_tensor_AlphaBeta(0,0)=(1-(normal_AlphaBeta_x*normal_AlphaBeta_x))*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(0,1)=normal_AlphaBeta_x*normal_AlphaBeta_y*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(0,2)=normal_AlphaBeta_x*normal_AlphaBeta_z*inv_distance_AlphaBeta;
//
//	project_tensor_AlphaBeta(1,0)=normal_AlphaBeta_y*normal_AlphaBeta_x*inv_distance_AlphaBeta ;
//	project_tensor_AlphaBeta(1,1)=(1-(normal_AlphaBeta_y*normal_AlphaBeta_y))*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(1,2)=normal_AlphaBeta_y*normal_AlphaBeta_z*inv_distance_AlphaBeta;
//
//	project_tensor_AlphaBeta(2,0)=normal_AlphaBeta_z*normal_AlphaBeta_x*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(2,1)=normal_AlphaBeta_z*normal_AlphaBeta_y*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(2,2)=(1-(normal_AlphaBeta_z*normal_AlphaBeta_z))*inv_distance_AlphaBeta;
//
//    //Alpha Gamma
//	project_tensor_AlphaGamma(0,0)=(1-(normal_AlphaGamma_x*normal_AlphaGamma_x))*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(0,1)=normal_AlphaGamma_x*normal_AlphaGamma_y*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(0,2)=normal_AlphaGamma_x*normal_AlphaGamma_z*inv_distance_AlphaGamma;
//
//	project_tensor_AlphaGamma(1,0)=normal_AlphaGamma_y*normal_AlphaGamma_x*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(1,1)=(1-(normal_AlphaGamma_y*normal_AlphaGamma_y))*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(1,2)=normal_AlphaGamma_y*normal_AlphaGamma_z*inv_distance_AlphaGamma;
//
//	project_tensor_AlphaGamma(2,0)=normal_AlphaGamma_z*normal_AlphaGamma_x*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(2,1)=normal_AlphaGamma_z*normal_AlphaGamma_y*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(2,2)=(1-(normal_AlphaGamma_z*normal_AlphaGamma_z))*inv_distance_AlphaGamma;
//
//	Point <dim> term1 (0.,0.,0.);
//
//	term1.SetXCoord(normal_AlphaGamma_x*
//			        project_tensor_AlphaBeta(0,0)+
//					normal_AlphaGamma_y*
//			        project_tensor_AlphaBeta(0,1)+
//					normal_AlphaGamma_z*
//					project_tensor_AlphaBeta(0,2));
//
//	term1.SetYCoord(normal_AlphaGamma_x*
//	                project_tensor_AlphaBeta(1,0)+
//					normal_AlphaGamma_y*
//	                project_tensor_AlphaBeta(1,1)+
//					normal_AlphaGamma_z*
//			        project_tensor_AlphaBeta(1,2));
//
//	term1.SetZCoord(normal_AlphaGamma_x*
//                    project_tensor_AlphaBeta(2,0)+
//					normal_AlphaGamma_y*
//                    project_tensor_AlphaBeta(2,1)+
//					normal_AlphaGamma_z*
//	                project_tensor_AlphaBeta(2,2));
//
//	Point <dim> term2 (0.,0.,0.);
//
//	term2.SetXCoord(normal_AlphaBeta_x*
//			        project_tensor_AlphaGamma(0,0)+
//					normal_AlphaBeta_y*
//			        project_tensor_AlphaGamma(0,1)+
//					normal_AlphaBeta_z*
//					project_tensor_AlphaGamma(0,2));
//
//	term2.SetYCoord(normal_AlphaBeta_x*
//	                project_tensor_AlphaGamma(1,0)+
//					normal_AlphaBeta_y*
//	                project_tensor_AlphaGamma(1,1)+
//					normal_AlphaBeta_z*
//			        project_tensor_AlphaGamma(1,2));
//
//	term2.SetZCoord(normal_AlphaBeta_x*
//                    project_tensor_AlphaGamma(2,0)+
//					normal_AlphaBeta_y*
//                    project_tensor_AlphaGamma(2,1)+
//					normal_AlphaBeta_z*
//	                project_tensor_AlphaGamma(2,2));

	Point <dim> threeBody_DeformForce_Alpha=(term1+term2)*force_coefficient ;
	Point <dim> threeBody_DeformForce_Beta=term1*force_coefficient ;
	Point <dim> threeBody_DeformForce_Gamma=term2*force_coefficient;

	atom_alpha.SetForce(atom_alpha.GetForce()+threeBody_DeformForce_Alpha);
	atom_beta.SetForce(atom_beta.GetForce()-threeBody_DeformForce_Beta);
	atom_gamma.SetForce(atom_gamma.GetForce()-threeBody_DeformForce_Gamma);

	return energy;

}

template <int dim>
double StillingerWeber <dim>::StillingerWeberTwoBody_Configurational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta,
		                                                             double sigma_AlphaBeta, double epsilon_AlphaBeta,
					                                                 double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
								                                     double J_AlphaBeta, double a_AlphaBeta)
{

	Bond <dim> bond;

	this -> atomi=atom_alpha;
	this -> atomj=atom_beta;

	double energy=0;

	double material_stretch=bond.MaterialBondStretch(atomi, atomj);
	double material_distance=bond.MaterialBondDistance(atomi, atomj);

	double sigma_0=sigma_AlphaBeta/material_distance;
	double epsilon_0=epsilon_AlphaBeta/material_distance;

	double term1=pow(sigma_0*material_stretch, I_AlphaBeta);
	double term2=pow(sigma_0*material_stretch, J_AlphaBeta);

	double term3=1/(1-material_stretch*a_AlphaBeta*sigma_0);
	double term4=exp(sigma_0*material_stretch*term3);

	energy=-epsilon_0*A_AlphaBeta*material_stretch*(B_AlphaBeta*term1-term2)*term4;

	//term1=pow(sigma_0*material_stretch, I_AlphaBeta)
	//term2=pow(sigma_0*material_stretch, J_AlphaBeta)
    //term3=1/(1-material_stretch*a_AlphaBeta*sigma_0);

	Point <dim> bond_normal_vec=bond.MaterialBondNormal(atomi,atomj);

	Point <dim> configurational_twobody_force(0.,0.,0.);

	double inv_material_stretch=1/material_stretch;

	double coefficient0=B_AlphaBeta*I_AlphaBeta;
	double coefficient1=B_AlphaBeta*material_stretch;

	configurational_twobody_force=bond_normal_vec*energy*
			                      (inv_material_stretch+
			                       ((coefficient0*term1-J_AlphaBeta*term2)/
			                        (coefficient1*term1-material_stretch*term2))
								    +sigma_0*term3*term3);

	atom_alpha.SetConfigForce(atom_alpha.GetConfigForce()+configurational_twobody_force);

	atom_beta.SetConfigForce(atom_beta.GetConfigForce()-configurational_twobody_force);

	return energy;

}

template <int dim>
double StillingerWeber <dim>:: StillingerWeberThreeBody_Configurational(Atom <dim> &atom_alpha, Atom <dim> &atom_beta, Atom <dim> &atom_gamma,
		                                                                double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0)
{

	Bond <dim> bond;

	this -> atomi=atom_alpha;
	this -> atomj=atom_beta;
	this -> atomk=atom_gamma;

	Point <dim> material_normal_AlphaBeta=bond.MaterialBondNormal(atomi,atomj);
	double material_normal_AlphaBeta_x=material_normal_AlphaBeta.GetXCoord();
	double material_normal_AlphaBeta_y=material_normal_AlphaBeta.GetYCoord();
	double material_normal_AlphaBeta_z=material_normal_AlphaBeta.GetZCoord();

	Point <dim> material_normal_AlphaGamma=bond.MaterialBondNormal(atomi,atomk);
	double material_normal_AlphaGamma_x=material_normal_AlphaGamma.GetXCoord();
	double material_normal_AlphaGamma_y=material_normal_AlphaGamma.GetYCoord();
	double material_normal_AlphaGamma_z=material_normal_AlphaGamma.GetZCoord();

	double cos_Theta_AlphaBetaGamma=material_normal_AlphaBeta_x*material_normal_AlphaGamma_x+
			                        material_normal_AlphaBeta_y*material_normal_AlphaGamma_y+
									material_normal_AlphaBeta_z*material_normal_AlphaGamma_z;

	Point <dim> spatial_normal_AlphaBeta=bond.SpatialBondNormal(atomi,atomj);
	double spatial_normal_AlphaBeta_x=spatial_normal_AlphaBeta.GetXCoord();
	double spatial_normal_AlphaBeta_y=spatial_normal_AlphaBeta.GetYCoord();
	double spatial_normal_AlphaBeta_z=spatial_normal_AlphaBeta.GetZCoord();

	Point <dim> spatial_normal_AlphaGamma=bond.SpatialBondNormal(atomi,atomk);
	double spatial_normal_AlphaGamma_x=spatial_normal_AlphaGamma.GetXCoord();
	double spatial_normal_AlphaGamma_y=spatial_normal_AlphaGamma.GetYCoord();
	double spatial_normal_AlphaGamma_z=spatial_normal_AlphaGamma.GetZCoord();

	double cos_theta_AlphaBetaGamma=spatial_normal_AlphaBeta_x*spatial_normal_AlphaGamma_x+
			                        spatial_normal_AlphaBeta_y*spatial_normal_AlphaGamma_y+
									spatial_normal_AlphaBeta_z*spatial_normal_AlphaGamma_z;

	double Omega_AlphaBetaGamma=cos_Theta_AlphaBetaGamma/cos_theta_AlphaBetaGamma;

	double inv_cos_Theta_AlphaBetaGamma=1/cos_Theta_AlphaBetaGamma;

	double epsilon0_AlphaBetaGamma=epsilon_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double omega0_AlphaBetaGamma=cosine_teta0*inv_cos_Theta_AlphaBetaGamma;
	double eta0_AlphaBetaGamma=eta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma;

	double inv_Omega_AlphaBetaGamma=1/Omega_AlphaBetaGamma;

	double diff_omegas=inv_Omega_AlphaBetaGamma-omega0_AlphaBetaGamma;

	double coefficient=epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma;

	double force_coefficient=-2*coefficient*inv_Omega_AlphaBetaGamma*diff_omegas
			                  +coefficient*diff_omegas*diff_omegas;

	cout << "force_coefficient: " << force_coefficient<< endl;

	double energy=coefficient*Omega_AlphaBetaGamma*diff_omegas*diff_omegas;

//	SecondTensor <double> project_tensor_AlphaBeta (3,3,0);
//	SecondTensor <double> project_tensor_AlphaGamma (3,3,0);

	Point <dim> bond_normal_vec_AlphaBeta=bond.MaterialBondNormal(atomi,atomj);
	double material_bond_distance_AlphaBeta=bond.MaterialBondDistance(atomi,atomj);

	double normal_AlphaBeta_x=bond_normal_vec_AlphaBeta.GetXCoord();
	double normal_AlphaBeta_y=bond_normal_vec_AlphaBeta.GetYCoord();
	double normal_AlphaBeta_z=bond_normal_vec_AlphaBeta.GetZCoord();
	double inv_distance_AlphaBeta=1/material_bond_distance_AlphaBeta;

	Point <dim> bond_normal_vec_AlphaGamma=bond.MaterialBondNormal(atomi,atomk);
	double material_bond_distance_AlphaGamma=bond.MaterialBondDistance(atomi,atomk);

	double normal_AlphaGamma_x=bond_normal_vec_AlphaGamma.GetXCoord();
	double normal_AlphaGamma_y=bond_normal_vec_AlphaGamma.GetYCoord();
	double normal_AlphaGamma_z=bond_normal_vec_AlphaGamma.GetZCoord();
	double inv_distance_AlphaGamma=1/material_bond_distance_AlphaGamma;

	Point <dim> term1(0.,0.,0.);
	Point <dim> term2(0.,0.,0.);

	term1.SetXCoord((inv_distance_AlphaGamma
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			         *normal_AlphaBeta_x);

	term1.SetYCoord((inv_distance_AlphaGamma
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			         *normal_AlphaBeta_y);

	term1.SetZCoord((inv_distance_AlphaGamma
			        -(cos_theta_AlphaBetaGamma*inv_distance_AlphaBeta))
			        *normal_AlphaBeta_z);

	term2.SetXCoord((inv_distance_AlphaBeta
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			         *normal_AlphaGamma_x);

	term2.SetYCoord((inv_distance_AlphaBeta
			         -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			         *normal_AlphaGamma_y);

	term2.SetZCoord((inv_distance_AlphaBeta
			        -(cos_theta_AlphaBetaGamma*inv_distance_AlphaGamma))
			        *normal_AlphaGamma_z);

//	//Alpha Beta
//	project_tensor_AlphaBeta(0,0)=(1-(normal_AlphaBeta_x*normal_AlphaBeta_x))*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(0,1)=normal_AlphaBeta_x*normal_AlphaBeta_y*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(0,2)=normal_AlphaBeta_x*normal_AlphaBeta_z*inv_distance_AlphaBeta;
//
//	project_tensor_AlphaBeta(1,0)=normal_AlphaBeta_y*normal_AlphaBeta_x*inv_distance_AlphaBeta ;
//	project_tensor_AlphaBeta(1,1)=(1-(normal_AlphaBeta_y*normal_AlphaBeta_y))*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(1,2)=normal_AlphaBeta_y*normal_AlphaBeta_z*inv_distance_AlphaBeta;
//
//	project_tensor_AlphaBeta(2,0)=normal_AlphaBeta_z*normal_AlphaBeta_x*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(2,1)=normal_AlphaBeta_z*normal_AlphaBeta_y*inv_distance_AlphaBeta;
//	project_tensor_AlphaBeta(2,2)=(1-(normal_AlphaBeta_z*normal_AlphaBeta_z))*inv_distance_AlphaBeta;
//
//    //Alpha Gamma
//	project_tensor_AlphaGamma(0,0)=(1-(normal_AlphaGamma_x*normal_AlphaGamma_x))*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(0,1)=normal_AlphaGamma_x*normal_AlphaGamma_y*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(0,2)=normal_AlphaGamma_x*normal_AlphaGamma_z*inv_distance_AlphaGamma;
//
//	project_tensor_AlphaGamma(1,0)=normal_AlphaGamma_y*normal_AlphaGamma_x*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(1,1)=(1-(normal_AlphaGamma_y*normal_AlphaGamma_y))*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(1,2)=normal_AlphaGamma_y*normal_AlphaGamma_z*inv_distance_AlphaGamma;
//
//	project_tensor_AlphaGamma(2,0)=normal_AlphaGamma_z*normal_AlphaGamma_x*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(2,1)=normal_AlphaGamma_z*normal_AlphaGamma_y*inv_distance_AlphaGamma;
//	project_tensor_AlphaGamma(2,2)=(1-(normal_AlphaGamma_z*normal_AlphaGamma_z))*inv_distance_AlphaGamma;
//
//	Point <dim> term1 (0.,0.,0.);
//
//	term1.SetXCoord(normal_AlphaGamma_x*
//			        (project_tensor_AlphaBeta(0,0)+
//			         project_tensor_AlphaBeta(0,1)+
//					 project_tensor_AlphaBeta(0,2)));
//
//	term1.SetYCoord(normal_AlphaGamma_y*
//	                (project_tensor_AlphaBeta(1,0)+
//	                 project_tensor_AlphaBeta(1,1)+
//			         project_tensor_AlphaBeta(1,2)));
//
//	term1.SetZCoord(normal_AlphaGamma_z*
//                    (project_tensor_AlphaBeta(2,0)+
//                     project_tensor_AlphaBeta(2,1)+
//	                 project_tensor_AlphaBeta(2,2)));
//
//	Point <dim> term2 (0.,0.,0.);
//
//	term2.SetXCoord(normal_AlphaBeta_x*
//			        (project_tensor_AlphaGamma(0,0)+
//			         project_tensor_AlphaGamma(0,1)+
//					 project_tensor_AlphaGamma(0,2)));
//
//	term2.SetYCoord(normal_AlphaBeta_y*
//	                (project_tensor_AlphaGamma(1,0)+
//	                 project_tensor_AlphaGamma(1,1)+
//			         project_tensor_AlphaGamma(1,2)));
//
//	term2.SetZCoord(normal_AlphaBeta_z*
//                    (project_tensor_AlphaGamma(2,0)+
//                     project_tensor_AlphaGamma(2,1)+
//	                 project_tensor_AlphaGamma(2,2)));

	Point <dim> threeBody_ConfigForce_Alpha=(term1+term2)*force_coefficient;
	Point <dim> threeBody_ConfigForce_Beta=term1*force_coefficient;
	Point <dim> threeBody_ConfigForce_Gamma=term2*force_coefficient;

	atom_alpha.SetConfigForce(atom_alpha.GetConfigForce()+threeBody_ConfigForce_Alpha);
	atom_beta.SetConfigForce(atom_beta.GetConfigForce()-threeBody_ConfigForce_Beta);
	atom_gamma.SetConfigForce(atom_gamma.GetConfigForce()-threeBody_ConfigForce_Gamma);

	return (energy);

}
