/*
 * Stillinger_Weber.cpp
 *
 *  Created on: Oct 19, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#include <vector>
#include <string>
#include <iostream>
//---
#include "atom.h"
#include "bond.h"
#include "point.h"
#include "Stillinger_Weber.h"
#include "math.h"
#include "matrix.h"

template<int dim> StillingerWeber<dim>::StillingerWeber() {}
template<int dim> StillingerWeber<dim>::~StillingerWeber() {}

template<int dim> double StillingerWeber<dim>::StillingerWeberTwoBody_Deformational(
    int i, int j, double sigma_AlphaBeta, double epsilon_AlphaBeta,
    double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta, double J_AlphaBeta, double a_AlphaBeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto bond_normal_vec = spatial_delta / spatial_distance;
    double spatial_stretch = spatial_distance / material_distance;
    double sigma_0=sigma_AlphaBeta / material_distance;
    double epsilon_0=epsilon_AlphaBeta / material_distance;
    double term1=1/(spatial_stretch-a_AlphaBeta*sigma_0);
    double exp_term1=exp(sigma_0*term1);
    double term2=1/spatial_stretch;
    double term3=pow(sigma_0*term2, I_AlphaBeta);
    double term4=pow(sigma_0*term2, J_AlphaBeta);
    double coefficient0 = I_AlphaBeta * B_AlphaBeta * term2;
    double coefficient1 = J_AlphaBeta * term2;
    double energy = -epsilon_0 * A_AlphaBeta * (B_AlphaBeta * term3 - term4) * exp_term1;

    Point <dim> deformational_twobody_force(0.0, 0.0, 0.0);
    deformational_twobody_force = bond_normal_vec * energy *
                                  (((coefficient0 * term3 - coefficient1 * term4) / (B_AlphaBeta * term3 - term4)) +
                                  (term1 * term1 * sigma_0));

    this->atoms->setForce(i, this->atoms->getForce(i) + deformational_twobody_force);
    this->atoms->setForce(j, this->atoms->getForce(j) - deformational_twobody_force);
    return energy;
}

template<int dim> double StillingerWeber<dim>::StillingerWeberThreeBody_Deformational(
    int i, int j, int k,
    double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0) {

	double mu = 1.2;
	double a = 1.8;
	double sigma = 2.0951;
    auto spatial_delta_AlphaBeta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_delta_AlphaGamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spatial_distance_AlphaBeta = spatial_delta_AlphaBeta.norm();
    auto spatial_distance_AlphaGamma = spatial_delta_AlphaGamma.norm();
    auto material_delta_AlphaBeta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_delta_AlphaGamma = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(k);
    auto material_distance_AlphaBeta = material_delta_AlphaBeta.norm();
    auto material_distance_AlphaGamma = material_delta_AlphaGamma.norm();
	double terma = exp((mu * sigma) / (spatial_distance_AlphaBeta - a * sigma) +
                       (mu * sigma) / (spatial_distance_AlphaGamma - a * sigma));

	Point<dim> material_normal_AlphaBeta = material_delta_AlphaBeta / material_distance_AlphaBeta;
	double material_normal_AlphaBeta_x = material_normal_AlphaBeta.GetXCoord();
	double material_normal_AlphaBeta_y = material_normal_AlphaBeta.GetYCoord();
	double material_normal_AlphaBeta_z = material_normal_AlphaBeta.GetZCoord();

	Point<dim> material_normal_AlphaGamma = material_delta_AlphaGamma / material_distance_AlphaGamma;
	double material_normal_AlphaGamma_x = material_normal_AlphaGamma.GetXCoord();
	double material_normal_AlphaGamma_y = material_normal_AlphaGamma.GetYCoord();
	double material_normal_AlphaGamma_z = material_normal_AlphaGamma.GetZCoord();

    // TODO: implement and use a dot product for Point instead
	double cos_Theta_AlphaBetaGamma = material_normal_AlphaBeta_x * material_normal_AlphaGamma_x +
			                          material_normal_AlphaBeta_y * material_normal_AlphaGamma_y +
									  material_normal_AlphaBeta_z * material_normal_AlphaGamma_z;

	Point<dim> spatial_normal_AlphaBeta = spatial_delta_AlphaBeta / spatial_distance_AlphaBeta;
	double spatial_normal_AlphaBeta_x = spatial_normal_AlphaBeta.GetXCoord();
	double spatial_normal_AlphaBeta_y = spatial_normal_AlphaBeta.GetYCoord();
	double spatial_normal_AlphaBeta_z = spatial_normal_AlphaBeta.GetZCoord();

	Point<dim> spatial_normal_AlphaGamma = spatial_delta_AlphaGamma / spatial_distance_AlphaGamma;
	double spatial_normal_AlphaGamma_x = spatial_normal_AlphaGamma.GetXCoord();
	double spatial_normal_AlphaGamma_y = spatial_normal_AlphaGamma.GetYCoord();
	double spatial_normal_AlphaGamma_z = spatial_normal_AlphaGamma.GetZCoord();

    // TODO: implement and use a dot product for Point instead
	double cos_theta_AlphaBetaGamma = spatial_normal_AlphaBeta_x * spatial_normal_AlphaGamma_x +
			                          spatial_normal_AlphaBeta_y * spatial_normal_AlphaGamma_y +
									  spatial_normal_AlphaBeta_z * spatial_normal_AlphaGamma_z;

	double inv_cos_Theta_AlphaBetaGamma=1/cos_Theta_AlphaBetaGamma;
	double omega_AlphaBetaGamma=cos_theta_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double omega0_AlphaBetaGamma=cosine_teta0*inv_cos_Theta_AlphaBetaGamma;
	double epsilon0_AlphaBetaGamma=epsilon_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double eta0_AlphaBetaGamma=eta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma;
	double diff_omegas=omega_AlphaBetaGamma-omega0_AlphaBetaGamma;
	double force_coefficient=2*epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma*diff_omegas;
	double energy=epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma*diff_omegas*diff_omegas;
	SecondTensor <double> project_tensor_AlphaBeta (3,3,0);
	SecondTensor <double> project_tensor_AlphaGamma (3,3,0);
	double inv_distance_AlphaBeta = 1.0 / spatial_distance_AlphaBeta;
	double inv_distance_AlphaGamma = 1.0 / spatial_distance_AlphaGamma;

	Point<dim> term1(0.0, 0.0, 0.0);
	Point<dim> term2(0.0, 0.0, 0.0);
    /* TODO:
    Point<dim> term1 = inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * spatial_normal_AlphaBeta;
    Point<dim> term2 = inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * spatial_normal_AlphaGamma;
    */
	term1.SetXCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * spatial_normal_AlphaBeta_x);
	term1.SetYCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * spatial_normal_AlphaBeta_y);
	term1.SetZCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * spatial_normal_AlphaBeta_z);
	term2.SetXCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * spatial_normal_AlphaGamma_x);
	term2.SetYCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * spatial_normal_AlphaGamma_y);
	term2.SetZCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * spatial_normal_AlphaGamma_z);

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

    Point<dim> threeBody_DeformForce_Alpha = (term1 + term2) * force_coefficient ;
    Point<dim> threeBody_DeformForce_Beta = term1 * force_coefficient;
    Point<dim> threeBody_DeformForce_Gamma = term2 * force_coefficient;
    this->atoms->setForce(i, this->atoms->getForce(i) + threeBody_DeformForce_Alpha);
    this->atoms->setForce(j, this->atoms->getForce(j) - threeBody_DeformForce_Beta);
    this->atoms->setForce(k, this->atoms->getForce(k) - threeBody_DeformForce_Gamma);
	return energy;

}

template<int dim> double StillingerWeber<dim>::StillingerWeberTwoBody_Configurational(
    int i, int j, double sigma_AlphaBeta, double epsilon_AlphaBeta,
    double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
    double J_AlphaBeta, double a_AlphaBeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto bond_normal_vec = spatial_delta / spatial_distance;
	double material_stretch = material_distance / spatial_distance;
	double sigma_0=sigma_AlphaBeta/material_distance;
	double epsilon_0=epsilon_AlphaBeta/material_distance;
	double term1=pow(sigma_0*material_stretch, I_AlphaBeta);
	double term2=pow(sigma_0*material_stretch, J_AlphaBeta);
	double term3=1/(1-material_stretch*a_AlphaBeta*sigma_0);
	double term4=exp(sigma_0*material_stretch*term3);
	double energy = -epsilon_0 * A_AlphaBeta * material_stretch * (B_AlphaBeta * term1 - term2) * term4;

	//term1=pow(sigma_0*material_stretch, I_AlphaBeta)
	//term2=pow(sigma_0*material_stretch, J_AlphaBeta)
    //term3=1/(1-material_stretch*a_AlphaBeta*sigma_0);

	double inv_material_stretch = 1.0 / material_stretch;
	double coefficient0 = B_AlphaBeta * I_AlphaBeta;
	double coefficient1 = B_AlphaBeta * material_stretch;
	Point<dim> config_force = bond_normal_vec * energy * (inv_material_stretch + ((coefficient0 * term1 - J_AlphaBeta * term2) /
			                  (coefficient1 * term1 - material_stretch * term2)) + sigma_0 * term3 * term3);
    this->atoms->setForce(i, this->atoms->getForce(i) + config_force);
    this->atoms->setForce(j, this->atoms->getForce(j) - config_force);
	return energy;

}

template<int dim> double StillingerWeber<dim>::StillingerWeberThreeBody_Configurational(
    int i, int j, int k,
	double epsilon_AlphaBetaGamma, double eta_AlphaBetaGamma, double cosine_teta0) {

    auto spatial_delta_AlphaBeta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_delta_AlphaGamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spatial_distance_AlphaBeta = spatial_delta_AlphaBeta.norm();
    auto spatial_distance_AlphaGamma = spatial_delta_AlphaGamma.norm();
    auto material_delta_AlphaBeta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_delta_AlphaGamma = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(k);
    auto material_distance_AlphaBeta = material_delta_AlphaBeta.norm();
    auto material_distance_AlphaGamma = material_delta_AlphaGamma.norm();
    Point<dim> material_normal_AlphaBeta = material_delta_AlphaBeta / material_distance_AlphaBeta;
	double material_normal_AlphaBeta_x = material_normal_AlphaBeta.GetXCoord();
	double material_normal_AlphaBeta_y = material_normal_AlphaBeta.GetYCoord();
	double material_normal_AlphaBeta_z = material_normal_AlphaBeta.GetZCoord();

    Point<dim> material_normal_AlphaGamma = material_delta_AlphaGamma / material_distance_AlphaGamma;
	double material_normal_AlphaGamma_x = material_normal_AlphaGamma.GetXCoord();
	double material_normal_AlphaGamma_y = material_normal_AlphaGamma.GetYCoord();
	double material_normal_AlphaGamma_z = material_normal_AlphaGamma.GetZCoord();

    // TODO: implement and use a dot product for Point instead
	double cos_Theta_AlphaBetaGamma = material_normal_AlphaBeta_x * material_normal_AlphaGamma_x +
			                          material_normal_AlphaBeta_y * material_normal_AlphaGamma_y +
									  material_normal_AlphaBeta_z * material_normal_AlphaGamma_z;

	Point<dim> spatial_normal_AlphaBeta = spatial_delta_AlphaBeta / spatial_distance_AlphaBeta;
	double spatial_normal_AlphaBeta_x = spatial_normal_AlphaBeta.GetXCoord();
	double spatial_normal_AlphaBeta_y = spatial_normal_AlphaBeta.GetYCoord();
	double spatial_normal_AlphaBeta_z = spatial_normal_AlphaBeta.GetZCoord();

	Point<dim> spatial_normal_AlphaGamma = spatial_delta_AlphaGamma / spatial_distance_AlphaGamma;
	double spatial_normal_AlphaGamma_x = spatial_normal_AlphaGamma.GetXCoord();
	double spatial_normal_AlphaGamma_y = spatial_normal_AlphaGamma.GetYCoord();
	double spatial_normal_AlphaGamma_z = spatial_normal_AlphaGamma.GetZCoord();

    // TODO: implement and use a dot product for Point instead
	double cos_theta_AlphaBetaGamma = spatial_normal_AlphaBeta_x * spatial_normal_AlphaGamma_x +
			                          spatial_normal_AlphaBeta_y * spatial_normal_AlphaGamma_y +
									  spatial_normal_AlphaBeta_z * spatial_normal_AlphaGamma_z;
	double Omega_AlphaBetaGamma = cos_Theta_AlphaBetaGamma / cos_theta_AlphaBetaGamma;
	double inv_cos_Theta_AlphaBetaGamma = 1.0 / cos_Theta_AlphaBetaGamma;
	double epsilon0_AlphaBetaGamma=epsilon_AlphaBetaGamma*inv_cos_Theta_AlphaBetaGamma;
	double omega0_AlphaBetaGamma=cosine_teta0*inv_cos_Theta_AlphaBetaGamma;
	double eta0_AlphaBetaGamma=eta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma*cos_Theta_AlphaBetaGamma;
	double inv_Omega_AlphaBetaGamma=1/Omega_AlphaBetaGamma;
	double diff_omegas=inv_Omega_AlphaBetaGamma-omega0_AlphaBetaGamma;
	double coefficient=epsilon0_AlphaBetaGamma*eta0_AlphaBetaGamma;
	double force_coefficient = -2 * coefficient * inv_Omega_AlphaBetaGamma * diff_omegas + coefficient * diff_omegas * diff_omegas;
	std::cout << "force_coefficient: " << force_coefficient << std::endl;
	double energy = coefficient * Omega_AlphaBetaGamma * diff_omegas * diff_omegas;
//	SecondTensor <double> project_tensor_AlphaBeta (3,3,0);
//	SecondTensor <double> project_tensor_AlphaGamma (3,3,0);
	double inv_distance_AlphaBeta = 1.0 / material_distance_AlphaBeta;
	double inv_distance_AlphaGamma = 1.0 / material_distance_AlphaGamma;

	Point<dim> term1(0.0, 0.0, 0.0);
	Point<dim> term2(0.0, 0.0, 0.0);
    /* TODO:
    Point<dim> term1 = inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * material_normal_AlphaBeta;
    Point<dim> term2 = inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * material_normal_AlphaGamma;
    */
	term1.SetXCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * material_normal_AlphaBeta_x);
	term1.SetYCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * material_normal_AlphaBeta_y);
	term1.SetZCoord((inv_distance_AlphaGamma - (cos_theta_AlphaBetaGamma * inv_distance_AlphaBeta)) * material_normal_AlphaBeta_z);
	term2.SetXCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * material_normal_AlphaGamma_x);
	term2.SetYCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * material_normal_AlphaGamma_y);
	term2.SetZCoord((inv_distance_AlphaBeta - (cos_theta_AlphaBetaGamma * inv_distance_AlphaGamma)) * material_normal_AlphaGamma_z);

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

	Point<dim> threeBody_ConfigForce_Alpha = (term1 + term2) * force_coefficient;
	Point<dim> threeBody_ConfigForce_Beta = term1 * force_coefficient;
	Point<dim> threeBody_ConfigForce_Gamma = term2 * force_coefficient;
    this->atoms->setForce(i, this->atoms->getForce(i) + threeBody_ConfigForce_Alpha);
    this->atoms->setForce(j, this->atoms->getForce(j) - threeBody_ConfigForce_Beta);
    this->atoms->setForce(k, this->atoms->getForce(k) - threeBody_ConfigForce_Gamma);
	return energy;
}
