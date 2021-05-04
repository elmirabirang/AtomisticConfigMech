/*
 * energy.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#include <cmath>
#include <vector>
#include <iostream>
//---
#include "bond.h"
#include "energy.h"
#include "matrix.h"
//#include <nlopt.hpp>


using namespace std;

template<int dim> Energy <dim>::Energy(Atoms<dim> *atoms) {
    this->atoms = atoms;
    this->bond = new Bond<dim>(atoms);
}

template<int dim> Energy <dim>::Energy(Atoms<dim> *atoms, double sigma, double epsilon) {
    this->atoms = atoms;
	this->sigma = sigma;
	this->epsilon = epsilon;
    this->bond = new Bond<dim>(atoms);
}

template<int dim> Energy <dim>::~Energy() { delete bond; }
template<int dim> double Energy<dim>::InteratomicEnergy(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto distance = spatial_delta.norm();
//	double interatomic_distance=bond.SpatialBondDistance(atomi,atomj);
    double distanceInv = 1.0 / distance;
    double energy_ij = 4.0 * 1.0 * (pow(distanceInv, 6) - pow(distanceInv, 3));
//	try displacement locally.
    return energy_ij;
}

template<int dim> double Energy<dim>::ConfigInteratomicEnergy(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    double spatial_distance = spatial_delta.norm();
    double material_distance = material_delta.norm();
	double material_stretch = material_distance / spatial_distance;
	double sigma0 = 1.0 / material_distance;
	double epsilon0 = 1.0 / material_distance;
	double energy_ij = 4 * epsilon0 * material_stretch * pow(sigma0 * material_stretch, 6) * (pow(sigma0 * material_stretch, 6) - 1);
	return energy_ij;
}

template<int dim> double Energy<dim>::CutoffEnergy(int i) {
    double cutoff_energy = 0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        cutoff_energy += InteratomicEnergy(i, j);
    )

    return cutoff_energy;
}

//double Energy <dim>::ClusterEnergy(Atom atoma, double clusterR)
//
//		{
//	this -> atomi=atoma;
//	this -> cluster_radius=clusterR;
//
//	double cluster_energy=0;
//
//	vector < Atom <dim> > cluster_atoms=atomi.Cluster(cluster_radius);
//	vector <Atom <dim> >::iterator cluster_atom;
//
//	for (cluster_atom=cluster_atoms.begin();cluster_atom=cluster_atoms.end();++cluster_atom)
//	{
//		double energy_ij=InteratomicEnergy(atomi, cluster_atom);
//		cluster_energy+=energy_ij
//
//	}
//
//		}

template<int dim> double Energy<dim>::Material_Interatomic_Engergy(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    double spatial_distance = spatial_delta.norm();
    double material_distance = material_delta.norm();
	double interatomic_stretch = material_distance / spatial_distance;
	double energy_ij = (4 * epsilon * interatomic_stretch) *
                       (pow(sigma * interatomic_stretch, 12) - pow(sigma * interatomic_stretch, 6));
	return energy_ij;
}

template<int dim> double Energy<dim>::Spatial_Interatomic_Energy(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    double spatial_distance = spatial_delta.norm();
    double material_distance = material_delta.norm();
	double interatomic_stretch = spatial_distance / material_distance;
	double energy_ij = 4 * epsilon * (pow((sigma / interatomic_stretch), 12) - pow((sigma / interatomic_stretch), 6));
	return energy_ij;
}

template<int dim> double Energy<dim>::TotPotentialEnergy() {
	double potential_energy = 0.0;

    LOOP_OVER_ATOMS(this->atoms, i,
        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
			potential_energy += InteratomicEnergy(i, j);
		)
    )

	return potential_energy;
}

template<int dim> double Energy<dim>::atomTotalPotentialEnergy(int i) {
	double potential_energy = 0.0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        potential_energy += InteratomicEnergy(i, j);
    )

	return potential_energy;
}

template<int dim> double Energy<dim>::ConfigTotPotentialEnergy() {
	double potential_energy = 0.0;

    LOOP_OVER_ATOMS(this->atoms, i,
        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
            auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
            double interatomic_distance = spatial_delta.norm();
			double energy_ij = ConfigInteratomicEnergy(i, j);
			potential_energy += energy_ij * interatomic_distance;
		)
	)

	return potential_energy;
}

//template <int dim>
//double Energy <dim>::embeddingEnergy(vector < Atom <dim>* > atoms_list, double fitting_parameter, double atomic_number)
//{
//
//	typedef typename vector < Atom <dim>* >::iterator At;
//	double total_embedding_energy=0;
//	Bond <dim> bond;
//
//	for (At atom=atoms_list.begin(); atom!=atoms_list.end(); ++atom)
//	{
//    	Atom <dim> *Atm;
//    	Atm=*atom;
//
//		double electron_density=electronDensity(*Atm, atomic_number);
//
//		double embedding_energy=-fitting_parameter*sqrt(electron_density);
//
//		total_embedding_energy+=embedding_energy;
//
//	}
//
//	return (total_embedding_energy);
//
//}

template<int dim> double Energy<dim>::electronDensity(int i, double atomic_number) {
	double electron_density = 0.0;
	double total_electron_density = 0.0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        double spatial_distance = spatial_delta.norm();
		electron_density = exp(-2 * atomic_number * spatial_distance);
		total_electron_density += electron_density;
	)

	return total_electron_density;
}

template<int dim> double Energy<dim>::tersoffEnergy(double lambda, double A) {
    double total_tersoff_energy = 0.0;

    LOOP_OVER_ATOMS(this->atoms, i,
        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
            auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
            double spatial_distance = spatial_delta.norm();
			total_tersoff_energy += A * exp(-1 * lambda * spatial_distance);
        )
    )

    return total_tersoff_energy;
}


template<int dim> double Energy<dim>::morseFunction(double alpha1, double alpha2,
                                                    double r_01, double r_02, double E1, double E2, double delta,
                                                    double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                                    double s1, double s2, double s3) {
	double inv_h = 1.0 / h;
	double psi = 0.0;
	double potential_eng_tot = 0.0;

    LOOP_OVER_ATOMS(this->atoms, i,
        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
            auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
            double interatomic_distance = spatial_delta.norm();
			double M1 = exp(-1*2*alpha1*(interatomic_distance-r_01))-(2*exp(-1*alpha1*((interatomic_distance-r_01))));
			double M2 = exp(-1*2*alpha2*(interatomic_distance-r_02))-(2*exp(-1*alpha2*((interatomic_distance-r_02))));
			double control_repulsion = 0.0;

            if(interatomic_distance < r_cut) {
                double psi_param_fourth = pow(((interatomic_distance - r_cut) * inv_h), 4);
                double denominator = 1.0 / (1.0 + psi_param_fourth);
                psi = psi_param_fourth * denominator;
            } else {
                psi = 0.0;
            }

		    double part1 = 0.0;
		    double part2 = 0.0;
		    double part3 = 0.0;

		    if(r_s1 >= interatomic_distance) {
//		    	cout << "here" << endl;
		    	part1 = s1 * pow((r_s1 - interatomic_distance), 4);
		    }

		    if(r_s2 >= interatomic_distance) {
//		    	cout << "here" << endl;
		    	part2 = s2 * pow((r_s2 - interatomic_distance), 4);
		    }

		    if(r_s3 >= interatomic_distance) {
//		    	cout << "here" << endl;
		    	part3 = s3 * pow((r_s3 - interatomic_distance), 4);
		    }

            control_repulsion = part1 + part2 + part3;
            potential_eng_tot += (E1 * M1 + E2 * M2 + delta) * psi - control_repulsion;
        )
    )

	return potential_eng_tot;
}


template<int dim> double Energy<dim>::embeddingEnergy(double a, double beta1,
                                                      double beta2, double r_03, double r_04, double F0, double F2,
                                                      double q1, double q2, double q3, double q4, double Q1, double Q2,
                                                      double h, double r_cut) {
    double potential_eng_tot = 0.0;
    double inv_h = 1.0 / h;

    LOOP_OVER_ATOMS(this->atoms, i,
        double atom_electron_density = 0.0;
        double part3 = 0.0;
        double embedding_function = 0.0;

        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
            auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
            double interatomic_distance = spatial_delta.norm();
            double electron_density_ij = 0;
            double psi = 0.0;

            if(interatomic_distance < r_cut) {
                double psi_param_fourth = pow(((interatomic_distance - r_cut) * inv_h), 4);
                double denominator = 1.0 / (1.0 + psi_param_fourth);
                psi = psi_param_fourth * denominator;
            }

            electron_density_ij = (a * exp(-1 * beta1 * pow((interatomic_distance - r_03), 2)) +
                                       exp(-1 * beta2 * (interatomic_distance - r_04)) ) * psi;
            atom_electron_density += electron_density_ij;
        )

        if(atom_electron_density < 1) {
            part3 = q1 * pow((atom_electron_density - 1), 3) +
                    q2 * pow((atom_electron_density - 1), 4) +
                    q3 * pow((atom_electron_density - 1), 5) +
                    q4 * pow((atom_electron_density - 1), 6);

            embedding_function = F0 + (0.5 * F2 * pow((atom_electron_density - 1), 2)) + part3;
        }

        if(atom_electron_density > 1) {
            double denominator = 1.0 / (1 + Q2 * pow((atom_electron_density - 1), 3));
            embedding_function = (F0 + (0.5 * F2 * pow((atom_electron_density - 1), 2)) +
                                  q1 * pow((atom_electron_density - 1), 3) +
                                  Q1 * pow((atom_electron_density - 1), 4)) * denominator;
        }

        potential_eng_tot += embedding_function;
    )

    return potential_eng_tot;
}

template<int dim> double Energy<dim>::stillingerWeberTwoBody(
    int i, int j, double sigma_AlphaBeta, double epsilon_AlphaBeta,
    double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta, double q_AlphaBeta, double a_AlphaBeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    double distance = spatial_delta.norm();
    double distanceInv = sigma_AlphaBeta / distance;
    double energy_ij = 0.0;

    if(distance < a_AlphaBeta * sigma_AlphaBeta) {
        energy_ij = A_AlphaBeta * epsilon_AlphaBeta * (B_AlphaBeta * pow(distanceInv, p_AlphaBeta) - pow(distanceInv, q_AlphaBeta)) *
                    exp(sigma_AlphaBeta / (distance - (a_AlphaBeta * sigma_AlphaBeta)));
    }

    return energy_ij;
}

template<int dim> double Energy<dim>::stillingerWeberThreeBody(
    int j, int i, int k, double sigma_AlphaBeta,
    double sigma_AlphaGamma,double epsilon_AlphaBetaGamma, double gamma, double lambda_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

    auto spa_delta_alphabeta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spa_delta_alphagamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spa_delta_bettagamma = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(k);
	double distance_AlphaBeta = spa_delta_alphabeta.norm();
	double distance_AlphaGamma = spa_delta_alphagamma.norm();
	double distance_BetaGamma = spa_delta_bettagamma.norm();
	double energy_ij = 0.0;

    if(distance_AlphaBeta < (sigma_AlphaBeta * a_AlphaBeta) && distance_AlphaGamma < (sigma_AlphaGamma * a_AlphaGamma)) {
        Point<dim> normalVectorAlphaBeta = spa_delta_alphabeta / distance_AlphaBeta;
        Point<dim> normalVectorAlphaGamma = spa_delta_alphagamma / distance_AlphaGamma;
//		double cosine_teta_BetaAlphaGamma1=(normalVectorAlphaBeta.GetXCoord()*normalVectorAlphaGamma.GetXCoord()
//										  +normalVectorAlphaBeta.GetYCoord()*normalVectorAlphaGamma.GetYCoord()
//										  +normalVectorAlphaBeta.GetZCoord()*normalVectorAlphaGamma.GetZCoord());
//		cout << "" << << "\n";

        float cosine_teta_BetaAlphaGamma = (pow(distance_AlphaBeta, 2) + pow(distance_AlphaGamma, 2) - pow(distance_BetaGamma, 2)) /
                                           (2 * distance_AlphaBeta * distance_AlphaGamma);
        double term1 = (gamma * sigma_AlphaBeta) / (distance_AlphaBeta - (a_AlphaBeta * sigma_AlphaBeta));
        double term2 = (gamma * sigma_AlphaGamma) / (distance_AlphaGamma - (a_AlphaGamma * sigma_AlphaGamma));
        double terms_sum = term1 + term2;
        double difference_cosines = cosine_teta_BetaAlphaGamma - cosine_teta0;
        energy_ij = epsilon_AlphaBetaGamma * lambda_AlphaBetaGamma * exp(terms_sum) * difference_cosines * difference_cosines;
    }

	return energy_ij;
}

template<int dim> double Energy<dim>::totalStillingerWeberEnergy(
    int i, double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
    double Gamma, double lambda_AlphaBetaGamma, double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma,
    double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta, double q_AlphaBeta) {

	double total_2BodyEng = 0.0;
	double total_3BodyEng = 0.0;

    LOOP_OVER_ATOM_HALF_NEIGHBORS(this->atoms, i, j,
		total_2BodyEng += stillingerWeberTwoBody(
            i, j, sigma_AlphaBeta, epsilon_AlphaBetaGamma,
            A_AlphaBeta, B_AlphaBeta, p_AlphaBeta, q_AlphaBeta, a_AlphaBeta);
    )

    LOOP_OVER_ATOM_NEIGHBORS_3BODY(atoms, i, j, k,
	    total_3BodyEng += stillingerWeberThreeBody(
            j, i, k, sigma_AlphaBeta, sigma_AlphaGamma,
            epsilon_AlphaBetaGamma, Gamma, lambda_AlphaBetaGamma, cosine_teta0, a_AlphaBeta, a_AlphaGamma);
	)

	return total_2BodyEng + total_3BodyEng;
}

template<int dim> double Energy<dim>::deformationalSWtwoBody(
    int i, int j, double sigma_alphabeta, double epsilon_alphabeta,
    double A_alphabeta, double B_alphabeta, double I_alphabeta, double J_alphabeta, double a_alphabeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    double spatial_distance = spatial_delta.norm();
    double material_bond_distance = material_delta.norm();
	double spatial_stretch = spatial_distance / material_bond_distance;
	double epsilon_0 = epsilon_alphabeta / material_bond_distance;
	double sigma_0 = sigma_alphabeta / material_bond_distance;
	double inv_stretch = sigma_0 / spatial_stretch;
	double energy = 0.0;

	if(spatial_stretch < a_alphabeta * sigma_0) {
		energy = epsilon_0 * A_alphabeta * (B_alphabeta * pow(inv_stretch, I_alphabeta) -
                 pow(inv_stretch, J_alphabeta)) * exp(sigma_0 / (spatial_stretch - (a_alphabeta * sigma_0)));
	}

	return energy;
}

template<int dim> double Energy<dim>::deformationalSWthreeBody(
    int i, int j, int k, double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
	double mu, double eta_AlphaBetaGamma ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

    auto spa_delta_alphabeta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spa_delta_alphagamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto mat_delta_alphabeta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto mat_delta_alphagamma = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(k);
    auto mat_delta_betagamma = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(k);
    double spa_dist_alphabeta = spa_delta_alphabeta.norm();
    double spa_dist_alphagamma = spa_delta_alphagamma.norm();
	double mat_dist_alphabeta = mat_delta_alphabeta.norm();
	double mat_dist_alphagamma = mat_delta_alphagamma.norm();
	double mat_dist_betagamma = mat_delta_betagamma.norm();
	double M = (mat_dist_alphabeta + mat_dist_alphagamma + mat_dist_betagamma) / 2;
	double area = sqrt(M * abs(M - mat_dist_alphabeta) * abs(M - mat_dist_alphagamma) * abs(M - mat_dist_betagamma));
	std::cout << area << endl;

	double epslion_0_AlphaBetaGamma = epsilon_AlphaBetaGamma / area;
	double sigma_0_AlphaBeta = sigma_AlphaBeta / mat_dist_alphabeta;
	double sigma_0_AlphaGamma = sigma_AlphaGamma / mat_dist_alphagamma;
	double spatial_stretch_AlphaBeta = spa_dist_alphabeta / mat_dist_alphabeta;
	double spatial_stretch_AlphaGamma = spa_dist_alphagamma / mat_dist_alphagamma;

	Point<dim> spatial_normal_AlphaBeta = spa_delta_alphabeta / spa_dist_alphabeta;
	Point<dim> spatial_normal_AlphaGamma = spa_delta_alphagamma / spa_dist_alphagamma;
	Point<dim> material_normal_AlphaBeta = mat_delta_alphabeta / mat_dist_alphabeta;
	Point<dim> material_normal_AlphaGamma = mat_delta_alphagamma / mat_dist_alphagamma;
	SecondTensor<double> rotation_mat_AlphaBeta = this->bond->RotationMatrix(spatial_normal_AlphaBeta, material_normal_AlphaBeta);
	SecondTensor<double> rotation_mat_AlphaGamma = this->bond->RotationMatrix(spatial_normal_AlphaGamma, material_normal_AlphaGamma);

	std::vector<double> material_normalvec_AlphaBeta;
	material_normalvec_AlphaBeta.push_back(material_normal_AlphaBeta.GetXCoord());
	material_normalvec_AlphaBeta.push_back(material_normal_AlphaBeta.GetYCoord());
	material_normalvec_AlphaBeta.push_back(material_normal_AlphaBeta.GetZCoord());

	std::vector<double> material_normalvec_AlphaGamma;
	material_normalvec_AlphaGamma.push_back(material_normal_AlphaGamma.GetXCoord());
	material_normalvec_AlphaGamma.push_back(material_normal_AlphaGamma.GetYCoord());
	material_normalvec_AlphaGamma.push_back(material_normal_AlphaGamma.GetZCoord());

	double energy = 0.0;

	std::vector<double> spatial_bondvec_AlphaBeta = rotation_mat_AlphaBeta * material_normalvec_AlphaBeta;
	std::vector<double> spatial_bondvec_AlphaGamma = rotation_mat_AlphaGamma * material_normalvec_AlphaGamma;

	double term1 = (spatial_bondvec_AlphaBeta[0] * spatial_bondvec_AlphaGamma[0] +
                    spatial_bondvec_AlphaBeta[1] * spatial_bondvec_AlphaGamma[1] +
                    spatial_bondvec_AlphaBeta[2] * spatial_bondvec_AlphaGamma[2]) - cosine_teta0;
	double term1_squared = pow(term1, 2);

	cout << "term1_squared: " << term1_squared << endl;

	if(spatial_stretch_AlphaBeta < a_AlphaBeta * sigma_0_AlphaBeta && spatial_stretch_AlphaGamma < a_AlphaGamma * sigma_0_AlphaGamma) {
	       energy=epslion_0_AlphaBetaGamma*eta_AlphaBetaGamma*term1_squared*
		   exp(mu*sigma_0_AlphaBeta/(spatial_stretch_AlphaBeta-a_AlphaBeta*sigma_0_AlphaBeta))
	       *exp(mu*sigma_0_AlphaGamma/(spatial_stretch_AlphaGamma-a_AlphaGamma*sigma_0_AlphaGamma));
	}

	return energy;
}


template<int dim> double Energy<dim>::configurationalSWtwoBody(
    int i, int j, double sigma_alphabeta, double epsilon_alphabeta,
    double A_alphabeta, double B_alphabeta, double I_alphabeta, double J_alphabeta, double a_alphabeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    double spatial_distance = spatial_delta.norm();
    double mat_dist_alphabeta = material_delta.norm();
	double material_stretch_AlphaBeta = mat_dist_alphabeta / spatial_distance;
	double sigma_0_AlphaBeta = sigma_alphabeta/mat_dist_alphabeta;
	double epsilon_0_AlphaBeta = epsilon_alphabeta/mat_dist_alphabeta;
	double term1 = B_alphabeta * (pow((sigma_0_AlphaBeta * material_stretch_AlphaBeta), I_alphabeta));
	double term2 = pow((sigma_0_AlphaBeta * material_stretch_AlphaBeta), J_alphabeta);
	double term3 = term1 - term2;
	double term4 = material_stretch_AlphaBeta * sigma_0_AlphaBeta /
                   (1 - (a_alphabeta * material_stretch_AlphaBeta * sigma_0_AlphaBeta));

	if(material_stretch_AlphaBeta < a_alphabeta * sigma_0_AlphaBeta) {
		return epsilon_0_AlphaBeta * A_alphabeta * material_stretch_AlphaBeta * term3 * exp(term4);
	}

	return 0.0;
}

template<int dim> double Energy<dim>::configurationalSWthreeBody(
    int i, int j, int k, double sigma_AlphaBeta, double sigma_AlphaGamma, double epsilon_AlphaBetaGamma,
	double mu, double eta_AlphaBetaGamma, double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	double configThreeBodyEng = 0.0;
    auto spatialbond_AlphaBeta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatialbond_AlphaGamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spatialbond_BetaGamma = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(k);
    auto materialbond_AlphaBeta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto materialbond_AlphaGamma = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(k);
    auto materialbond_BetaGamma = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(k);
    double spa_dist_alphabeta = spatialbond_AlphaBeta.norm();
    double spa_dist_alphagamma = spatialbond_AlphaGamma.norm();
    double spa_dist_betagamma = spatialbond_BetaGamma.norm();
	double mat_dist_alphabeta = materialbond_AlphaBeta.norm();
	double mat_dist_alphagamma = materialbond_AlphaGamma.norm();
	double mat_dist_betagamma = materialbond_BetaGamma.norm();
	double material_stretch_AlphaBeta = mat_dist_alphabeta / spa_dist_alphabeta;
	double material_stretch_AlphaGamma = mat_dist_alphagamma / spa_dist_alphagamma;

	double m=(spa_dist_alphabeta+spa_dist_alphagamma+spa_dist_betagamma)/2;
	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;

	double area=sqrt(m*(m-spa_dist_alphabeta)*(m-spa_dist_alphagamma)*(m-spa_dist_betagamma));
	double Area=sqrt(M*(M-mat_dist_alphabeta)*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma));

	double ratio_areas=Area/area;

	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaBeta/mat_dist_alphagamma;

	double epsilon_0_AlphaBetaGamma=epsilon_AlphaBetaGamma/Area;

	std::vector <double> spa_bondvec_AlphaBeta;
	spa_bondvec_AlphaBeta.push_back(spatialbond_AlphaBeta.GetXCoord());
	spa_bondvec_AlphaBeta.push_back(spatialbond_AlphaBeta.GetYCoord());
	spa_bondvec_AlphaBeta.push_back(spatialbond_AlphaBeta.GetZCoord());

	std::vector <double> spa_bondvec_AlphaGamma;
	spa_bondvec_AlphaGamma.push_back(spatialbond_AlphaGamma.GetXCoord());
	spa_bondvec_AlphaGamma.push_back(spatialbond_AlphaGamma.GetYCoord());
	spa_bondvec_AlphaGamma.push_back(spatialbond_AlphaGamma.GetZCoord());

    auto spa_normal_alphabeta = spatialbond_AlphaBeta / spa_dist_alphabeta;
    auto spa_normal_alphagamma = spatialbond_AlphaGamma / spa_dist_alphagamma;
    auto mat_normal_alphabeta = materialbond_AlphaBeta / mat_dist_alphabeta;
    auto mat_normal_alphagamma = materialbond_AlphaGamma / mat_dist_alphagamma;
	SecondTensor<double> RotationMat_AlphaBeta = this->bond->RotationMatrix(spa_normal_alphabeta, mat_normal_alphabeta);
	SecondTensor<double> RotationMat_AlphaGamma = this->bond->RotationMatrix(spa_normal_alphagamma, mat_normal_alphagamma);

	SecondTensor<double> RotationMat_AlphaBetaTranspose=RotationMat_AlphaBeta.transpose();
	SecondTensor<double> RotationMat_AlphaGammaTranspose=RotationMat_AlphaGamma.transpose();

	std::vector <double> mat_bondvec_AlphaBeta=RotationMat_AlphaBetaTranspose*spa_bondvec_AlphaBeta;
	std::vector <double> mat_bondvec_AlphaGamma=RotationMat_AlphaGammaTranspose*spa_bondvec_AlphaGamma;

	double magnitude_mat_bondvec_AlphaBeta=sqrt(pow(mat_bondvec_AlphaBeta[0],2)+
			                                    pow(mat_bondvec_AlphaBeta[1],2)+
												pow(mat_bondvec_AlphaBeta[2],2));

	double magnitude_mat_bondvec_AlphaGamma=sqrt(pow(mat_bondvec_AlphaGamma[0],2)+
			                                     pow(mat_bondvec_AlphaGamma[1],2)+
												 pow(mat_bondvec_AlphaGamma[2],2));

	double term1=(mat_bondvec_AlphaBeta[0]*mat_bondvec_AlphaGamma[0]+
			      mat_bondvec_AlphaBeta[1]*mat_bondvec_AlphaGamma[1]+
				  mat_bondvec_AlphaBeta[2]*mat_bondvec_AlphaGamma[2])/
				 (magnitude_mat_bondvec_AlphaBeta*magnitude_mat_bondvec_AlphaGamma);

	double diff_cosines=term1-cosine_teta0;

//	if (diff_cosines < 1e-2)
//	{
//
//		diff_cosines=0.;
//
//	}

	double term2=pow((diff_cosines),2);
	double term3=material_stretch_AlphaBeta*mu*sigma_0_AlphaBeta
			     /(1-(material_stretch_AlphaBeta*a_AlphaBeta*sigma_0_AlphaBeta));
	double term4=material_stretch_AlphaGamma*mu*sigma_0_AlphaGamma
			     /(1-(material_stretch_AlphaGamma*a_AlphaGamma*sigma_0_AlphaGamma));

	if(material_stretch_AlphaBeta < a_AlphaBeta*sigma_0_AlphaBeta && material_stretch_AlphaGamma < a_AlphaGamma*sigma_0_AlphaGamma) {
	    configThreeBodyEng = epsilon_0_AlphaBetaGamma * eta_AlphaBetaGamma * ratio_areas * term2 * exp(term3) * exp(term4);
	}

	return configThreeBodyEng;
}

template<int dim> double Energy<dim>::LennrdJonesPotential(double sigma, double epsilon, double distance_ab) {
	double term = pow((sigma / distance_ab), 6);
	double energy= 4 * epsilon * (pow(term, 2) - 2 * term);
	return energy;
}

template<int dim> double Energy<dim>::phaseAverageLJ(
    double sigma, double epsilon,
    double boltzman_const, double plancks_const, int num_quad_points) {

    LOOP_OVER_ATOMS(this->atoms, i,
		double atoma_temperature = this->atoms->getTemperature(i);
		double sigma_a = sqrt(boltzman_const * atoma_temperature);
		double frequency_a = this->atoms->getFrequency(i);
		double coefficient_a = sqrt(2) * (sigma_a / frequency_a);
		double LJ_potential = 0.0;
		Point<dim> atoma_mean_position = this->atoms->getSpatialMeanPosition(i);

        LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
			double atomb_temperature = this->atoms->getTemperature(j);
			double sigma_b = sqrt(boltzman_const * atomb_temperature);
			double frequency_b = this->atoms->getFrequency(j);
			double coefficient_b = sqrt(2) * (sigma_b / frequency_b);
		    Point<dim> atomb_mean_position = this->atoms->getSpatialMeanPosition(j);

			for(int n = 1; n < num_quad_points; ++n) {
				double weight=pow(M_PI,(n/2))/(2*n);
				double QP=sqrt(n/2);
				Point <dim> quadrature_pa;
				Point <dim> quadrature_pb;

				if (n<4) {
					if (n==1) {
						quadrature_pa.SetXCoord(coefficient_a*QP);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						Point <dim> dist_vec;
						double distance=0.;
						double weighted_energy=0.;

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;

						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(-coefficient_a*QP);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						dist_vec.SetXCoord(0.);
						dist_vec.SetYCoord(0.);
						dist_vec.SetZCoord(0.);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;

					}

					if (n==2) {
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(coefficient_a*QP);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						Point <dim> dist_vec(0., 0., 0.);
						double distance=0.;
						double weighted_energy=0.;

						dist_vec=quadrature_pa-quadrature_pb
								 +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;

						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(-coefficient_a*QP);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;

					}

					if (n==3) {
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(coefficient_a*QP);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						Point <dim> dist_vec(0., 0., 0.);
						double distance=0.;
						double weighted_energy=0.;

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;

						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(-coefficient_a*QP);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;
					}
				}

				if (n>4) {
					if (n==4) {
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(coefficient_b*QP);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						Point <dim> dist_vec(0., 0., 0.);
						double distance=0.;
						double weighted_energy=0.;

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(-coefficient_b*QP);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;

					}

					if (n==5) {
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.0);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(-coefficient_b*QP);
						quadrature_pb.SetXCoord(0.);

						Point <dim> dist_vec(0., 0., 0.);
						double distance=0.;
						double weighted_energy=0.;

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;

						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(-coefficient_b*QP);
						quadrature_pb.SetXCoord(0.);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;

						LJ_potential=LJ_potential+weighted_energy;

					}

					if (n==6) {
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(coefficient_b*QP);

						Point <dim> dist_vec(0., 0., 0.);
						double distance=0;
						double weighted_energy=0;

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;

						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;

						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);
						quadrature_pa.SetXCoord(0.);

						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(0.);
						quadrature_pb.SetXCoord(-coefficient_b*QP);

						dist_vec=quadrature_pa-quadrature_pb
								             +atoma_mean_position-atomb_mean_position;



						distance=dist_vec.norm();
						weighted_energy=LennrdJonesPotential(sigma, epsilon, distance)*weight;
						LJ_potential=LJ_potential+weighted_energy;
					}
				}
			}
		)
	)
}

template<int dim> double Energy<dim>::internalEnergyLJ(
    double temperature, double boltzman_const, double planck_const,
    double sigma, double epsilon) {

	int num_atoms = this->atoms->getNumberOfAtoms();
	double internal_energy_term1 = 0.0;

    LOOP_OVER_ATOMS(this->atoms, i,
		double frequency = this->atoms->getFrequency(i);
		double entropy = this->atoms->getEntropy(i);

		this->atoms->setTemperature(i, temperature);

		double internal_energy_atoma_term1=planck_const*frequency*
										   exp((entropy/(3*boltzman_const))-(1.333333333)
												+(0.33333333)*log(num_atoms));

		internal_energy_term1=internal_energy_term1+1.5*internal_energy_atoma_term1;
		double coefficient=pow((1/sqrt(M_PI)),6);
	)
}

template class Energy<2>;
template class Energy<3>;
