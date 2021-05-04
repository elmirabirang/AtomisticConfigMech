/*
 * force.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: Elmira Birang
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
//---
#include "atom.h"
#include "bond.h"
#include "energy.h"
#include "force.h"
#include "point.h"

using namespace std;

template<int dim> Force <dim>::Force(Atoms<dim> *atoms, double sigma, double epsilon) {
    this->atoms = atoms;
    this->sigma = sigma;
    this->epsilon = epsilon;
    this->bond = new Bond<dim>(atoms);
}

template<int dim> Force<dim>::~Force() {}

// TODO: optimize these calculations
template<int dim> Point<dim> Force<dim>::BondForce(int i, int j) {
	// this function calculates interaction force between
	// a pair of atoms using Lennard-Jones potential energy.
    auto delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto interatomic_distance = delta.norm();
    auto interatomic_normal = delta / interatomic_distance;
    double distanceInv = 1.0 / interatomic_distance;
    return interatomic_normal * (-4 * distanceInv) * (12 * pow(distanceInv, 12) - 6 * pow(distanceInv, 6));
}

// TODO: optimize these calculations
// initial distances can be calculated only once
// material distances too?
template<int dim> Point<dim> Force<dim>::MaterialBondForce(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto interatomic_normal = material_delta / material_distance;
    auto initial_delta = this->atoms->getInitialPosition(i) - this->atoms->getInitialPosition(j);
    auto initial_distance = initial_delta.norm();
    auto interatomic_stretch = material_distance / spatial_distance;
    auto interatomic_stretch_0 = material_distance / initial_distance;
	double sigma0 = 1.0 / (material_distance * initial_distance);
	double epsilon0 = 1.0 / (material_distance * initial_distance);
	return interatomic_normal * 4 * epsilon0 * (13 * pow((sigma0 * interatomic_stretch * interatomic_stretch_0), 12)
	                                            -7 * pow((sigma0 * interatomic_stretch * interatomic_stretch_0), 6));
}

// TODO: I would mix these forces calculation when possible because there are many terms that can be reused!
template<int dim> double Force<dim>::interatomicConfigForceValue(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto interatomic_normal = material_delta / material_distance;
    auto initial_delta = this->atoms->getInitialPosition(i) - this->atoms->getInitialPosition(j);
    auto initial_distance = initial_delta.norm();
    auto interatomic_stretch = material_distance / spatial_distance;
    auto interatomic_stretch_0 = material_distance / initial_distance;
	double sigma0 = 1.0 / (material_distance * initial_distance);
	double epsilon0 = 1.0 / (material_distance * initial_distance);
	double term = (sigma0 * interatomic_stretch * interatomic_stretch_0);
	return 4 * epsilon0 * (13 * pow(term, 12) -7 * pow(term, 6));
}

template<int dim> Point<dim> Force<dim>::CriticalMaterialBondForce(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto interatomic_normal = material_delta / material_distance;
    auto initial_delta = this->atoms->getInitialPosition(i) - this->atoms->getInitialPosition(j);
    auto initial_distance = initial_delta.norm();
    auto interatomic_stretch = material_distance / spatial_distance;
    auto interatomic_stretch_0 = material_distance / initial_distance;
	double stretch_critical = 98.0 / 26.0;
	double sigma0 = 1.0 / (material_distance * initial_distance);
	double epsilon0 = 1.0 / (material_distance * initial_distance);

//	Point <dim> force_ij= interatomic_normal*4*epsilon0*
//			              (13*pow((sigma0*(stretch_critical)*interatomic_stretch_0),12)
//	                      -7*pow((sigma0*(stretch_critical)*interatomic_stretch_0),6));

	return interatomic_normal * stretch_critical * epsilon0;
}

template<int dim> Point<dim> Force<dim>::SpatialBondForce(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto interatomic_normal = material_delta / material_distance;
    auto interatomic_stretch = material_distance / spatial_distance;
	return interatomic_normal * (4 * (this->epsilon / interatomic_stretch)) *
			                    (12 * pow((this->sigma / interatomic_stretch), 12) - 6 * pow((this->sigma / interatomic_stretch), 6));
}

template<int dim> Point<dim> Force<dim>::ConfigurationalForce(int i) {
	Point<dim> config_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		config_force += MaterialBondForce(i, j);
	)

	return config_force;
}

template<int dim> double Force<dim>::EnergyRelease(int i, double delta_x) {
	double config_force = 0.0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		auto config_interaction_force = MaterialBondForce(i, j);
		config_force += config_interaction_force.GetXCoord() * delta_x;
    )

	return config_force;
}

//calculate the interaction force for each atom with neighbors.
//The neighbor list is in a way that atoms only interact with atoms that have larger id number.
//example: atom0 interacts with atom1,atom2,.. , atom1 interacts with atom2 and atom3

template<int dim> Point<dim> Force<dim>::ResultantForce(int i) {
    Point<dim> resultant_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
         resultant_force += BondForce(i, j);
    )

    return resultant_force;
}

template<int dim> Point<dim> Force<dim>::ResultantConfigForce(int i) {
	Point<dim> config_force(0.0, 0.0, 0.0);

    // TODO: It may be interesting to calculate all these stages together to reuse computations
    // I am not sure if the compiler can do it (maybe with -ipo?)
    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        auto spatial_distance = spatial_delta.norm();
        auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
        auto material_distance = material_delta.norm();
        double sigma0 = 1.0 / material_distance;
        double epsilon0 = 1.0 / material_distance;
		double deriv_interatomic_config_eng = DerivMaterialEnergy(i, j);
		double material_bond_stretch = material_distance / spatial_distance;
		Point<dim> spatial_normal = material_delta / material_distance;
	    Energy<dim> energy(this->atoms, sigma0, epsilon0);
		double material_energy = energy.ConfigInteratomicEnergy(i, j);
		config_force += spatial_normal * (deriv_interatomic_config_eng * material_bond_stretch + material_energy) * -1;
    )

	return config_force;
}

template<int dim> double Force<dim>::DerivMaterialEnergy(int i, int j) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
	double sigma0 = 1.0 / material_distance;
	double epsilon0 = 1.0 / material_distance;
	double interatomic_stretch = material_distance / spatial_distance;
	return 4 * epsilon0 * (13 * pow((sigma0 * interatomic_stretch), 12) - 7 * pow((sigma0 * interatomic_stretch), 6));
}

double repulsionControlDerivative(double r_sn, double interatomic_distance) {
    return -1 * 4 * pow((r_sn - interatomic_distance), 3);
}

double Psi(double r_cut, double h, double interatomic_distance) {
	double psi = 0.0;

	if(interatomic_distance < r_cut) {
        double inv_h = 1.0 / h;
        double psi_param_fourth = pow(((interatomic_distance - r_cut) * inv_h), 4);
        double denominator = 1.0 / (1 + psi_param_fourth);
        psi = psi_param_fourth * denominator;
	}

    return psi;
}

double psiDerivative(double r_cut, double h, double interatomic_distance) {
//	double inv_h=1.0/h;
//	double psi_param_third=pow( ((interatomic_distance-r_cut)*inv_h),3 );
//	double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );
//
//	double denominator=1.0/pow((1+ psi_param_fourth),2);
//
//	double psi_deriv= (4*psi_param_third)*denominator;
	double psi_deriv = 0.0;

	if(interatomic_distance < r_cut) {
        double nominator = 4 * pow(h, 4) * pow((interatomic_distance - r_cut), 3);
        double denominator = 1.0 / (pow((pow((interatomic_distance - r_cut), 4) + pow(h, 4)), 2));
        psi_deriv = nominator * denominator;
	}

	return psi_deriv;
}

double calculateM(double r_0, double interatomic_distance, double alpha) {
	return exp(-1 * 2 * alpha * (interatomic_distance - r_0)) - (2 * exp(-1 * alpha * ((interatomic_distance - r_0))));
}

double derivM(double alpha, double r_0, double interatomic_distance) {
	return -1 * 2 * alpha * exp(-1 * 2 * alpha * (interatomic_distance - r_0)) + 2 *
                    alpha * exp(-1 * alpha * ((interatomic_distance - r_0))) ;
}

template<int dim> Point<dim> Force<dim>::derivativeMorseFunction(int i, double alpha1, double alpha2,
                                                                 double r_01, double r_02, double E1, double E2, double delta,
                                                                 double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                                                 double s1, double s2, double s3) {
	Point <dim> atoma_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        double interatomic_bond = spatial_delta.norm();
        auto bond_direction = spatial_delta / interatomic_bond;
		double part1 = 0.0;
		double part2_1 = 0;
		double part2_2 = 0.0;
		double part2_3 = 0.0;
		double morse_deriv_ij = 0.0;
		Point<dim> morse_force_ij(0.0, 0.0, 0.0);

	    if(interatomic_bond < r_cut) {
		     part1 = (E1 * derivM(alpha1, r_01, interatomic_bond) +
                      E2 * derivM(alpha2, r_02, interatomic_bond)) * Psi(r_cut, h, interatomic_bond) +
		    	     (E1 * calculateM(r_01, interatomic_bond, alpha1) +
                      E2 * calculateM(r_02, interatomic_bond, alpha2) + delta) * psiDerivative(r_cut, h, interatomic_bond);
	    }

	    if(r_s1 > interatomic_bond) {
//	    	cout << "here" << endl;
	    	part2_1= s1*repulsionControlDerivative (r_s1, interatomic_bond);
	    }

	    if(r_s2 > interatomic_bond) {
//	    	cout << "here" << endl;
	    	part2_2=s2*repulsionControlDerivative (r_s2, interatomic_bond);
	    }

	    if(r_s3 > interatomic_bond) {
//	    	cout << "here" << endl;
	    	part2_3=s3*repulsionControlDerivative (r_s3, interatomic_bond);
	    }

	    morse_deriv_ij = part1 - part2_1 - part2_2 - part2_3;
	    morse_force_ij = bond_direction * morse_deriv_ij;
	    atoma_force = atoma_force + morse_force_ij;
	)

	return atoma_force;
}


template<int dim> double Force<dim>::electronDensity(int i, double a, double beta1,
                                                     double beta2, double r_03, double r_04,
                                                     double h, double r_cut) {
	Point<dim> atoma_force;
	double atom_electron_density=0;
	double inv_h=1.0/h;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        auto interatomic_distance = delta.norm();
		double psi = 0.0;

	    if(interatomic_distance < r_cut) {
	       double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );
	       double denominator=1.0/(1+psi_param_fourth);
	       psi=psi_param_fourth*denominator;
	    }

		double electron_density_ij = (a * exp (-1*beta1 * pow((interatomic_distance-r_03),2)) + exp (-1*beta2 *(interatomic_distance-r_04)) ) * psi;
		atom_electron_density += electron_density_ij;
	)

	return atom_electron_density;
}

double derivativeElectronDensity(double interatomic_distance, double a, double beta1,
                                 double beta2, double r_03, double r_04,
                                 double h, double r_cut) {
	double deriv_electron_density=((-2*a*beta1*(interatomic_distance-r_03))*exp(-1*beta1*pow ((interatomic_distance-r_03),2))
			                      -beta2*exp(-1*beta2*(interatomic_distance-r_04)) )*Psi(r_cut, h, interatomic_distance) +
			                      (a*exp(-1*beta1*pow ((interatomic_distance-r_03),2))+exp(-1*beta2*(interatomic_distance-r_04)))*psiDerivative(r_cut, h, interatomic_distance);
	return deriv_electron_density;
}

//F'(rho)
double derivativeEmbeddingEnergy(double electron_density,
                                 double F0, double F2,
                                 double q1, double q2,
								 double q3, double q4,
								 double Q1, double Q2,
                                 double h, double r_cut) {
	double derivative_embedding_energy=0.0;

	if (electron_density < 1) {
        derivative_embedding_energy=F2*(electron_density-1)
                                    + 3*q1*pow ((electron_density-1),2)
                                    + 4*q2*pow ((electron_density-1),3)
                                    + 5*q3*pow ((electron_density-1),4)
                                    + 6*q4*pow ((electron_density-1),5);
	}

	if (electron_density > 1) {
        double denominator= 1.0/(pow((1+Q2*pow((electron_density-1),3)),2));
        double nominator= (F2*(electron_density-1)+3*q1*pow((electron_density-1),2)+4*Q1*pow((electron_density-1),3)) * (1+Q2*pow((electron_density-1),3))
        	              - (F0+0.5*F2*pow((electron_density-1),2)+q1*pow((electron_density-1),3)+Q1*pow((electron_density-1),4))*(3*Q2*pow((electron_density-1),2));
		derivative_embedding_energy=nominator*denominator;
	}

	return derivative_embedding_energy;
}

template<int dim> Point<dim> Force<dim>::embeddingForce(int i, double a, double beta1,
                                                        double beta2, double r_03, double r_04, double F0, double F2,
                                                        double q1, double q2, double q3, double q4, double Q1, double Q2,
		                                                double h, double r_cut) {
	Point<dim> atoma_force(0.0, 0.0, 0.0);
	double atom_electron_density = this->electronDensity(i, a, beta1, beta2, r_03, r_04, h, r_cut);
	double atom_deriv_embedding_eng=derivativeEmbeddingEnergy(atom_electron_density,
                                                     F0, F2,
                                                     q1, q2,
			                                         q3, q4,
			                                         Q1, Q2,
                                                     h, r_cut);
	double inv_h=1.0/h;
	double psi=0.0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		Point<dim> Embedding_Force(0.0, 0.0, 0.0);
        auto delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
		auto interatomic_distance = delta.norm();
		auto bond_direction = delta / interatomic_distance;
		double neighbor_electron_density = this->electronDensity(j, a, beta1, beta2, r_03, r_04, h, r_cut);
		double neighbor_deriv_embedding_eng = derivativeEmbeddingEnergy(neighbor_electron_density, F0, F2, q1, q2, q3, q4, Q1, Q2, h, r_cut);
		double deriv_electron_density = derivativeElectronDensity(interatomic_distance, a, beta1, beta2, r_03, r_04, h, r_cut);
		Embedding_Force = bond_direction * (neighbor_deriv_embedding_eng + atom_deriv_embedding_eng) * deriv_electron_density;
		atoma_force = atoma_force + Embedding_Force;
	)

	return atoma_force;

}

template<int dim> Point<dim> Force<dim>::configForceEamPair(int i, double alpha1, double alpha2,
                                                            double r_01, double r_02, double E1, double E2, double delta,
                                                            double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                                            double s1, double s2, double s3) {
	Point<dim> atoma_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        auto interatomic_bond = spatial_delta.norm();
        auto bond_direction = spatial_delta / interatomic_bond;
        auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
        auto material_bond_distance = material_delta.norm();
        auto material_bond_stretch = material_bond_distance / interatomic_bond;
		double inv_material_bond_stretch = 1.0 / material_bond_stretch;

		double alpha1_mod=alpha1*material_bond_distance;
		double alpha2_mod=material_bond_distance*alpha2;
		double r_01_mod=1.0/material_bond_distance*r_01;
		double r_02_mod=1.0/material_bond_distance*r_02;
		double r_cut_mod=1.0/material_bond_distance*r_cut;
		double h_mod=1.0/material_bond_distance*h;
		double r_s1_mod=1.0/material_bond_distance*r_s1;
		double r_s2_mod=1.0/material_bond_distance*r_s2;
		double r_s3_mod=1.0/material_bond_distance*r_s3;
		double s1_mod=pow(material_bond_distance,4)*s1;
		double s2_mod=pow(material_bond_distance,4)*s2;
		double s3_mod=pow(material_bond_distance,4)*s3;

		double part1=0.0;
		double part2_1=0;
		double part2_2=0.0;
		double part2_3=0.0;
		double morse_deriv_ij=0;
		Point <dim> morse_force_ij(0.0, 0.0, 0.0);

	    if(interatomic_bond < r_cut_mod) {
		     part1 = (E1 * derivM(alpha1_mod, r_01_mod, inv_material_bond_stretch) +
                      E2 * derivM(alpha2_mod, r_02_mod, inv_material_bond_stretch)) * Psi(r_cut_mod, h_mod, inv_material_bond_stretch) +
                     (E1 * calculateM(r_01_mod, inv_material_bond_stretch, alpha1_mod)
                    + E2 * calculateM(r_02_mod, inv_material_bond_stretch, alpha2_mod) + delta) * psiDerivative(r_cut_mod, h_mod, inv_material_bond_stretch);
	    }

	    if(r_s1_mod > interatomic_bond) {
	    	part2_1 = s1_mod*repulsionControlDerivative(r_s1_mod, inv_material_bond_stretch);
	    }

	    if(r_s2_mod > interatomic_bond) {
	    	part2_2 = s2_mod*repulsionControlDerivative(r_s2_mod, inv_material_bond_stretch);
	    }

	    if(r_s3_mod > interatomic_bond) {
	    	part2_3 = s3_mod*repulsionControlDerivative(r_s3_mod, inv_material_bond_stretch);
	    }

	    morse_deriv_ij = part1 - material_bond_distance * (part2_1 + part2_2 + part2_3);
	    morse_force_ij = bond_direction * morse_deriv_ij;
	    atoma_force += morse_force_ij;
	)

	return atoma_force;
}

template<int dim> Point<dim> Force<dim>::configForceEamEmbedding(int i, double a, double beta1,
                                                                 double beta2, double r_03, double r_04, double F0, double F2,
                                                                 double q1, double q2, double q3, double q4, double Q1, double Q2,
                                                                 double h, double r_cut) {
	Point<dim> atoma_force(0.0, 0.0, 0.0);
	double atom_electron_density = this->electronDensity(i, a, beta1, beta2, r_03, r_04, h, r_cut);
	double atom_deriv_embedding_eng = derivativeEmbeddingEnergy(atom_electron_density, F0, F2, q1, q2, q3, q4, Q1, Q2, h, r_cut);
	double inv_h = 1.0 / h;
	double psi = 0.0;

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		Point<dim> Embedding_Force(0.0, 0.0, 0.0);
		double neighbor_electron_density=0;
		double neighbor_deriv_embedding_eng=0;
		double deriv_electron_density=0;
        auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
        auto material_bond_distance = material_delta.norm();
		double beta1_mod=pow(material_bond_distance,2)*beta1;
		double beta2_mod=material_bond_distance*beta2;
		double r_03_mod=(1.0/material_bond_distance)*r_03;
		double r_04_mod=(1.0/material_bond_distance)*r_04;
		double h_mod=(1.0/material_bond_distance)*h;
		double r_cut_mod=(1.0/material_bond_distance)*r_cut;
        auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
        double interatomic_distance = spatial_delta.norm();
		double material_bond_stretch = material_bond_distance / interatomic_distance;
		double inv_material_bond_stretch = 1.0 / material_bond_stretch;

		neighbor_electron_density=this->electronDensity(j, a, beta1_mod, beta2_mod, r_03_mod, r_04_mod, h_mod, r_cut_mod);
		neighbor_deriv_embedding_eng=derivativeEmbeddingEnergy(neighbor_electron_density,
		                                                       F0, F2,
		                                                       q1, q2,
					                                           q3, q4,
					                                           Q1, Q2,
		                                                       h_mod, r_cut_mod);


		deriv_electron_density= derivativeElectronDensity(inv_material_bond_stretch, a, beta1_mod,
		                                                  beta2_mod, r_03_mod, r_04_mod,
		                                                  h_mod, r_cut_mod);


		Point<dim> bond_direction = material_delta / material_bond_distance;
		Embedding_Force=bond_direction*(neighbor_deriv_embedding_eng+atom_deriv_embedding_eng)*deriv_electron_density;
		atoma_force=atoma_force+Embedding_Force;
	)

	return atoma_force;
}

template<int dim> Point<dim> Force<dim>::StillingerWeberThreeBodyForceH2H3(
    int i, int j, int k,
    double sigma_AlphaBeta, double sigma_AlphaGamma, double epsilon_AlphaBetaGamma,
    double gamma, double lambda_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

    auto interatomic_distance_AlphaBeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto interatomic_distance_AlphaGamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(i);
    auto interatomic_distance_BetaGamma = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(k);
    double distance_AlphaBeta = interatomic_distance_AlphaBeta.norm();
    double distance_AlphaGamma = interatomic_distance_AlphaGamma.norm();
    double distance_BetaGamma = interatomic_distance_BetaGamma.norm();

	Point<dim> SWThreeBodyForce(0.0, 0.0, 0.0);
	cout << distance_AlphaBeta<< "\n";

	auto normalVectorAlphaBeta = interatomic_distance_AlphaBeta / distance_AlphaBeta;
	auto normalVectorAlphaGamma = interatomic_distance_AlphaGamma / distance_AlphaGamma;
	double cosine_teta_BetaAlphaGamma = (pow(distance_AlphaBeta, 2) + pow(distance_AlphaGamma, 2) - pow(distance_BetaGamma, 2)) /
                                        (2 * distance_AlphaBeta * distance_AlphaGamma);
	double term1 = (gamma * sigma_AlphaBeta) / (distance_AlphaBeta - (a_AlphaBeta * sigma_AlphaBeta));
	double term2 = (gamma * sigma_AlphaGamma) / (distance_AlphaGamma - (a_AlphaGamma * sigma_AlphaGamma));
	double terms_sum = term1 + term2;
	double exponent = exp(terms_sum);
	double difference_cosines = (cosine_teta_BetaAlphaGamma - cosine_teta0) / 100;
	double difference_cosines_squared = difference_cosines * difference_cosines;
	
	Point<dim> term3(0.0, 0.0, 0.0);
	term3.SetXCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetXCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta));
	term3.SetYCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta));
	term3.SetZCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta));

	Point<dim> term4(0.0, 0.0, 0.0);
	term4.SetXCoord(term1 * normalVectorAlphaBeta.GetXCoord());
	term4.SetYCoord(term1 * normalVectorAlphaBeta.GetYCoord());
	term4.SetZCoord(term1 * normalVectorAlphaBeta.GetZCoord());

	SWThreeBodyForce.SetXCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetXCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetXCoord());
	SWThreeBodyForce.SetYCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetYCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetYCoord());
	SWThreeBodyForce.SetZCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetZCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetZCoord());

	return SWThreeBodyForce;
}

template<int dim> Point<dim> Force <dim>::StillingerWeberThreeBodyForceH1(
    int j, int i, int k,
    double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
    double gamma, double lambda_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

    auto interatomic_distance_AlphaBeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto interatomic_distance_AlphaGamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(i);
    auto interatomic_distance_BetaGamma = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(k);
    double distance_AlphaBeta = interatomic_distance_AlphaBeta.norm();
    double distance_AlphaGamma = interatomic_distance_AlphaGamma.norm();
    double distance_BetaGamma = interatomic_distance_BetaGamma.norm();
	Point<dim> SWThreeBodyForce(0.0, 0.0, 0.0);

	if(distance_AlphaBeta < (sigma_AlphaBeta*a_AlphaBeta) && distance_AlphaGamma < (sigma_AlphaGamma*a_AlphaGamma)) {
        auto normalVectorAlphaBeta = interatomic_distance_AlphaBeta / distance_AlphaBeta;
        auto normalVectorAlphaGamma = interatomic_distance_AlphaGamma / distance_AlphaGamma;
        auto normalVectorBetaGamma = interatomic_distance_BetaGamma / distance_BetaGamma;

		double cosine_teta_BetaAlphaGamma=(pow(distance_AlphaBeta,2)+pow(distance_AlphaGamma,2)-pow(distance_BetaGamma,2))
						                   /(2*distance_AlphaBeta*distance_AlphaGamma);
		double term1=(gamma*sigma_AlphaBeta)/(distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta));
		double term2=(gamma*sigma_AlphaGamma)/(distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma));
		double exponent=exp(term1)*exp(term2);
		float difference_cosines=cosine_teta_BetaAlphaGamma-cosine_teta0;
		double difference_cosines_squared=difference_cosines*difference_cosines;
		Point <dim> term3(0.,0.,0.);
		
		term3.SetXCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetXCoord()/distance_AlphaBeta
				         -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetXCoord()/distance_AlphaGamma)
						+(normalVectorAlphaBeta.GetXCoord()/distance_AlphaGamma+
						  normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta));
	
		term3.SetYCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta
				        -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetYCoord()/distance_AlphaGamma)
						+(normalVectorAlphaBeta.GetYCoord()/distance_AlphaGamma+
						normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta));
	
	
		term3.SetZCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta
				        -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetZCoord()/distance_AlphaGamma)
						+(normalVectorAlphaBeta.GetZCoord()/distance_AlphaGamma+
						normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta));

		Point <dim> term4(0.,0.,0.);
		double prefactor1=0.;
		double prefactor2=0.;
	
		prefactor1=sigma_AlphaBeta/((distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta))*
				                      (distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta)));
	
		prefactor2=sigma_AlphaGamma/((distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma))*
				                       (distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma)));
		
		term4.SetXCoord(prefactor1*normalVectorAlphaBeta.GetXCoord()+prefactor2*normalVectorAlphaGamma.GetXCoord());
		term4.SetYCoord(prefactor1*normalVectorAlphaBeta.GetYCoord()+prefactor2*normalVectorAlphaGamma.GetYCoord());
		term4.SetZCoord(prefactor1*normalVectorAlphaBeta.GetZCoord()+prefactor2*normalVectorAlphaGamma.GetZCoord());

		Energy<dim> energy(this->atoms);
		double threebodyEng=energy.stillingerWeberThreeBody(j, i, k,
                                                            sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma, gamma,
                                                            lambda_AlphaBetaGamma, cosine_teta0, a_AlphaBeta, a_AlphaGamma);
		
		SWThreeBodyForce.SetXCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetXCoord()+
				                   gamma*threebodyEng*term4.GetXCoord());
		SWThreeBodyForce.SetYCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetYCoord()+
				                    gamma*threebodyEng*term4.GetYCoord());
		SWThreeBodyForce.SetZCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetZCoord()+
								   gamma*threebodyEng*term4.GetZCoord());

		Point <dim> force_atom_beta (0.,0.,0.);
		force_atom_beta.SetXCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
                                  *exponent*difference_cosines*
								  normalVectorAlphaBeta.GetXCoord()*(1/distance_AlphaGamma-(cosine_teta_BetaAlphaGamma/distance_AlphaBeta))
	                              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetXCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								  normalVectorBetaGamma.GetXCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		force_atom_beta.SetYCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta
					              +normalVectorAlphaBeta.GetYCoord()/distance_AlphaGamma)
					              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetYCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetYCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		force_atom_beta.SetZCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta
					              +normalVectorAlphaBeta.GetZCoord()/distance_AlphaGamma)
					              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetZCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetZCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));

        this->atoms->setForce(j, this->atoms->getForce(j) + force_atom_beta);

		Point <dim> force_atom_gamma (0.,0.,0.);
		force_atom_gamma.SetXCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetXCoord()/distance_AlphaGamma
					              +normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta)
					              -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetXCoord())
								  +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetXCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		force_atom_gamma.SetYCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                               *exponent*difference_cosines*
				                   (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetYCoord()/distance_AlphaGamma
					               +normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta)
					               -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetYCoord())
								   +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   	normalVectorBetaGamma.GetYCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		force_atom_gamma.SetZCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetZCoord()/distance_AlphaGamma
					              +normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta)
					              -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetZCoord())
								  +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetZCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));

        this->atoms->setForce(k, this->atoms->getForce(k) + force_atom_gamma);
	}

	return(SWThreeBodyForce);
}

template<int dim> Point<dim> Force<dim>::StillingerWeberTwoBodyForce(
    int i, int j,
    double sigma_AlphaBeta, double epsilon_AlphaBeta,
    double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
    double q_AlphaBeta, double a_AlphaBeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto distance = spatial_delta.norm();
	double term1 = sigma_AlphaBeta / pow((distance - (a_AlphaBeta * sigma_AlphaBeta)), 2);
	Energy<dim> energy(this->atoms);
	double energy_ij = energy.stillingerWeberTwoBody(i, j, sigma_AlphaBeta, epsilon_AlphaBeta,
                                                     A_AlphaBeta, B_AlphaBeta, p_AlphaBeta, 
												     q_AlphaBeta, a_AlphaBeta);

    Point<dim> normalVectorAlphaBeta = spatial_delta / distance;
    double term2_nominator = 0.0;
    double term2_denominator = 0.0;
    term2_nominator=(p_AlphaBeta*B_AlphaBeta*pow((sigma_AlphaBeta/distance),(p_AlphaBeta)))/distance-
    		        (q_AlphaBeta*pow((sigma_AlphaBeta/distance),(q_AlphaBeta)))/distance;
    term2_denominator=B_AlphaBeta*pow((sigma_AlphaBeta/distance),p_AlphaBeta)
                      -pow((sigma_AlphaBeta/distance),q_AlphaBeta);
    double term2 = term2_nominator/term2_denominator;

    Point <dim> SWTwoBodyForce(0.0, 0.0, 0.0);
	if (distance < (sigma_AlphaBeta*a_AlphaBeta)) {
		SWTwoBodyForce.SetXCoord(energy_ij*normalVectorAlphaBeta.GetXCoord()*(term1+term2));
		SWTwoBodyForce.SetYCoord(energy_ij*normalVectorAlphaBeta.GetYCoord()*(term1+term2));
		SWTwoBodyForce.SetZCoord(energy_ij*normalVectorAlphaBeta.GetZCoord()*(term1+term2));
	}

	return(SWTwoBodyForce);
}

template<int dim> Point<dim> Force<dim>::ResultantSWTwoBodyForce(int i, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                                 double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                                                 double q_AlphaBeta, double a_AlphaBeta) {
	Point<dim> atom_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		atom_force += StillingerWeberTwoBodyForce(i, j,
		                                          sigma_AlphaBeta, epsilon_AlphaBeta,
		                                          A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                          q_AlphaBeta, a_AlphaBeta);
    )

	return atom_force;
}

template<int dim> Point<dim> Force<dim>::ResultantSWThreeBodyForce(
    int i, double sigma_AlphaBeta, double sigma_AlphaGamma, double epsilon_AlphaBetaGamma, double mu,
    double lambda_AlphaBetaGamma, double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point <dim> atomaThreeBodyForce(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS_3BODY(this->atoms, i, j, k,
        atomaThreeBodyForce += StillingerWeberThreeBodyForceH1(k, i, j,
                                                               sigma_AlphaBeta,
                                                               sigma_AlphaGamma,epsilon_AlphaBetaGamma,
                                                               mu, lambda_AlphaBetaGamma,
                                                               cosine_teta0, a_AlphaBeta, a_AlphaGamma);
    )

    auto force = atoms->getForce(i) + atomaThreeBodyForce;
    atoms->setForce(i, force);
	return force;
}

template<int dim> Point<dim> Force<dim>::deformationalBondSWForce( int i, int j,
                                                                   double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                                   double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
                                                                   double J_AlphaBeta, double a_AlphaBeta) {

    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_bond_distance = material_delta.norm();
	double sigma0_AlphaBeta=sigma_AlphaBeta/material_bond_distance;
	double epsilon0_AlphaBeta=epsilon_AlphaBeta/material_bond_distance;
	Point<dim> spatial_bond_normal = spatial_delta / distance;
	double spatial_stretch = distance / material_bond_distance;
	Energy<dim> energy(this->atoms);
	double energy_AlphaBeta = energy.deformationalSWtwoBody(i, j,
                                                            sigma_AlphaBeta, epsilon_AlphaBeta,
                                                            A_AlphaBeta, B_AlphaBeta, I_AlphaBeta,
                                                            J_AlphaBeta, a_AlphaBeta);

	double term1=(I_AlphaBeta*B_AlphaBeta/spatial_stretch)*pow((sigma0_AlphaBeta/spatial_stretch),I_AlphaBeta)
			     -(J_AlphaBeta/spatial_stretch)*pow((sigma0_AlphaBeta/spatial_stretch),J_AlphaBeta);
	double term2=B_AlphaBeta*pow((sigma0_AlphaBeta/spatial_stretch),I_AlphaBeta)
	             -pow((sigma0_AlphaBeta/spatial_stretch),J_AlphaBeta);
	double term3=term1/term2;
	double term4=sigma0_AlphaBeta/pow((spatial_stretch-a_AlphaBeta*sigma0_AlphaBeta),2);

	Point<dim> BondSWForce(0.,0.,0.);
	BondSWForce.SetXCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetXCoord());
	BondSWForce.SetYCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetYCoord());
	BondSWForce.SetZCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetZCoord());
	return BondSWForce;
}

template<int dim> Point<dim> Force<dim>::deformationalprSW3BodyForce(
    int i, int j, int k,
    double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
	double mu, double eta_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point<dim> parallelForce(0.0, 0.0, 0.0);
	Energy<dim> energy(this->atoms);
	double deformational_energy = energy.deformationalSWthreeBody(
        i, j, k, sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
        mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
    auto mat_delta_alphabeta = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(i);
    auto mat_delta_alphagamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(i);
    auto mat_delta_betagamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(j);
	double mat_dist_alphabeta = mat_delta_alphabeta.norm();
	double mat_dist_alphagamma = mat_delta_alphagamma.norm();
	double mat_dist_betagamma = mat_delta_betagamma.norm();
	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;
	double area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;
    auto spa_alphabeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto spa_alphagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(i);
    auto spa_betagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(j);
	double spa_dist_alphabeta = spa_alphabeta.norm();
	double spa_dist_alphagamma = spa_alphagamma.norm();
	double spa_dist_betagamma = spa_betagamma.norm();
	double spatial_stretch_AlphaBeta = spa_dist_alphabeta / mat_dist_alphabeta;
	double spatial_stretch_AlphaGamma = spa_dist_alphagamma / mat_dist_alphagamma;
	
	if(spatial_stretch_AlphaBeta<(a_AlphaBeta*sigma_0_AlphaBeta) && spatial_stretch_AlphaGamma<(a_AlphaGamma*sigma_0_AlphaGamma)) {
        double term1=(mu*sigma_0_AlphaBeta/pow((spatial_stretch_AlphaBeta-(a_AlphaBeta*sigma_0_AlphaBeta)),2))*(1/mat_dist_alphabeta);
        double term2=(mu*sigma_0_AlphaGamma/pow((spatial_stretch_AlphaGamma-(a_AlphaGamma*sigma_0_AlphaGamma)),2))*(1/mat_dist_alphagamma);

        Point<dim> spatial_norm_AlphaBeta = mat_delta_alphabeta / mat_dist_alphabeta;
        Point<dim> spatial_norm_AlphaGamma = mat_delta_alphagamma / mat_dist_alphagamma;

        parallelForce.SetXCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetXCoord()
                                 +term2*spatial_norm_AlphaGamma.GetXCoord()));
        parallelForce.SetYCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetYCoord()
                                +term2*spatial_norm_AlphaGamma.GetYCoord()));
        parallelForce.SetZCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetZCoord()
                                +term2*spatial_norm_AlphaGamma.GetZCoord()));
	}
	
	return parallelForce;
	
}

template<int dim> Point<dim> Force<dim>::deformationalcpSW3BodyForce(
    int i, int j, int k,
    double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
	double mu, double eta_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point<dim> Result(0.0, 0.0, 0.0);
    auto spatial_bond_vec_AlphaBeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto spatial_bond_vec_AlphaGamma = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spatial_bond_vec_BetaGamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(j);
	double spatial_bond_distance_alphabeta = spatial_bond_vec_AlphaBeta.norm();
	double spatial_bond_distance_alphagamma = spatial_bond_vec_AlphaGamma.norm();
	double spatial_bond_distance_betagamma = spatial_bond_vec_BetaGamma.norm();
    auto material_bond_vec_AlphaBeta = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(i);
    auto material_bond_vec_AlphaGamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(i);
    auto material_bond_vec_BetaGamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(j);
	double mat_dist_alphabeta = material_bond_vec_AlphaBeta.norm();
	double mat_dist_alphagamma = material_bond_vec_AlphaGamma.norm();
	double mat_dist_betagamma = material_bond_vec_BetaGamma.norm();
	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;
	double area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));
	double epslion_0_AlphaBetaGamma=epsilon_AlphaBetaGamma/area;
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;
	double spatial_stretch_AlphaBeta = spatial_bond_distance_alphabeta / mat_dist_alphabeta; 
	double spatial_stretch_AlphaGamma = spatial_bond_distance_alphagamma / mat_dist_alphagamma;

	Point<dim> normalVectorAlphaBeta = spatial_bond_vec_AlphaBeta / spatial_bond_distance_alphabeta;
	double spa_AlphaBeta_x=normalVectorAlphaBeta.GetXCoord();
	double spa_AlphaBeta_y=normalVectorAlphaBeta.GetYCoord();
	double spa_AlphaBeta_z=normalVectorAlphaBeta.GetZCoord();

	Point<dim> normalVectorAlphaGamma = spatial_bond_vec_AlphaGamma / spatial_bond_distance_alphagamma;
	double spa_AlphaGamma_x=normalVectorAlphaGamma.GetXCoord();
	double spa_AlphaGamma_y=normalVectorAlphaGamma.GetYCoord();
	double spa_AlphaGamma_z=normalVectorAlphaGamma.GetZCoord();

	Point<dim> normalVectorBetaGamma = spatial_bond_vec_BetaGamma / spatial_bond_distance_betagamma;
	double spa_BetaGamma_x=normalVectorBetaGamma.GetXCoord();
	double spa_BetaGamma_y=normalVectorBetaGamma.GetYCoord();
	double spa_BetaGamma_z=normalVectorBetaGamma.GetZCoord();

	Point<dim> normalVectorAlphaBetaMat = material_bond_vec_AlphaBeta / mat_dist_alphabeta;
	double trace_matAlphaBeta=0;
	trace_matAlphaBeta=normalVectorAlphaBetaMat.GetXCoord()+normalVectorAlphaBetaMat.GetYCoord()+normalVectorAlphaBetaMat.GetZCoord();

	Point<dim> normalVectorAlphaGammaMat = material_bond_vec_AlphaGamma / mat_dist_alphagamma;
	double trace_matAlphaGamma=0;
	trace_matAlphaGamma=normalVectorAlphaGammaMat.GetXCoord()+normalVectorAlphaGammaMat.GetYCoord()+normalVectorAlphaGammaMat.GetZCoord();

	double cosine_teta_BetaAlphaGamma=normalVectorAlphaBeta.GetXCoord()*normalVectorAlphaGamma.GetXCoord()
			                          +normalVectorAlphaBeta.GetYCoord()*normalVectorAlphaGamma.GetYCoord()
									  +normalVectorAlphaBeta.GetZCoord()*normalVectorAlphaGamma.GetZCoord();
	
	double exponet_AlphaBeta=0.;
	exponet_AlphaBeta=exp(mu*sigma_0_AlphaBeta/(spatial_stretch_AlphaBeta-(a_AlphaBeta*sigma_0_AlphaBeta)));
	
	double exponet_AlphaGamma=0.;
	exponet_AlphaGamma=exp(mu*sigma_0_AlphaGamma/(spatial_stretch_AlphaGamma-(a_AlphaGamma*sigma_0_AlphaGamma)));

	double spatial_bond_vec_AlphaGamma_x = spatial_bond_vec_AlphaGamma.GetXCoord();
	double spatial_bond_vec_AlphaGamma_y = spatial_bond_vec_AlphaGamma.GetYCoord();
	double spatial_bond_vec_AlphaGamma_z = spatial_bond_vec_AlphaGamma.GetZCoord();

	double spatial_bond_vec_AlphaBeta_x = spatial_bond_vec_AlphaBeta.GetXCoord();
	double spatial_bond_vec_AlphaBeta_y = spatial_bond_vec_AlphaBeta.GetYCoord();
	double spatial_bond_vec_AlphaBeta_z = spatial_bond_vec_AlphaBeta.GetZCoord();

	double material_bond_vec_AlphaBeta_x=material_bond_vec_AlphaBeta.GetXCoord();
	double material_bond_vec_AlphaBeta_y=material_bond_vec_AlphaBeta.GetYCoord();
	double material_bond_vec_AlphaBeta_z=material_bond_vec_AlphaBeta.GetZCoord();
	
	double material_bond_vec_AlphaGamma_x=material_bond_vec_AlphaGamma.GetXCoord();
	double material_bond_vec_AlphaGamma_y=material_bond_vec_AlphaGamma.GetYCoord();
	double material_bond_vec_AlphaGamma_z=material_bond_vec_AlphaGamma.GetZCoord();
	
	double coefficient=2*eta_AlphaBetaGamma*epslion_0_AlphaBetaGamma*area*
			            exponet_AlphaBeta*exponet_AlphaGamma
						*(cosine_teta_BetaAlphaGamma-cosine_teta0);

	if(spatial_stretch_AlphaBeta < a_AlphaBeta*sigma_0_AlphaBeta && spatial_stretch_AlphaGamma < a_AlphaGamma*sigma_0_AlphaGamma) {
		Result.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_x/spatial_bond_distance_alphabeta+
				                      spa_AlphaGamma_x/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_x/spatial_bond_distance_alphabeta)
				                      +(spa_AlphaBeta_x/spatial_bond_distance_alphagamma)));
				
		
		Result.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_y/spatial_bond_distance_alphabeta+
                                      spa_AlphaGamma_y/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_y/spatial_bond_distance_alphabeta)
                                      +(spa_AlphaBeta_y/spatial_bond_distance_alphagamma)));
		
		Result.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_z/spatial_bond_distance_alphabeta+
                                      spa_AlphaGamma_z/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_z/spatial_bond_distance_alphabeta)
                                      +(spa_AlphaBeta_z/spatial_bond_distance_alphagamma)));

		Point <dim> deform_force_atom_beta(0.,0.,0.);

		deform_force_atom_beta.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_x/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_x/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_x);

		deform_force_atom_beta.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_y/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_y/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_y);

		deform_force_atom_beta.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_z/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_z/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_z);

        this->atoms->setDeformForce(j , this->atoms->getDeformForce(j) + deform_force_atom_beta);

		Point <dim> deform_force_atom_gamma(0.,0.,0.);
		deform_force_atom_gamma.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_x/spatial_bond_distance_alphagamma)
				                          +(spa_AlphaGamma_x/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);
		
		deform_force_atom_gamma.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_y/spatial_bond_distance_alphagamma)
                                          +(spa_AlphaGamma_y/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);
		
		deform_force_atom_gamma.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_z/spatial_bond_distance_alphagamma)
                                          +(spa_AlphaGamma_z/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);

        this->atoms->setDeformForce(k, this->atoms->getDeformForce(k) + deform_force_atom_gamma);
	}
	
	return Result;
}


template<int dim> Point<dim> Force<dim>::configurationalBondSWForce(
    int i, int j,
    double sigma_AlphaBeta, double epsilon_AlphaBeta,
    double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
    double J_AlphaBeta, double a_AlphaBeta) {

    auto spa_delta_alphabeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
	double spa_dist_alphabeta = spa_delta_alphabeta.norm();
    auto mat_delta_alphabeta = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(i);
	double mat_dist_alphabeta = mat_delta_alphabeta.norm();
	double sigma_0_AlphaBeta = sigma_AlphaBeta / mat_dist_alphabeta;
	double epsilon_0_AlphaBeta = epsilon_AlphaBeta / mat_dist_alphabeta;
	double material_AlphaBeta = mat_dist_alphabeta / spa_dist_alphabeta;

	Point<dim> material_bond_normal = mat_delta_alphabeta / spa_dist_alphabeta;
	Energy<dim> energy(this->atoms);
	double Config2BodyEng = energy.configurationalSWtwoBody(i, j, sigma_AlphaBeta, epsilon_AlphaBeta,
                                                            A_AlphaBeta,  B_AlphaBeta, I_AlphaBeta, J_AlphaBeta, a_AlphaBeta);

	double term0=1/material_AlphaBeta;
	double term1=(B_AlphaBeta*I_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),I_AlphaBeta)
	             -J_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),J_AlphaBeta))
	              /(B_AlphaBeta*material_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),I_AlphaBeta)
                    -material_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),J_AlphaBeta));
	double term2=sigma_0_AlphaBeta/pow((1-(material_AlphaBeta*a_AlphaBeta*sigma_0_AlphaBeta)),2);
	double term3=-Config2BodyEng*(term0+term1+term2);

	Point<dim> config_BondSWForce(0.,0.,0.);
	config_BondSWForce.SetXCoord(term3*material_bond_normal.GetXCoord());
	config_BondSWForce.SetYCoord(term3*material_bond_normal.GetYCoord());
	config_BondSWForce.SetZCoord(term3*material_bond_normal.GetZCoord());
	return config_BondSWForce;
	
}


template<int dim> Point<dim> Force<dim>::configurationalcpSW3BodyForce(
    int i, int j, int k,
	double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
    double mu, double eta_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point<dim> config_cpSW3BodyForce(0.0, 0.0, 0.0);
    
    auto spa_alphabeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto spa_alphagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(i);
    auto spa_betagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(j);
	double spa_dist_alphabeta = spa_alphabeta.norm();
	double spa_dist_alphagamma = spa_alphagamma.norm();
	double spa_dist_betagamma = spa_betagamma.norm();
	Point<dim> spa_bond_normal_alphabeta = spa_alphabeta / spa_dist_alphabeta;
	Point<dim> spa_bond_normal_alphagamma = spa_alphagamma / spa_dist_alphagamma;
    auto mat_alphabeta = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(i);
    auto mat_alphagamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(i);
    auto mat_betagamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(j);
	double mat_dist_alphabeta = mat_alphabeta.norm();
	double mat_dist_alphagamma = mat_alphagamma.norm();
	double mat_dist_betagamma = mat_betagamma.norm();
	Point<dim> mat_bond_normal_alphabeta = mat_alphabeta / mat_dist_alphabeta;
	Point<dim> mat_bond_normal_alphagamma = mat_alphagamma / mat_dist_alphagamma;
	Point<dim> mat_bond_normal_betagamma = mat_betagamma / mat_dist_betagamma;
	double m=(spa_dist_alphabeta+spa_dist_alphagamma+spa_dist_betagamma)/2;
	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;
	double area=sqrt(m*abs(m-spa_dist_alphabeta)*abs(m-spa_dist_alphagamma)*abs(m-spa_dist_betagamma));
	double Area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));
	double ratio_areas=Area/area;
	Energy<dim> energy(this->atoms);
	double config_energy=energy.configurationalSWthreeBody(i, j, k, sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                           mu, eta_AlphaBetaGamma,
                                                           cosine_teta0, a_AlphaBeta, a_AlphaGamma);

	double epsilon_AlphaBetaGamma_t=epsilon_AlphaBetaGamma/area;
	double epsilon_AlphaBetaGamma_0=epsilon_AlphaBetaGamma/Area;
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;
	double material_stretch_AlphaBeta = mat_dist_alphabeta / spa_dist_alphabeta;
	double material_stretch_AlphaGamma = mat_dist_alphagamma / spa_dist_alphagamma;
	
	SecondTensor<double> rotation_mat_alphabeta(3,3,0.);
	rotation_mat_alphabeta = this->bond->RotationMatrix(spa_bond_normal_alphabeta, mat_bond_normal_alphabeta);
	SecondTensor<double> rotation_mat_alphabeta_transpose(3,3,0.);
	rotation_mat_alphabeta_transpose=rotation_mat_alphabeta.transpose();
	
	SecondTensor<double> rotation_mat_alphagamma(3,3,0.);
	rotation_mat_alphabeta = this->bond->RotationMatrix(spa_bond_normal_alphagamma, mat_bond_normal_alphagamma);
	SecondTensor<double> rotation_mat_alphagamma_transpose(3,3,0.);
	rotation_mat_alphagamma_transpose=rotation_mat_alphagamma.transpose();
	
	std::vector <double> spa_alphabeta_vec;
	spa_alphabeta_vec.push_back(spa_alphabeta.GetXCoord());
	spa_alphabeta_vec.push_back(spa_alphabeta.GetYCoord());
	spa_alphabeta_vec.push_back(spa_alphabeta.GetZCoord());
	
	std::vector <double> spa_alphagamma_vec;
	spa_alphagamma_vec.push_back(spa_alphagamma.GetXCoord());
	spa_alphagamma_vec.push_back(spa_alphagamma.GetYCoord());
	spa_alphagamma_vec.push_back(spa_alphagamma.GetZCoord());
	
	std::vector <double> mat_alphabeta_vec=rotation_mat_alphabeta_transpose*spa_alphabeta_vec;
	std::vector <double> mat_alphagamma_vec=rotation_mat_alphagamma_transpose*spa_alphagamma_vec;
	
	double term1 =  mat_alphabeta_vec[0] * mat_alphagamma_vec[0] +
                    mat_alphabeta_vec[1] * mat_alphagamma_vec[1] +
                    mat_alphabeta_vec[2] * mat_alphagamma_vec[2];
	
	double term2 = sqrt(pow(mat_alphabeta_vec[0], 2) + pow(mat_alphabeta_vec[1], 2) + pow(mat_alphabeta_vec[2], 2)) *
				   sqrt(pow(mat_alphagamma_vec[0], 2) + pow(mat_alphagamma_vec[1], 2) + pow(mat_alphagamma_vec[2], 2));

	double cosine_theta=term1 / term2;
	double diff_cosines = cosine_theta - cosine_teta0;
			
	double term3=material_stretch_AlphaBeta*mu*sigma_0_AlphaBeta
				 /(1-(material_stretch_AlphaBeta*a_AlphaBeta*sigma_0_AlphaBeta));

	double term4=material_stretch_AlphaGamma*mu*sigma_0_AlphaGamma
				 /(1-(material_stretch_AlphaGamma*a_AlphaGamma*sigma_0_AlphaGamma));
	
	if(material_stretch_AlphaBeta < a_AlphaBeta*sigma_0_AlphaBeta && material_stretch_AlphaGamma < a_AlphaGamma*sigma_0_AlphaGamma) {
		config_cpSW3BodyForce.SetXCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphagamma)))
										-(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetXCoord()
										+((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetXCoord()));

		config_cpSW3BodyForce.SetYCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphagamma)))
										 -(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetYCoord()
										 +((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetYCoord()));

		config_cpSW3BodyForce.SetZCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphagamma)))
										-(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetZCoord()
										+((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetZCoord()));
		
		Point<dim> config_force_atom_beta(0.0, 0.0, 0.0);
		config_force_atom_beta.SetXCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetXCoord())
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetXCoord()));	
		
		config_force_atom_beta.SetYCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetYCoord())
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetYCoord()));
		
		config_force_atom_beta.SetZCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetXCoord())
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetZCoord()));

        this->atoms->setForce(j, this->atoms->getForce(j) + config_force_atom_beta);

		Point<dim> config_force_atom_gamma(0.0, 0.0, 0.0);
		config_force_atom_gamma.SetXCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
				                          ((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
				                          -cosine_theta*(mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphagamma))
				                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
				                          /(2*Area))*mat_bond_normal_alphagamma.GetXCoord()));	

		config_force_atom_gamma.SetYCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
                                          ((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
                                          -cosine_theta*(mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphagamma)
										  -(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										  *mat_bond_normal_betagamma.GetYCoord())
                                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
                                          /(2*Area))*mat_bond_normal_alphagamma.GetYCoord()));

		config_force_atom_gamma.SetZCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
                                          ((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
                                          -cosine_theta*(mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphagamma))
										  -(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										  *mat_bond_normal_betagamma.GetZCoord()
                                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
                                          /(2*Area))*mat_bond_normal_alphagamma.GetZCoord()));

        this->atoms->setForce(k, this->atoms->getForce(k) + config_force_atom_gamma);
	}

	return config_cpSW3BodyForce;
}
																   
template<int dim> Point<dim> Force<dim>::configurationalprSW3BodyForce(
    int i, int j, int k,
    double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
    double mu, double eta_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point<dim> config_prSW3BodyForce(0.0, 0.0, 0.0);
	Energy<dim> energy(this->atoms);

	double configurational_energy = energy.configurationalSWthreeBody(
        i, j, k, sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
        mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
    auto spa_alphabeta = this->atoms->getSpatialPosition(j) - this->atoms->getSpatialPosition(i);
    auto spa_alphagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(i);
    auto spa_betagamma = this->atoms->getSpatialPosition(k) - this->atoms->getSpatialPosition(j);
	double spa_dist_alphabeta = spa_alphabeta.norm();
	double spa_dist_alphagamma = spa_alphagamma.norm();
	double spa_dist_betagamma = spa_betagamma.norm();
    auto mat_alphabeta = this->atoms->getMaterialPosition(j) - this->atoms->getMaterialPosition(i);
    auto mat_alphagamma = this->atoms->getMaterialPosition(k) - this->atoms->getMaterialPosition(i);
	double mat_dist_alphabeta = mat_alphabeta.norm();
	double mat_dist_alphagamma = mat_alphagamma.norm();
	Point <dim> material_bond_normal_AlphaBeta = mat_alphabeta / mat_dist_alphabeta;
	Point <dim> material_bond_normal_AlphaGamma = mat_alphagamma / mat_dist_alphagamma;
	double m=(spa_dist_alphabeta+spa_dist_alphagamma+spa_dist_betagamma)/2;
	double area=sqrt(m*abs(m-spa_dist_alphabeta)*abs(m-spa_dist_alphagamma)*abs(m-spa_dist_betagamma));
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;
	double material_stretch_AlphaBeta = mat_dist_alphabeta / spa_dist_alphabeta;
	double material_stretch_AlphaGamma = mat_dist_alphagamma / spa_dist_alphagamma;
	double term1=(mu*sigma_0_AlphaBeta)/((1-(material_stretch_AlphaBeta)*a_AlphaBeta*sigma_0_AlphaBeta));
	double term2=(mu*sigma_0_AlphaGamma)/((1-(material_stretch_AlphaGamma)*a_AlphaGamma*sigma_0_AlphaGamma));

	double config_energy=energy.configurationalSWthreeBody(i, j, k,
                                                           sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                           mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);

	config_prSW3BodyForce.SetXCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetXCoord()
			                                            +term1*material_bond_normal_AlphaGamma.GetXCoord()));
	config_prSW3BodyForce.SetYCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetYCoord()
                                                         +term2*material_bond_normal_AlphaGamma.GetYCoord()));
	config_prSW3BodyForce.SetZCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetZCoord()
                                                        +term2*material_bond_normal_AlphaGamma.GetZCoord()));

	Point <dim> config_forcepr_atom_beta(0.0, 0.0, 0.0);
	config_forcepr_atom_beta.SetXCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetXCoord()));
	config_forcepr_atom_beta.SetYCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetYCoord()));
	config_forcepr_atom_beta.SetZCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetZCoord()));
    this->atoms->setForce(j, this->atoms->getForce(j) + config_forcepr_atom_beta);

	Point <dim> config_forcepr_atom_gamma(0.0, 0.0, 0.0);
	config_forcepr_atom_gamma.SetXCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetXCoord()));
	config_forcepr_atom_gamma.SetYCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetYCoord()));
	config_forcepr_atom_gamma.SetZCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetZCoord()));
    this->atoms->setForce(k, this->atoms->getForce(k) + config_forcepr_atom_gamma);

	return config_prSW3BodyForce;
}

template<int dim> Point<dim> Force<dim>::ResultantConfigSWTwoBodyForce(int i, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                                       double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                                                       double q_AlphaBeta, double a_AlphaBeta) {
	Point<dim> atoma_force(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
		atoma_force += configurationalBondSWForce(i, j, sigma_AlphaBeta, epsilon_AlphaBeta,
		                                          A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                          q_AlphaBeta, a_AlphaBeta);
	)

	return atoma_force;
}

template<int dim> Point<dim> Force <dim>::ResultantConfigSWThreeBodyForce(
    int i, double sigma_AlphaBeta,
    double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
    double Gamma, double lambda_AlphaBetaGamma,
    double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma) {

	Point <dim> atomaThreeBodyForce(0.0, 0.0, 0.0);

    LOOP_OVER_ATOM_NEIGHBORS_3BODY(this->atoms, i, j, k,
        Point <dim> ThreeBodyForceH1(0.,0.,0.);
        ThreeBodyForceH1 = configurationalcpSW3BodyForce(i, j, k, sigma_AlphaBeta, sigma_AlphaGamma,epsilon_AlphaBetaGamma,
                                                         Gamma, lambda_AlphaBetaGamma, cosine_teta0, a_AlphaBeta, a_AlphaGamma) +
                           configurationalprSW3BodyForce(i, j, k, sigma_AlphaBeta, sigma_AlphaGamma,epsilon_AlphaBetaGamma,
                                                         Gamma, lambda_AlphaBetaGamma, cosine_teta0, a_AlphaBeta, a_AlphaGamma);

        atomaThreeBodyForce.SetXCoord(atomaThreeBodyForce.GetXCoord() + ThreeBodyForceH1.GetXCoord());
        atomaThreeBodyForce.SetYCoord(atomaThreeBodyForce.GetYCoord() + ThreeBodyForceH1.GetYCoord());
        atomaThreeBodyForce.SetZCoord(atomaThreeBodyForce.GetZCoord() + ThreeBodyForceH1.GetZCoord());
    )

    auto force = this->atoms->getConfigForce(i) + atomaThreeBodyForce;
    this->atoms->setConfigForce(i, force);
	return force;
}

template class Force<2>;
template class Force<3>;
