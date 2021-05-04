/*
 * LennardJones.cpp
 *
 *  Created on: Oct 20, 2020
 *      Author: S. Elmira Birang O.
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nürnberg
 */

#include <string>
#include <vector>
//---
#include "atom.h"
#include "bond.h"
#include "point.h"
#include "Lennard_Jones.h"
#include "math.h"
#include "matrix.h"

template<int dim> LennardJones<dim>::LennardJones() {}
template<int dim> LennardJones<dim>::~LennardJones() {}

// TODO: Provide the distance/normals/strecth terms as parameters and just compute the forces/energy in these kernels
template<int dim> double LennardJones<dim>::LennardJones_Classic(int i, int j, double epsilon, double sigma) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto spatial_bond_vec = spatial_delta / spatial_distance;
    double term1 = sigma / spatial_distance;
    double term2 = pow(term1, 6);
    double term3 = term2 * term2;
    double energy = 4 * epsilon * (term3 - term2);

    Point<dim> interatomic_force = spatial_bond_vec * (4 * epsilon / spatial_distance) * (-12 * term3 + 6 * term2);
    this->atoms->setForce(i, this->atoms->getForce(i) + interatomic_force);
    this->atoms->setForce(j, this->atoms->getForce(j) - interatomic_force);
    return energy;
}

template<int dim> double LennardJones<dim>::LennardJones_Configurational(int i, int j, double epsilon, double sigma) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto material_bond_vec = material_delta / material_distance;
    auto material_stretch = material_distance / spatial_distance;
    double sigma0 = sigma / material_distance;
    double epsilon0 = epsilon / material_distance;
    double term1 = sigma0 * material_stretch;
    double term2 = pow(term1, 6);
    double term3 = term2 * term2;
    double energy = 4 * epsilon0 * material_stretch * (term3 - term2);

    Point<dim> configurational_interatomic_force = material_bond_vec * 4 * epsilon0 * (13 * term3 - 7 * term2);
    this->atoms->setForce(i, this->atoms->getForce(i) + configurational_interatomic_force);
    this->atoms->setForce(j, this->atoms->getForce(j) - configurational_interatomic_force);
    return energy;
}

template<int dim> double LennardJones<dim>::LennardJones_Deformational(int i, int j, double epsilon, double sigma) {
    auto spatial_delta = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_distance = spatial_delta.norm();
    auto spatial_bond_vec = spatial_delta / spatial_distance;
    auto material_delta = this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
    auto material_distance = material_delta.norm();
    auto spatial_stretch = spatial_distance / material_distance;
    double sigma0 = sigma / material_distance;
    double epsilon0 = epsilon / material_distance;
    double term1 = sigma0 / spatial_stretch;
    double term2 = pow(term1, 6);
    double term3 = term2 * term2;
    double energy = 4 * epsilon * (term3 - term2);

    Point<dim> interatomic_force = spatial_bond_vec * (4 * epsilon / spatial_stretch) * (-12 * term3 + 6 * term2);
    this->atoms->setForce(i, this->atoms->getForce(i) + interatomic_force);
    this->atoms->setForce(j, this->atoms->getForce(j) - interatomic_force);
    return energy;
}
