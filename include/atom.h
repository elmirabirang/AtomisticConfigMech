/*
 * atom.h
 *
 *  Created on: Nov 16, 2018
 *      Author: S.Elmira Birang.O
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#pragma once

#include <vector>
//---
#include "point.h"

using namespace std;

template <int dim>
class Atom {

public:

    Atom() {}
    ~Atom() {}
    Atom(Point<dim> material_pos, Point<dim> current_pos, int id) :
        material_position(material_pos),
        spatial_position(current_pos),
        ID(id) {}

    void SetMaterialPosition(Point<dim> position) { this->material_position = position; }
    Point <dim> GetMaterialPosition() { return this->material_position; }
    void SetSpatialPosition(Point<dim> position) { this->spatial_position = position; }
    Point <dim> GetSpatialPosition() { return this->spatial_position; }
    void SetInitialPosition(Point<dim> position) { this->initial_position = position; }
    Point <dim> GetInitialPosition() { return this->initial_position; }
    void setTemperature(double temperature) { this->temperature = temperature; }
    double getTemperature() { return this->temperature; }
    void setFrequency(double temperature) { this->frequency = frequency; }
    double getFrequency() { return this->frequency; }
    void setEntropy(double temperature) { this->entropy = entropy; }
    double getEntropy() { return this->entropy; }

    void SetID(int id) { this->ID = id; }
    int GetID() { return this->ID; }

    //This list of neighbors to calculate total potential energy of the system.
    vector<Atom<dim> *> Neighbor() { return this->atom_neighbor; }
    //The list of neighbors to calculate ResultantBondVec and ResultantBondVec_Spatial
    vector<Atom<dim> *> BondNeighbor() { return this->bond_neighbor; }
    void SetNeighbors(vector<Atom<dim> *> neighbors) { this->atom_neighbor = neighbors; }
    void SetBondNeighbors(vector<Atom<dim> *> &bond_neighbors) { this->bond_neighbor = bond_neighbors; }

    void SetAtomRegion(int bID) { this->boundary_id = bID; }
    int GetAtomRegion() { return this->boundary_id; }
    void SetForce(Point<dim> force) { this->force = force; }
    Point<dim> GetForce() { return this->force; }
    void  SetCellID(int cell_id) { this->Cell_ID = cell_id; }
    int GetCellID() { return this->Cell_ID; }
    void SetConfigForce(Point<dim> force) { this->config_force = force; }
    Point <dim> GetConfigForce() { return this->config_force; }
    void SetDeformForce(Point<dim> force) { this->deform_force = force; }
    Point <dim> GetDeformForce() { return this->deform_force; }
    void setAtomZone(int zone) { this->zone = zone; }
    int getAtomZone() { return this->zone; }
    void setAtomicStress(vector<double> stress) { this->atomic_stress = stress; }
    vector<double> getAtomicStress() { return this->atomic_stress; }

    // TODO: Properly define this!
    Point<dim> getSpatialMeanPosition() { return this->spatial_position; }

 private:

   Point <dim> material_position;
   Point <dim> spatial_position;
   Point <dim> initial_position;

   int ID;
   int Cell_ID;
   double cutRadius;
   int boundary_id;

   double temperature;
   double frequency;
   double entropy;

   vector < Atom <dim>* > atom_neighbor;
   vector < Atom <dim>* > bond_neighbor;

   Point <dim> force;
   Point <dim> config_force;
   Point <dim> deform_force;

   int zone;

   vector <double> atomic_stress;

};

