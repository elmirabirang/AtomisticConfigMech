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

#define LOOP_OVER_ATOMS(atoms, it, CODE) \
    { \
        for(int it = 0; it < atoms->getNumberOfAtoms(); ++it) { \
            CODE \
        } \
    }

#define LOOP_OVER_ATOMS_REVERSE(atoms, it, CODE) \
    { \
        for(int it = atoms->getNumberOfAtoms() - 1; it >= 0; --it) { \
            CODE \
        } \
    }

#define LOOP_OVER_ATOMS_PARALLEL(omp, atoms, it, CODE) \
    { \
        __pragma(omp) \
        for(int it = 0; it < atoms->getNumberOfAtoms(); ++it) { \
            CODE \
        } \
    }

#define LOOP_OVER_ATOM_NEIGHBORS(atoms, i, it, CODE) \
    { \
        auto neighbors = atoms->neighbors(i); \
        for(auto &it: neighbors) { \
            CODE \
        } \
        auto bond_neighbors = atoms->bondNeighbors(i); \
        for(auto &it: bond_neighbors) { \
            CODE \
        } \
    }

#define LOOP_OVER_ATOM_NEIGHBORS_CUSTOM(atoms, array, it, i, v, CODE) \
    { \
        auto neighbors = atoms->neighbors(i); \
        auto array = neighbors; \
        for(auto it = neighbors.begin(); it != neighbors.end();) { \
            auto v = *it; \
            CODE \
        } \
        auto bond_neighbors = atoms->bondNeighbors(i); \
        array = bond_neighbors; \
        for(auto it = bond_neighbors.begin(); it != bond_neighbors.end();) { \
            auto v = *it; \
            CODE \
        } \
    }

#define LOOP_OVER_ATOM_NEIGHBORS_3BODY(atoms, i, it1, it2, CODE) \
    { \
        auto neighbors = atoms->neighbors(i); \
        auto bond_neighbors = atoms->bondNeighbors(i); \
        for(auto &it1: neighbors) { \
            for(auto &it2: neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
            for(auto &it2: bond_neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
        } \
        for(auto &it1: bond_neighbors) { \
            for(auto &it2: bond_neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
        } \
    }

#define LOOP_OVER_ATOM_HALF_NEIGHBORS(atoms, i, it, CODE) \
    { \
        auto neighbors = atoms->neighbors(i); \
        for(auto &it: neighbors) { \
            CODE \
        } \
    }

#define LOOP_OVER_ATOM_HALF_NEIGHBORS_3BODY(atoms, i, it1, it2, CODE) \
    { \
        auto neighbors = atoms->neighbors(i); \
        auto bond_neighbors = atoms->bondNeighbors(i); \
        for(auto &it1: neighbors) { \
            for(auto &it2: neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
            for(auto &it2: bond_neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
        } \
        for(auto &it1: bond_neighbors) { \
            for(auto &it2: bond_neighbors) { \
                if(it1 < it2) { \
                    CODE \
                } \
            } \
        } \
    }

template<int dim> class Atoms {
public:
    Atoms() { this->capacity = 0; }
    ~Atoms() {}
    Atoms(size_t _capacity) { this->resizeCapacity(_capacity); }

    void setMaterialPosition(int atom, Point<dim> position) { this->material_positions[atom] = position; }
    Point<dim> getMaterialPosition(int atom) { return this->material_positions[atom]; }
    void setSpatialPosition(int atom, Point<dim> position) { this->spatial_positions[atom] = position; }
    Point<dim> getSpatialPosition(int atom) { return this->spatial_positions[atom]; }
    void setInitialPosition(int atom, Point<dim> position) { this->initial_positions[atom] = position; }
    Point<dim> getInitialPosition(int atom) { return this->initial_positions[atom]; }
    void setTemperature(int atom, double temperature) { this->temperatures[atom] = temperature; }
    double getTemperature(int atom) { return this->temperatures[atom]; }
    void setFrequency(int atom, double frequency) { this->frequencies[atom] = frequency; }
    double getFrequency(int atom) { return this->frequencies[atom]; }
    void setEntropy(int atom, double entropy) { this->entropies[atom] = entropy; }
    double getEntropy(int atom) { return this->entropies[atom]; }

    //This list of neighbors to calculate total potential energy of the system.
    void setNeighbors(int atom, std::vector<int> &neighbors) { this->atom_neighbors[atom] = neighbors; }
    const std::vector<int> &neighbors(int atom) const { return this->atom_neighbors[atom]; }
    void appendNeighbor(int atom, int neighbor) { this->atom_neighbors[atom].push_back(neighbor); }
    //The list of neighbors to calculate ResultantBondVec and ResultantBondVec_Spatial
    void setBondNeighbors(int atom, std::vector<int> &bond_neighbors) { this->bond_neighbors[atom] = bond_neighbors; }
    const std::vector<int> &bondNeighbors(int atom) const { return this->bond_neighbors[atom]; }
    void appendBondNeighbor(int atom, int neighbor) { this->bond_neighbors[atom].push_back(neighbor); }

    void setRegion(int atom, int boundary_id) { this->boundaries[atom] = boundary_id; }
    int getRegion(int atom) { return this->boundaries[atom]; }
    void setForce(int atom, Point<dim> force) { this->forces[atom] = force; }
    Point<dim> getForce(int atom) { return this->forces[atom]; }
    void setCell(int atom, int cell) { this->cells[atom] = cell; }
    int getCell(int atom) { return this->cells[atom]; }
    void setConfigForce(int atom, Point<dim> force) { this->config_forces[atom] = force; }
    Point<dim> getConfigForce(int atom) { return this->config_forces[atom]; }
    void setDeformForce(int atom, Point<dim> force) { this->deform_forces[atom] = force; }
    Point<dim> getDeformForce(int atom) { return this->deform_forces[atom]; }
    void setZone(int atom, int zone) { this->zones[atom] = zone; }
    int getZone(int atom) { return this->zones[atom]; }
    void setAtomicStress(int atom, std::vector<double> stress) { this->atomic_stresses[atom] = stress; }
    std::vector<double> getAtomicStress(int atom) { return this->atomic_stresses[atom]; }

    // TODO: Properly define this!
    Point<dim> getSpatialMeanPosition(int atom) { return this->spatial_positions[atom]; }

    int addAtom() {
        int atom_id = this->natoms++;

        if(this->natoms >= this->capacity) {
            this->resizeCapacity(this->natoms * 2);
        }

        return atom_id;
    }
    size_t getNumberOfAtoms() { return natoms; }
    size_t getCapacity() { return capacity; }

    void resizeCapacity(size_t _capacity) {
        this->capacity = _capacity;
        material_positions.resize(_capacity);
        spatial_positions.resize(_capacity);
        initial_positions.resize(_capacity);
        forces.resize(_capacity);
        config_forces.resize(_capacity);
        deform_forces.resize(_capacity);
        cells.resize(_capacity);
        boundaries.resize(_capacity);
        temperatures.resize(_capacity);
        frequencies.resize(_capacity);
        entropies.resize(_capacity);
        atom_neighbors.resize(_capacity);
        bond_neighbors.resize(_capacity);
        zones.resize(_capacity);
        atomic_stresses.resize(_capacity);
    }

private:
    size_t natoms;
    size_t capacity;
    std::vector<Point<dim>> material_positions;
    std::vector<Point<dim>> spatial_positions;
    std::vector<Point<dim>> initial_positions;
    std::vector<Point<dim>> forces;
    std::vector<Point<dim>> config_forces;
    std::vector<Point<dim>> deform_forces;
    std::vector<int> cells;
    std::vector<int> boundaries;
    std::vector<double> temperatures;
    std::vector<double> frequencies;
    std::vector<double> entropies;
    std::vector<std::vector<int>> atom_neighbors;
    std::vector<std::vector<int>> bond_neighbors;
    std::vector<int> zones;
    std::vector<std::vector<double>> atomic_stresses;
};
