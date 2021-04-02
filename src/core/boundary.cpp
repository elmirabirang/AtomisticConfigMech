/*
 * boundary.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 *
 */


#include "boundary.h"
#include <vector>

template<int dim> Boundary<dim>::Boundary() {}
template<int dim> Boundary<dim>::~Boundary() {}

template <int dim> void Boundary <dim>::SetBoundaryRegion(double x0, double y0, double x1, double y1 , int id) {
    this -> x0_limit=x0;
    this -> y0_limit=y0;
    this -> x1_limit=x1;
    this -> y1_limit=y1;
    this -> Bid=id;
}


template <int dim> void Boundary <dim>::SetBoundaryRegion(double x0, double y0, double z0,double x1, double y1, double z1, int id) {
    this -> x0_limit=x0;
    this -> y0_limit=y0;
    this -> z0_limit=z0;
    this -> x1_limit=x1;
    this -> y1_limit=y1;
    this -> z1_limit=z1;
    this -> Bid=id;
}

template <int dim> void Boundary <dim>::SetTriangBoundaryRegion(double x0, double y0, double x1, double y1, double x2, double y2, int id) {
	this -> x0_limit=x0;
	this -> y0_limit=y0;
	this -> x1_limit=x1;
	this -> y1_limit=y1;
	this -> x2_limit=x2;
	this-> y2_limit=y2;
	this ->Bid=id;
}

template <int dim> void Boundary <dim>:: SetTriangBoundaryRegion(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, int id) {
	this -> x0_limit=x0;
	this -> y0_limit=y0;
	this -> z0_limit=z0;
	this -> x1_limit=x1;
	this -> y1_limit=y1;
	this -> z1_limit=z1;
	this -> x2_limit=x2;
	this-> y2_limit=y2;
	this -> z2_limit=z2;
	this -> Bid=id;
}

template <int dim> double Boundary <dim>::GetBoundaryX0() {
	return(x0_limit);
}

template <int dim> double Boundary <dim>::GetBoundaryY0() {
    return (y0_limit);
}

template <int dim> double Boundary <dim>::GetBoundaryZ0() {
    return (z0_limit);
}

template <int dim> double Boundary <dim>::GetBoundaryX1() {
   return(x1_limit);

}

template <int dim> double Boundary <dim>::GetBoundaryY1() {
   return(y1_limit);

}

template <int dim> double Boundary <dim>::GetBoundaryZ1() {
   return(z1_limit);
}


template <int dim> double Boundary <dim>::GetBoundaryX2() {
	return(x2_limit);
}

template <int dim> double Boundary <dim>::GetBoundaryY2() {
    return (y2_limit);

}

template <int dim> double Boundary <dim>::GetBoundaryZ2() {
    return (z2_limit);
}

template <int dim> int Boundary <dim> ::GetBid() {
    return (Bid);
}

template <int dim> void Boundary<dim>:: ApplyBC(int id, char C, double f) {
	// TODO: develop it
	this -> Bid=id;
	this -> coordinate=C;
	this -> BC=f;
}

template class Boundary<2>;
template class Boundary<3>;
