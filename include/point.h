/*
 * point.h
 *
 *
 *  Created on: Nov 29, 2018
 *      Author: S.Elmira Birang. O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#pragma once

#include <math.h>

template <int dim> class Point {
public:
    Point() {}
	Point(double x_value, double y_value, double z_value) : x(x_value), y(y_value), z(z_value) {}
	Point(double x_value, double y_value) : x(x_value), y(y_value) {}
	~Point() {}

	double GetXCoord() { return this->x; }
	double GetYCoord() { return this->y; }
	double GetZCoord() { return this->z; }
	void SetXCoord(double x_value) { this->x = x_value; }
	void SetYCoord(double y_value) { this->y = y_value; }
	void SetZCoord(double z_value) { this->z = z_value; }

	Point<dim> operator-(Point<dim> p) { 
        if(dim == 2) {
            return Point<dim>(p.GetXCoord() - this->GetXCoord(), p.GetYCoord() - this->GetYCoord());
        }

        return Point<dim>(p.GetXCoord() - this->GetXCoord(), p.GetYCoord() - this->GetYCoord(), p.GetZCoord() - this->GetZCoord());
    }

	Point<dim> operator+(Point<dim> p) { 
        if(dim == 2) {
            return Point<dim>(p.GetXCoord() + this->GetXCoord(), p.GetYCoord() + this->GetYCoord());
        }

        return Point<dim>(p.GetXCoord() + this->GetXCoord(), p.GetYCoord() + this->GetYCoord(), p.GetZCoord() + this->GetZCoord());
    }

	Point<dim> operator*(double factor) {
        if(dim == 2) {
            return Point<dim>(this->GetXCoord() * factor, this->GetYCoord() * factor);
        }

        return Point<dim>(this->GetXCoord() * factor, this->GetYCoord() * factor, this->GetZCoord() * factor);
    }

	double PointNorm() {
        double normsq = this->x * this->x + this->y * this->y;
        if(dim == 3) {
            normsq += this->z * this->z;
        }

        return sqrt(normsq);
    }

	int GetDim() { return dim; }

private:
	double x;
	double y;
	double z;
};
