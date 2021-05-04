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

template<int dim> class Point {
public:
    Point() {}
	Point(double x_value, double y_value, double z_value) : x(x_value), y(y_value), z(z_value) {}
	Point(double x_value, double y_value) : x(x_value), y(y_value) {}
	~Point() {}

	double GetXCoord() const { return this->x; }
	double GetYCoord() const { return this->y; }
	double GetZCoord() const { return this->z; }
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

	Point<dim> operator/(double factor) {
        double inv_factor = 1.0 / factor;
        if(dim == 2) {
            return Point<dim>(this->GetXCoord() * inv_factor, this->GetYCoord() * inv_factor);
        }

        return Point<dim>(this->GetXCoord() * inv_factor, this->GetYCoord() * inv_factor, this->GetZCoord() * inv_factor);
    }

    Point<dim> operator+=(double value) {
        this->x += value;
        this->y += value;
        if(dim == 3) {
            this->z += value;
        }

        return *this;
    }

	Point<dim> operator+=(const Point<dim>& p) {
        this->x += p.GetXCoord();
        this->y += p.GetYCoord();
        if(dim == 3) {
            this->z += p.GetZCoord();
        }

        return *this;
    }

	double norm() {
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

template<int dim> Point<dim> operator+(const double lhs, const Point<dim>& rhs) {
    return Point<dim>(lhs + rhs.GetXCoord(), lhs + rhs.GetYCoord(), lhs + rhs.GetZCoord());
} 
