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

template <int dim>
class Point
{
public:

	Point();
	Point(double x_value, double y_value, double z_value);
	Point(double x_value, double y_value);
	~Point();


	double GetXCoord();
	double GetYCoord();
	double GetZCoord();

	void SetXCoord(double x_value);
	void SetYCoord(double y_value);
	void SetZCoord(double z_value);

	Point <dim> operator-(Point <dim> p);
	Point <dim> operator*(double factor);
	Point <dim> operator+(Point <dim> P);



	double PointNorm();

	int GetDim();


private:

	double x;
	double y;
	double z;
};
