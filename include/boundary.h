/*
 * boundary.h
 *
 *  Created on: Dec 5, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */

#pragma once

#include <vector>

template <int dim>
class Boundary
{

public:

	Boundary();
	~Boundary();

	void SetBoundaryRegion(double x0,double y0, double x1, double y1, int id);
	void SetBoundaryRegion(double x0, double y0, double z0, double x1, double y1, double z1, int id);

	void SetTriangBoundaryRegion(double x0, double y0, double x1, double y1, double x2, double y2, int id);
	void SetTriangBoundaryRegion(double x0, double y0,double z0, double x1, double y1,double z1 ,
			                     double x2, double y2, double z2 ,int id);

	double GetBoundaryX0();
	double GetBoundaryY0();
	double GetBoundaryZ0();

	double GetBoundaryX1();
	double GetBoundaryY1();
	double GetBoundaryZ1();

	double GetBoundaryX2();
	double GetBoundaryY2();
	double GetBoundaryZ2();

	void ApplyBC(int id, char C, double f);

	int GetBid();

private:

	double x0_limit;
	double y0_limit;
	double z0_limit;

	double x1_limit;
	double y1_limit;
	double z1_limit;

	double x2_limit;
	double y2_limit;
	double z2_limit;

	char coordinate;
	double BC;
	int Bid;


};
