/*
 * source.cpp
 * This is the driving file of Point class
 *
 *  Created on: Nov 30, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 *
 */


#include "point.h"
#include <cmath>

template <int dim>
Point <dim>::Point()
{}

template <int dim>
Point <dim>:: Point(double x_value, double y_value, double z_value)
{
	this -> x=x_value;
	this -> y=y_value;
	this ->z=z_value;

}

template <int dim>
Point <dim>::Point (double x_value,double y_value)

		{

	this -> x=x_value;
	this -> y=y_value;

		}

template <int dim>
Point <dim>::~Point ()
{}

template <int dim>
double Point <dim>::GetXCoord()

		{
	return x;
		}

template <int dim>
double Point <dim>::GetYCoord()

		{
	return y;
		}

template <int dim>
double Point <dim>::GetZCoord()

		{
	return z;
		}

template <int dim>
void Point <dim>::SetXCoord(double x_value)
		{
	this -> x = x_value;
		}

template <int dim>
void Point <dim>::SetYCoord(double y_value)
		{
	this -> y = y_value;
		}

template <int dim>
void Point <dim>::SetZCoord(double z_value)
		{
	this -> z = z_value;
		}



template <int dim>
int Point <dim> :: GetDim()
{
  return(dim);
}

template <int dim>
double Point <dim> :: PointNorm()
{
	int dimension = this->GetDim();

	if (dimension==2)
	{
		x=this -> GetXCoord();
		y=this -> GetYCoord();

		double norm= sqrt(pow(x,2)+pow(y,2));
		return(norm);

	}

	else
	{

		double x=this -> GetXCoord();
		double y=this -> GetYCoord();
		double z=this -> GetZCoord();

		double norm=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
		return(norm);


	}

}


template <int dim>
Point <dim> Point <dim>::operator -(Point <dim> p)
{
	int dimension=this->GetDim();

	if (dimension==3){

		double x1=this->GetXCoord();
		double y1=this->GetYCoord();
		double z1=this->GetZCoord();

		double x2=p.GetXCoord();
		double y2=p.GetYCoord();
		double z2=p.GetZCoord();

		double x=x2-x1;
		double y=y2-y1;
		double z=z2-z1;

		Point <dim> point;
		point.SetXCoord(x);
		point.SetYCoord(y);
		point.SetZCoord(z);

		return (point);

	}

	else {

		double x1=this->GetXCoord();
		double y1=this->GetYCoord();

		double x2=p.GetXCoord();
		double y2=p.GetYCoord();

		double x=x2-x1;
		double y=y2-y1;

		Point <dim> point;
		point.SetXCoord(x);
		point.SetYCoord(y);

		return (point);

	}


}

template <int dim>
Point <dim> Point <dim> :: operator *(double factor)
{
int dimension=this->GetDim();

if (dimension==2)
{
	double x=this->GetXCoord();
	double y=this->GetYCoord();

	double xtimsfac=x*factor;
	double ytimsfac=y*factor;

	Point <dim> point;

	point.SetXCoord(xtimsfac);
	point.SetYCoord(ytimsfac);

	return(point);
}

else
{
	double x=this -> GetXCoord();
	double y=this -> GetYCoord();
	double z=this -> GetZCoord();

	double xtimsfac=x*factor;
	double ytimsfac=y*factor;
	double ztimsfac=z*factor;

	Point <dim> point;
	point.SetXCoord(xtimsfac);
	point.SetYCoord(ytimsfac);
	point.SetZCoord(ztimsfac);

	return(point);


}

}

template <int dim>
Point <dim> Point <dim> :: operator +(Point <dim> P)
{
	int dimension=this -> GetDim();
	
if (dimension==2)
{
	
	double x1=this -> GetXCoord();
	double y1=this -> GetYCoord();
	
	double x2=P.GetXCoord();
	double y2=P.GetYCoord();
	
	double x=x1+x2;
	double y=y1+y2;
	
	Point <dim> point;

	point.SetXCoord(x);
	point.SetYCoord(y);
	
	return(point);
	

}


else 
{
	
	double x1=this -> GetXCoord();
	double y1= this -> GetYCoord();
	double z1=this -> GetZCoord();
	
	double x2=P.GetXCoord();
	double y2=P.GetYCoord();
	double z2=P.GetZCoord();
	
	
	double x=x1+x2;
	double y=y1+y2;
	double z=z1+z2;
	
	Point <dim> point;

	point.SetXCoord(x);
	point.SetYCoord(y);
	point.SetZCoord(z);
	
	return(point);
	
}

}

template class Point<2>;
template class Point<3>;
