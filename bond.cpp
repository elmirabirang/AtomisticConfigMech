/*
 * bond.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: S.Elmira Birang.O
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */

#include "bond.h"
#include "atom.h"
#include "point.h"
#include "math.h"
#include "matrix.h"
#include "third_order_tensor.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

template <int dim>
Bond<dim>::Bond()
{}

template <int dim>
Bond<dim>::~Bond()
{}

template <int dim>
Point <dim> Bond <dim>::MaterialBondVec(Atom <dim> atoma, Atom <dim> atomb)
		{
	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> materialPj=atomj.GetMaterialPosition() ;
	Point <dim> materialPi= atomi.GetMaterialPosition();

	Point <dim> material_bondij=materialPi-materialPj;

	return (material_bondij);

		}

template <int dim>
Point <dim> Bond <dim> :: SpatialBondVec(Atom <dim> atoma, Atom <dim> atomb)
		{
	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> spatialPi=atomi.GetSpatialPosition();
	Point <dim> spatialPj=atomj.GetSpatialPosition();

	Point <dim> spatial_bondij=spatialPi-spatialPj;

	return(spatial_bondij);

		}

template <int dim>
Point <dim> Bond <dim> :: InitialBondVec(Atom <dim> atoma, Atom <dim> atomb)
		{
	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> initialPi=atomi.GetInitialPosition();
	Point <dim> initialPj=atomj.GetInitialPosition();

	Point <dim> initial_bondij=initialPi-initialPj;

	return(initial_bondij);

		}




template <int dim>
double Bond <dim> :: MaterialBondDistance (Atom <dim> atoma, Atom <dim> atomb)
		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> material_bond_ij=MaterialBondVec(atomi,atomj);
	double material_bond_distance=material_bond_ij.PointNorm();

	return(material_bond_distance);


		}



template <int dim>
double Bond <dim> :: SpatialBondDistance(Atom <dim> atoma, Atom <dim> atomb)
		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> spatial_bond_ij=SpatialBondVec(atomi,atomj);

	double spatial_bond_distance=spatial_bond_ij.PointNorm();

	return(spatial_bond_distance);

		}

template <int dim>
double Bond <dim> :: InitialBondDistance(Atom <dim> atoma, Atom <dim> atomb)
		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> initial_bond_ij=InitialBondVec(atomi,atomj);

	double initial_bond_distance=initial_bond_ij.PointNorm();

	return(initial_bond_distance);

		}



template <int dim>
Point <dim> Bond <dim>::MaterialBondNormal(Atom <dim> atoma,Atom <dim> atomb)
		{
	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> material_bond_ij=MaterialBondVec(atomi,atomj);
	double material_bond_distance=MaterialBondDistance(atomi,atomj);

	Point <dim> material_bond_normal=material_bond_ij*(1/material_bond_distance);

	return(material_bond_normal);
		}



template <int dim>
Point <dim> Bond <dim>::SpatialBondNormal(Atom <dim> atoma,Atom <dim> atomb)
		{
	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> spatial_bond_ij=SpatialBondVec(atomi,atomj);
	double spatial_bond_distance=SpatialBondDistance(atomi,atomj);

	Point <dim> spatial_bond_normal=spatial_bond_ij*(1/spatial_bond_distance);

	return(spatial_bond_normal);
		}



template <int dim>
double Bond <dim>::MaterialBondStretch(Atom <dim> atoma, Atom <dim> atomb)

{

	this -> atomi=atoma;
	this -> atomj=atomb;

	double material_bond=MaterialBondDistance(atomi,atomj);
	double spatial_bond=SpatialBondDistance(atomi,atomj);

	double material_bond_stretch=material_bond/spatial_bond;

	return(material_bond_stretch);

}

template <int dim>
double Bond <dim>::SpatialBondStretch(Atom <dim> atoma, Atom <dim> atomb)

{

	this -> atomi=atoma;
	this -> atomj=atomb;

	double material_bond=MaterialBondDistance(atomi,atomj);
	double spatial_bond=SpatialBondDistance(atomi,atomj);

	double spatial_bond_stretch=spatial_bond/material_bond;

	return(spatial_bond_stretch);

}

template <int dim>
double Bond <dim>::InitialBondStretch(Atom <dim> atoma, Atom <dim> atomb)

{

	this -> atomi=atoma;
	this -> atomj=atomb;

	double material_bond=MaterialBondDistance(atomi,atomj);
	double initial_bond=InitialBondDistance(atomi,atomj);

	double spatial_bond_stretch=material_bond/initial_bond;

	return(spatial_bond_stretch);

}




template <int dim>
Point <dim> Bond <dim>:: ResultantBondVec_Material(Atom <dim> atoma)
		{
	
	this -> atomi=atoma;

	Point <dim> bonds_i;
	Point <dim> material_position=atomi.GetMaterialPosition();


	int dimension=material_position.GetDim();

	if(dimension==2)
	{
		bonds_i.SetXCoord(0.0);
		bonds_i.SetYCoord(0.0);
	}
	
	else
	{
		bonds_i.SetXCoord(0.0);
		bonds_i.SetYCoord(0.0);
		bonds_i.SetZCoord(0.0);
		
	}


	vector < Atom <dim> > neighbors=atomi.BondNeighbor();
	typedef typename vector < Atom <dim> >::iterator Nghbr;

	for(Nghbr neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{
		Point <dim> bond_ij=MaterialBondVec(atomi,(*neighbor));
		bonds_i=bonds_i+bond_ij;

	}

	return(bonds_i);

		}



template <int dim>
Point <dim> Bond <dim>:: ResultantBondVec_Spatial(Atom <dim> atoma)
		{
	this -> atomi=atoma;

	Point <dim> bondspa_i;
	Point <dim> spatial_position=atomi.GetSpatialPosition();


	int dimension=spatial_position.GetDim();

	if(dimension==2)
	{
		bondspa_i.SetXCoord(0.0);
		bondspa_i.SetYCoord(0.0);
	}

	else
	{
		bondspa_i.SetXCoord(0.0);
		bondspa_i.SetYCoord(0.0);
		bondspa_i.SetZCoord(0.0);

	}

	vector < Atom <dim> > bneighbors=atomi.BondNeighbor();
	typedef typename vector < Atom <dim> >::iterator Neighbr;

	for(Neighbr neighbr=bneighbors.begin();neighbr!=bneighbors.end();++neighbr)
	{
		Point <dim> bondspa_ij=SpatialBondVec(atomi,(*neighbr));
		bondspa_i=bondspa_i+bondspa_ij;

	}

	return(bondspa_i);

		}

template <int dim>
Point <dim> Bond <dim>::RotationAxis(Atom <dim> atoma, Atom <dim> atomb)
{

	this->atomi=atoma;
	this->atomj=atomb;

	Point <dim> material_bond_vec=MaterialBondNormal(atomi, atomj);
	double material_bond_vec_x=material_bond_vec.GetXCoord();
	double material_bond_vec_y=material_bond_vec.GetYCoord();
	double material_bond_vec_z=material_bond_vec.GetZCoord();

	Point <dim> spatial_bond_vec=SpatialBondNormal(atomi, atomj);
	double spatial_bond_vec_x=spatial_bond_vec.GetXCoord();
	double spatial_bond_vec_y=spatial_bond_vec.GetYCoord();
	double spatial_bond_vec_z=spatial_bond_vec.GetZCoord();

	double rotation_axis_x=material_bond_vec_y*spatial_bond_vec_z-material_bond_vec_z*spatial_bond_vec_y;
	double rotation_axis_y=material_bond_vec_z*spatial_bond_vec_x-material_bond_vec_x*spatial_bond_vec_z;
	double rotation_axis_z=material_bond_vec_x*spatial_bond_vec_y-material_bond_vec_y*spatial_bond_vec_x;

	Point <dim> rotation_axis(0.,0.,0.);

	rotation_axis.SetXCoord(rotation_axis_x);
	rotation_axis.SetYCoord(rotation_axis_y);
	rotation_axis.SetZCoord(rotation_axis_z);

	return(rotation_axis);

}

template <int dim>
double Bond<dim>::RotationAngle(Atom <dim> atoma, Atom <dim> atomb)
{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> material_bond_vec(0.,0.,0.);
	material_bond_vec=MaterialBondNormal(atomi, atomj);
	double material_bond_vec_x=material_bond_vec.GetXCoord();
	double material_bond_vec_y=material_bond_vec.GetYCoord();
	double material_bond_vec_z=material_bond_vec.GetZCoord();

	Point <dim> spatial_bond_vec(0.,0.,0.);
	spatial_bond_vec=SpatialBondNormal(atomi, atomj);
	double spatial_bond_vec_x=spatial_bond_vec.GetXCoord();
	double spatial_bond_vec_y=spatial_bond_vec.GetYCoord();
	double spatial_bond_vec_z=spatial_bond_vec.GetZCoord();

	double rot_angle=0.;
	float project_normals=0.;

	project_normals=material_bond_vec_x*spatial_bond_vec_x
				    +material_bond_vec_y*spatial_bond_vec_y
					+material_bond_vec_z*spatial_bond_vec_z+0.0;
	
	rot_angle=acos(project_normals);

	return (rot_angle);

}

template <int dim>
SecondTensor <double> Bond<dim>::RotationMatrix(Atom <dim> atoma, Atom <dim> atomb)
{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Point <dim> rotation_axis=RotationAxis(atomi,atomj);
	double rotation_axis_x=rotation_axis.GetXCoord();
	double rotation_axis_y=rotation_axis.GetYCoord();
	double rotation_axis_z=rotation_axis.GetZCoord();

	SecondTensor <double> outer_product_rotation_axis(3,3,0.);

	outer_product_rotation_axis(0,0)=rotation_axis_x*rotation_axis_x;

	outer_product_rotation_axis(0,1)=rotation_axis_x*rotation_axis_y;
	outer_product_rotation_axis(0,2)=rotation_axis_x*rotation_axis_z;


	outer_product_rotation_axis(1,0)=rotation_axis_y*rotation_axis_x;
	outer_product_rotation_axis(1,1)=rotation_axis_y*rotation_axis_y;
	outer_product_rotation_axis(1,2)=rotation_axis_y*rotation_axis_z;


	outer_product_rotation_axis(2,0)=rotation_axis_z*rotation_axis_x;
	outer_product_rotation_axis(2,1)=rotation_axis_z*rotation_axis_y;
	outer_product_rotation_axis(2,2)=rotation_axis_z*rotation_axis_z;


	SecondTensor <double> cross_product_matrix(3,3,0.);

	cross_product_matrix(0,0)=0.;
	cross_product_matrix(0,1)=rotation_axis_z;
	cross_product_matrix(0,2)=-1*rotation_axis_y;

	cross_product_matrix(1,0)=-1*rotation_axis_z;
	cross_product_matrix(1,1)=0.;
	cross_product_matrix(1,2)=rotation_axis_x;

	cross_product_matrix(2,0)=rotation_axis_y;
	cross_product_matrix(2,1)=-1*rotation_axis_x;
	cross_product_matrix(2,2)=0.;

	double rotation_angle=RotationAngle(atomi, atomj);
	double cosine_rotation_angle=cos(rotation_angle);
	double sine_rotation_angle=sin(rotation_angle);


	SecondTensor <double> RotationMatrix(3,3,0.);

	RotationMatrix(0,0)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_x,2);
	RotationMatrix(0,1)=(1-cosine_rotation_angle)*rotation_axis_x*rotation_axis_y-sine_rotation_angle*rotation_axis_z;
	RotationMatrix(0,2)=(1-cosine_rotation_angle)*rotation_axis_x*rotation_axis_z+sine_rotation_angle*rotation_axis_y;

	RotationMatrix(1,0)=(1-cosine_rotation_angle)*rotation_axis_y*rotation_axis_x+sine_rotation_angle*rotation_axis_z;
	RotationMatrix(1,1)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_y,2);
	RotationMatrix(1,2)=(1-cosine_rotation_angle)*rotation_axis_y*rotation_axis_z-sine_rotation_angle*rotation_axis_x;

	RotationMatrix(2,0)=(1-cosine_rotation_angle)*rotation_axis_z*rotation_axis_x-sine_rotation_angle*rotation_axis_y;
	RotationMatrix(2,1)=(1-cosine_rotation_angle)*rotation_axis_z*rotation_axis_y+sine_rotation_angle*rotation_axis_x;
	RotationMatrix(2,2)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_z,2);


	return(RotationMatrix);

}

template <int dim>
ThirdTensor <double> Bond <dim>::DerivRotMat_RotAxis(Atom <dim> atoma, Atom <dim> atomb)
{

	this -> atomi=atoma;
	this -> atomj=atomb;

	double rotation_angle=RotationAngle(atomi, atomj);
	Point <dim> rotation_axis=RotationAxis(atomi,atomj);

	double rotation_axis_x=rotation_axis.GetXCoord();
	double rotation_axis_y=rotation_axis.GetYCoord();
	double rotation_axis_z=rotation_axis.GetZCoord();

	double sine_rotation_angle=sin(rotation_angle);
	double cosine_rotation_angle=cos(rotation_angle);

	ThirdTensor <double> result(3,3,3,0.);

	ThirdTensor <double> Id_skew(3,3,3,0.);

	Id_skew(1,2,0)=1.*(-sine_rotation_angle);
	Id_skew(2,1,0)=-1.*(-sine_rotation_angle);

	Id_skew(0,2,1)=-1.*(-sine_rotation_angle);
	Id_skew(2,0,1)=1.*(-sine_rotation_angle);

	Id_skew(0,1,2)=1.*(-sine_rotation_angle);
	Id_skew(1,0,2)=-1.*(-sine_rotation_angle);

	ThirdTensor <double> deriv_rot_axis(3,3,3,0.);

	deriv_rot_axis(0,0,0)=2*rotation_axis_x*(1-cosine_rotation_angle);
	deriv_rot_axis(1,0,0)=rotation_axis_y*(1-cosine_rotation_angle);
	deriv_rot_axis(2,0,0)=rotation_axis_z*(1-cosine_rotation_angle);
	deriv_rot_axis(1,1,0)=rotation_axis_x*(1-cosine_rotation_angle);
	deriv_rot_axis(2,2,0)=rotation_axis_x*(1-cosine_rotation_angle);

	deriv_rot_axis(0,0,1)=rotation_axis_y*(1-cosine_rotation_angle);
	deriv_rot_axis(0,1,1)=rotation_axis_x*(1-cosine_rotation_angle);
	deriv_rot_axis(1,1,1)=2*rotation_axis_y*(1-cosine_rotation_angle);
	deriv_rot_axis(2,1,1)=rotation_axis_z*(1-cosine_rotation_angle);
	deriv_rot_axis(2,2,1)=rotation_axis_y*(1-cosine_rotation_angle);

	deriv_rot_axis(0,0,2)=rotation_axis_z*(1-cosine_rotation_angle);
	deriv_rot_axis(1,1,2)=rotation_axis_z*(1-cosine_rotation_angle);
	deriv_rot_axis(0,2,2)=rotation_axis_x*(1-cosine_rotation_angle);
	deriv_rot_axis(1,2,2)=rotation_axis_y*(1-cosine_rotation_angle);
	deriv_rot_axis(2,2,2)=2*rotation_axis_z*(1-cosine_rotation_angle);

	result=Id_skew+deriv_rot_axis;

	SecondTensor <double> pn(3,3,0.);

	Point <dim> mat_norm=MaterialBondNormal(atomi,atomj);
	double mat_norm_x=mat_norm.GetXCoord();
	double mat_norm_y=mat_norm.GetYCoord();
	double mat_norm_z=mat_norm.GetZCoord();

	pn(0,1)=-mat_norm_z;
	pn(0,2)=mat_norm_y;
	pn(1,0)=mat_norm_z;
	pn(1,2)=-mat_norm_x;
	pn(2,0)=-mat_norm_y;
	pn(2,1)=mat_norm_x;

	ThirdTensor <double> R(3,3,3,0.);

	for (int k=0; k<3; k++)
	{
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{

				R(i,j,k)+=result(i,j,k)*pn(j,i);

			}


		}

	}


	return R;

}

template <int dim>
ThirdTensor <double> Bond<dim>::DerivRotMat_RotAngle(Atom <dim> atoma, Atom <dim> atomb)
{

	this -> atomi=atoma;
	this -> atomj=atomb;

	SecondTensor <double> identity_mat(3,3,0.);

	identity_mat(0,0)=1.;
	identity_mat(1,1)=1.;
	identity_mat(2,2)=1.;

	Point <dim> rotation_axis=RotationAxis(atomi,atomj);

	double rotation_axis_x=rotation_axis.GetXCoord();
	double rotation_axis_y=rotation_axis.GetYCoord();
	double rotation_axis_z=rotation_axis.GetZCoord();

	SecondTensor <double> outer_product_rotation_axis(3,3,0.);

	outer_product_rotation_axis(0,0)=rotation_axis_x*rotation_axis_x;

	outer_product_rotation_axis(0,1)=rotation_axis_x*rotation_axis_y;
	outer_product_rotation_axis(0,2)=rotation_axis_x*rotation_axis_z;


	outer_product_rotation_axis(1,0)=rotation_axis_y*rotation_axis_x;
	outer_product_rotation_axis(1,1)=rotation_axis_y*rotation_axis_y;
	outer_product_rotation_axis(1,2)=rotation_axis_y*rotation_axis_z;


	outer_product_rotation_axis(2,0)=rotation_axis_z*rotation_axis_x;
	outer_product_rotation_axis(2,1)=rotation_axis_z*rotation_axis_y;
	outer_product_rotation_axis(2,2)=rotation_axis_z*rotation_axis_z;

	SecondTensor <double> cross_product_matrix(3,3,0.);

	cross_product_matrix(0,0)=0.;
	cross_product_matrix(0,1)=rotation_axis_z;
	cross_product_matrix(0,2)=-1*rotation_axis_y;

	cross_product_matrix(1,0)=-1*rotation_axis_z;
	cross_product_matrix(1,1)=0.;
	cross_product_matrix(1,2)=rotation_axis_x;

	cross_product_matrix(2,0)=rotation_axis_y;
	cross_product_matrix(2,1)=-1*rotation_axis_x;
	cross_product_matrix(2,2)=0.;

	SecondTensor <double> term1(3,3,0);
	SecondTensor <double> term2(3,3,0);

	double rotation_angle=RotationAngle(atomi, atomj);
	double sine_rotation_angle=sin(rotation_angle);
	double cosine_rotation_angle=cos(rotation_angle);


	term1=(identity_mat-outer_product_rotation_axis)*sine_rotation_angle;
	term2=cross_product_matrix*cosine_rotation_angle;

	Point <dim> mat_norm=MaterialBondNormal(atomi,atomj);
	double mat_norm_x=mat_norm.GetXCoord();
	double mat_norm_y=mat_norm.GetYCoord();
	double mat_norm_z=mat_norm.GetZCoord();

	std::vector<double> mat_norm_vec;
	mat_norm_vec.push_back(mat_norm_x);
	mat_norm_vec.push_back(mat_norm_y);
	mat_norm_vec.push_back(mat_norm_z);


	Point <dim> spa_norm_vec=SpatialBondNormal(atomi,atomj);

	double dot_prod_norms_x=mat_norm_x*spa_norm_vec.GetXCoord();
	double dot_prod_norms_y=mat_norm_y*spa_norm_vec.GetYCoord();
	double dot_prod_norms_z=mat_norm_z*spa_norm_vec.GetZCoord();

	float sum_dot_norms=0.0;
	sum_dot_norms=dot_prod_norms_x+dot_prod_norms_y+dot_prod_norms_z;
	
	ThirdTensor <double> Result(3,3,3,0.);
	
	if (sum_dot_norms!=1.0)
	{
	double term3=sqrt(1-pow(sum_dot_norms,2));

	SecondTensor <double> term4=(term1+term2)/term3;


	for (int k=0; k<3; k++)
	{
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{

				Result(i,j,k)=term4(i,j)*mat_norm_vec[k];


			}


		}

	}
	}

	return Result;

}


































