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

template<int dim> Point<dim> Bond<dim>::MaterialBondVec(int i, int j) {
	return this->atoms->getMaterialPosition(i) - this->atoms->getMaterialPosition(j);
}

template<int dim> Point<dim> Bond<dim>::SpatialBondVec(int i, int j) {
	return this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
}

template<int dim> Point<dim> Bond<dim>::InitialBondVec(int i, int j) {
	return this->atoms->getInitialPosition(i) - this->atoms->getInitialPosition(j);
}

template<int dim> double Bond<dim>::MaterialBondDistance(int i, int j) {
    return this->MaterialBondVec(i, j).norm();
}

template<int dim> double Bond<dim>::SpatialBondDistance(int i, int j) {
	return this->SpatialBondVec(i, j).norm();
}

template<int dim> double Bond<dim>::InitialBondDistance(int i, int j) {
	return this->InitialBondVec(i, j).norm();
}

template<int dim> Point<dim> Bond<dim>::MaterialBondNormal(int i, int j) {
	Point<dim> material_bond_ij = this->MaterialBondVec(i, j);
	return material_bond_ij / material_bond_ij.norm();
}

template<int dim> Point<dim> Bond<dim>::SpatialBondNormal(int i, int j) {
	Point<dim> spatial_bond_ij = this->SpatialBondVec(i, j);
    return spatial_bond_ij / spatial_bond_ij.norm();
}

template<int dim> double Bond<dim>::MaterialBondStretch(int i, int j) {
	return this->MaterialBondDistance(i, j) / this->SpatialBondDistance(i, j);
}

template<int dim> double Bond<dim>::SpatialBondStretch(int i, int j) {
	return this->SpatialBondDistance(i, j) / this->MaterialBondDistance(i, j);
}

template<int dim> double Bond<dim>::InitialBondStretch(int i, int j) {
	return this->MaterialBondDistance(i, j) / this->InitialBondDistance(i, j);
}

template<int dim> Point<dim> Bond<dim>::ResultantBondVec_Material(int i) {
    Point<dim> bonds_i(0.0, 0.0, 0.0);
    Point<dim> material_position = this->atoms->getMaterialPosition(i);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        Point<dim> bond_ij = material_position - this->atoms->getMaterialPosition(j);
        bonds_i = bonds_i + bond_ij;
    )

    return bonds_i;
}

template<int dim> Point<dim> Bond<dim>::ResultantBondVec_Spatial(int i) {
	Point<dim> bondspa_i(0.0, 0.0, 0.0);
	Point<dim> spatial_position = this->atoms->getSpatialPosition(i);

    LOOP_OVER_ATOM_NEIGHBORS(this->atoms, i, j,
        Point<dim> bondspa_ij = spatial_position - this->atoms->getSpatialPosition(j);
        bondspa_i = bondspa_i + bondspa_ij;
    )

    return bondspa_i;
}

template<int dim> Point<dim> Bond<dim>::RotationAxis(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec) {
	double spatial_bond_vec_x = spatial_bond_vec.GetXCoord();
	double spatial_bond_vec_y = spatial_bond_vec.GetYCoord();
	double spatial_bond_vec_z = spatial_bond_vec.GetZCoord();
	double material_bond_vec_x = material_bond_vec.GetXCoord();
	double material_bond_vec_y = material_bond_vec.GetYCoord();
	double material_bond_vec_z = material_bond_vec.GetZCoord();
	double rotation_axis_x = material_bond_vec_y * spatial_bond_vec_z - material_bond_vec_z * spatial_bond_vec_y;
	double rotation_axis_y = material_bond_vec_z * spatial_bond_vec_x - material_bond_vec_x * spatial_bond_vec_z;
	double rotation_axis_z = material_bond_vec_x * spatial_bond_vec_y - material_bond_vec_y * spatial_bond_vec_x;
	return Point<dim>(rotation_axis_x, rotation_axis_y, rotation_axis_z);
}

template<int dim> double Bond<dim>::RotationAngle(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec) {
	double spatial_bond_vec_x = spatial_bond_vec.GetXCoord();
	double spatial_bond_vec_y = spatial_bond_vec.GetYCoord();
	double spatial_bond_vec_z = spatial_bond_vec.GetZCoord();
	double material_bond_vec_x = material_bond_vec.GetXCoord();
	double material_bond_vec_y = material_bond_vec.GetYCoord();
	double material_bond_vec_z = material_bond_vec.GetZCoord();
	float project_normals = material_bond_vec_x * spatial_bond_vec_x +
                            material_bond_vec_y * spatial_bond_vec_y +
                            material_bond_vec_z * spatial_bond_vec_z;
	return acos(project_normals);
}

template<int dim> SecondTensor<double> Bond<dim>::RotationMatrix(Point<dim> spatial_bond_vec, Point<dim> material_bond_vec) {
	Point<dim> rotation_axis = RotationAxis(spatial_bond_vec, material_bond_vec);
	double rotation_axis_x = rotation_axis.GetXCoord();
	double rotation_axis_y = rotation_axis.GetYCoord();
	double rotation_axis_z = rotation_axis.GetZCoord();
	double rotation_angle = RotationAngle(spatial_bond_vec, material_bond_vec);

	SecondTensor<double> outer_product_rotation_axis(3, 3, 0.0);
	outer_product_rotation_axis(0,0)=rotation_axis_x*rotation_axis_x;
	outer_product_rotation_axis(0,1)=rotation_axis_x*rotation_axis_y;
	outer_product_rotation_axis(0,2)=rotation_axis_x*rotation_axis_z;
	outer_product_rotation_axis(1,0)=rotation_axis_y*rotation_axis_x;
	outer_product_rotation_axis(1,1)=rotation_axis_y*rotation_axis_y;
	outer_product_rotation_axis(1,2)=rotation_axis_y*rotation_axis_z;
	outer_product_rotation_axis(2,0)=rotation_axis_z*rotation_axis_x;
	outer_product_rotation_axis(2,1)=rotation_axis_z*rotation_axis_y;
	outer_product_rotation_axis(2,2)=rotation_axis_z*rotation_axis_z;

	SecondTensor<double> cross_product_matrix(3,3,0.);
	cross_product_matrix(0,0)=0.;
	cross_product_matrix(0,1)=rotation_axis_z;
	cross_product_matrix(0,2)=-1*rotation_axis_y;
	cross_product_matrix(1,0)=-1*rotation_axis_z;
	cross_product_matrix(1,1)=0.;
	cross_product_matrix(1,2)=rotation_axis_x;
	cross_product_matrix(2,0)=rotation_axis_y;
	cross_product_matrix(2,1)=-1*rotation_axis_x;
	cross_product_matrix(2,2)=0.;

	double cosine_rotation_angle=cos(rotation_angle);
	double sine_rotation_angle=sin(rotation_angle);

	SecondTensor<double> RotationMatrix(3,3,0.);
	RotationMatrix(0,0)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_x,2);
	RotationMatrix(0,1)=(1-cosine_rotation_angle)*rotation_axis_x*rotation_axis_y-sine_rotation_angle*rotation_axis_z;
	RotationMatrix(0,2)=(1-cosine_rotation_angle)*rotation_axis_x*rotation_axis_z+sine_rotation_angle*rotation_axis_y;
	RotationMatrix(1,0)=(1-cosine_rotation_angle)*rotation_axis_y*rotation_axis_x+sine_rotation_angle*rotation_axis_z;
	RotationMatrix(1,1)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_y,2);
	RotationMatrix(1,2)=(1-cosine_rotation_angle)*rotation_axis_y*rotation_axis_z-sine_rotation_angle*rotation_axis_x;
	RotationMatrix(2,0)=(1-cosine_rotation_angle)*rotation_axis_z*rotation_axis_x-sine_rotation_angle*rotation_axis_y;
	RotationMatrix(2,1)=(1-cosine_rotation_angle)*rotation_axis_z*rotation_axis_y+sine_rotation_angle*rotation_axis_x;
	RotationMatrix(2,2)=cosine_rotation_angle+(1-cosine_rotation_angle)*pow(rotation_axis_z,2);
	return RotationMatrix;
}

template<int dim> ThirdTensor<double> Bond<dim>::DerivRotMat_RotAxis(int i, int j) {
	Point<dim> spa_norm = this->SpatialBondNormal(i, j);
	Point<dim> mat_norm = this->MaterialBondNormal(i, j);
	double mat_norm_x = mat_norm.GetXCoord();
	double mat_norm_y = mat_norm.GetYCoord();
	double mat_norm_z = mat_norm.GetZCoord();

	double rotation_angle = RotationAngle(spa_norm, mat_norm);
	Point<dim> rotation_axis = RotationAxis(spa_norm, mat_norm);
	double rotation_axis_x = rotation_axis.GetXCoord();
	double rotation_axis_y = rotation_axis.GetYCoord();
	double rotation_axis_z = rotation_axis.GetZCoord();
	double sine_rotation_angle = sin(rotation_angle);
	double cosine_rotation_angle = cos(rotation_angle);

	ThirdTensor<double> result(3,3,3,0.);
	ThirdTensor<double> Id_skew(3,3,3,0.);
	Id_skew(1,2,0)=1.*(-sine_rotation_angle);
	Id_skew(2,1,0)=-1.*(-sine_rotation_angle);
	Id_skew(0,2,1)=-1.*(-sine_rotation_angle);
	Id_skew(2,0,1)=1.*(-sine_rotation_angle);
	Id_skew(0,1,2)=1.*(-sine_rotation_angle);
	Id_skew(1,0,2)=-1.*(-sine_rotation_angle);

	ThirdTensor<double> deriv_rot_axis(3,3,3,0.);
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

	result = Id_skew + deriv_rot_axis;

	SecondTensor<double> pn(3, 3, 0.);
	pn(0,1) = -mat_norm_z;
	pn(0,2) = mat_norm_y;
	pn(1,0) = mat_norm_z;
	pn(1,2) = -mat_norm_x;
	pn(2,0) = -mat_norm_y;
	pn(2,1) = mat_norm_x;

	ThirdTensor <double> R(3,3,3,0.);
	for(int k = 0; k < 3; k++) {
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				R(i, j, k) += result(i, j, k) * pn(j, i);
			}
		}
	}

	return R;
}

template<int dim> ThirdTensor<double> Bond<dim>::DerivRotMat_RotAngle(int i, int j) {
	SecondTensor<double> identity_mat(3,3,0.);
	identity_mat(0,0)=1.;
	identity_mat(1,1)=1.;
	identity_mat(2,2)=1.;

	Point<dim> spa_norm_vec = this->SpatialBondNormal(i, j);
	Point<dim> mat_norm = this->MaterialBondNormal(i, j);
	double mat_norm_x = mat_norm.GetXCoord();
	double mat_norm_y = mat_norm.GetYCoord();
	double mat_norm_z = mat_norm.GetZCoord();

	std::vector<double> mat_norm_vec;
	mat_norm_vec.push_back(mat_norm_x);
	mat_norm_vec.push_back(mat_norm_y);
	mat_norm_vec.push_back(mat_norm_z);
	double dot_prod_norms_x = mat_norm_x*spa_norm_vec.GetXCoord();
	double dot_prod_norms_y = mat_norm_y*spa_norm_vec.GetYCoord();
	double dot_prod_norms_z = mat_norm_z*spa_norm_vec.GetZCoord();

	Point<dim> rotation_axis=RotationAxis(spa_norm_vec, mat_norm);
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

	double rotation_angle = RotationAngle(spa_norm_vec, mat_norm);
	double sine_rotation_angle = sin(rotation_angle);
	double cosine_rotation_angle = cos(rotation_angle);
	term1 = (identity_mat - outer_product_rotation_axis) * sine_rotation_angle;
	term2 = cross_product_matrix * cosine_rotation_angle;

	float sum_dot_norms = 0.0;
	sum_dot_norms = dot_prod_norms_x + dot_prod_norms_y + dot_prod_norms_z;

	ThirdTensor<double> result(3, 3, 3, 0.);
    if(sum_dot_norms != 1.0) {
        double term3 = sqrt(1 - pow(sum_dot_norms, 2));
        SecondTensor<double> term4 = (term1 + term2) / term3;

        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                for(int k = 0; k < 3; k++) {
                    result(i, j, k) = term4(i, j) * mat_norm_vec[k];
                }
            }
        }
    }

    return result;
}

template class Bond<2>;
template class Bond<3>;
