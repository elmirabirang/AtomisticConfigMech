#include <vector>
#include <string>
//---
#include "atom.h"
#include "bond.h"
#include "point.h"
#include "modified_tersoff.h"
#include "math.h"
//#include "matrixx.h"

#define beta 3.0
#define alpha 1.80536502
#define h -0.38136087
#define eta 2.16152496
#define lam1 3.18011795
#define lam2 1.39343356
#define B 117.78072440
#define D 0.33090566
#define R 2.87478837
#define A 3198.51383127
#define n 1.98633876
#define c1 0.20123243
#define c2 614230.04310619
#define c3 996439.09714140
#define c4 3.33560562
#define c5 25.20963770

template<int dim> double ModifiedTersoff<dim>::ModiTers(int i, int j, int k){
    auto spatial_vector_ij = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(j);
    auto spatial_vector_ik = this->atoms->getSpatialPosition(i) - this->atoms->getSpatialPosition(k);
    auto spatial_distance_ij = spatial_vector_ij.norm();
    auto spatial_distance_ik = spatial_vector_ik.norm();
    double spatial_vector_ij_x = spatial_vector_ij.GetXCoord();
    double spatial_vector_ij_y = spatial_vector_ij.GetYCoord();
    double spatial_vector_ij_z = spatial_vector_ij.GetZCoord();
    double spatial_vector_ik_x = spatial_vector_ik.GetXCoord();
    double spatial_vector_ik_y = spatial_vector_ik.GetYCoord();
    double spatial_vector_ik_z = spatial_vector_ik.GetZCoord();
    double cos_theta = (spatial_vector_ij_x * spatial_vector_ik_x +
                        spatial_vector_ij_y * spatial_vector_ik_y +
                        spatial_vector_ij_z * spatial_vector_ik_z) /
                        spatial_distance_ij * spatial_distance_ik;

    double R_cons_ij = (spatial_distance_ij - R)/D;			//spatial_distance = r(ij) //R_cons = (r-R)/D
    double f_C_ij = 0.5 - (9/16)*sin(M_PI-R_cons_ij) - (1/16)*sin(((3*M_PI)/2) - R_cons_ij);
    double R_cons_ik = (spatial_distance_ik - R)/D;
    double f_C_ik = 0.5 - (9/16)*sin(M_PI-R_cons_ik) - (1/16)*sin(((3*M_PI)/2) - R_cons_ik);
    double f_R_ij = A*exp(-lam1*spatial_distance_ij);
    double f_A_ij = -B*exp(-lam2*spatial_distance_ij);
    double cons  = pow(h - cos_theta,2);
    double ga_theta_ijk = 1+c4*exp(-c5*cons); 
    double g0_theta_ijk = (c2* cons)/(c3+cons);    
    double g_ijk = c1 + g0_theta_ijk*ga_theta_ijk;
    double spdistance = spatial_distance_ij - spatial_distance_ik;
    double sai_ij = f_C_ik* g_ijk* exp(alpha*pow(spdistance,beta));
    double b_ij = pow((1+pow(sai_ij, eta)), -(1/2*n));
    double V_ij = f_C_ij * ( f_R_ij + b_ij * f_A_ij );
    //double energy;
    return V_ij; 
}
