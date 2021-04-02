/*
 * force.cpp
 *
 *  Created on: Nov 16, 2018
 *      Author: Elmira Birang
 *              Mechanical Engineering PhD candidate
 *              Chair of Applied Mechanics
 *              University of Erlangen-Nürnberg
 */


#include"atom.h"
#include "bond.h"
#include "point.h"
#include "force.h"
#include "energy.h"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

template <int dim>
Force <dim>::Force()
{}

template <int dim>
Force <dim>::~Force()
{}


template <int dim>
Force <dim>::Force(double sigma, double epsilon)

		{
	this -> sigma=sigma;
	this -> epsilon=epsilon;
		}

template <int dim>
Point <dim> Force <dim>::BondForce(Atom <dim> atoma, Atom <dim> atomb)
		{

	// this function calculates interaction force between
	// a pair of atoms using Lennard-Jones potential energy.

	this -> atomi=atoma;
	this ->atomj=atomb;

	Bond <dim> bond;

	Point <dim> interatomic_normal=bond.SpatialBondNormal(atomi,atomj);
	double interatomic_distance=bond.SpatialBondDistance(atomi,atomj);
	double distanceInv=1.0/interatomic_distance;

	Point <dim> force_ij=interatomic_normal*(-4*distanceInv)*
			             (12*pow(distanceInv,12)-6*pow(distanceInv,6));

	return(force_ij);

		}



template <int dim>
Point <dim> Force <dim>::MaterialBondForce(Atom <dim> atoma, Atom <dim> atomb)

		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Bond <dim> bond;

	Point <dim> interatomic_normal=bond.MaterialBondNormal(atomi,atomj);

	double interatomic_stretch=bond.MaterialBondStretch(atomi,atomj);
	double interatomic_stretch_0=bond.InitialBondStretch(atomi,atomj);

	double material_distance=bond.MaterialBondDistance(atomi,atomj);
	double initial_distance=bond.InitialBondDistance(atomi,atomj);

	double sigma0=1./(material_distance*initial_distance);
	double epsilon0=1./(material_distance*initial_distance);

	Point <dim> force_ij= interatomic_normal*4*epsilon0*
			              (13*pow((sigma0*interatomic_stretch*interatomic_stretch_0),12)
	                      -7*pow((sigma0*interatomic_stretch*interatomic_stretch_0),6));

	return (force_ij);

		}


template <int dim>
double Force <dim>::interatomicConfigForceValue(Atom <dim> atoma, Atom <dim> atomb)

		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Bond <dim> bond;

	double interatomic_stretch=bond.MaterialBondStretch(atomi,atomj);
	double interatomic_stretch_0=bond.InitialBondStretch(atomi,atomj);

	double material_distance=bond.MaterialBondDistance(atomi,atomj);
	double initial_distance=bond.InitialBondDistance(atomi,atomj);

	double sigma0=1.0/(material_distance*initial_distance);
	double epsilon0=1.0/(material_distance*initial_distance);

	double term=(sigma0*interatomic_stretch*interatomic_stretch_0);

	double force_ij= 4*epsilon0*(13*pow(term,12)-7*pow(term,6));

	return (force_ij);

		}




template <int dim>
Point <dim> Force <dim>::CriticalMaterialBondForce(Atom <dim> atoma, Atom <dim> atomb,double cutR)

		{


	this -> atomi=atoma;
	this -> atomj=atomb;

	Bond <dim> bond;

	Point <dim> interatomic_normal=bond.MaterialBondNormal(atomi,atomj);

	double interatomic_stretch=bond.MaterialBondStretch(atomi,atomj);
	double interatomic_stretch_0=bond.InitialBondStretch(atomi,atomj);

	double material_distance=bond.MaterialBondDistance(atomi,atomj);
	double initial_distance=bond.InitialBondDistance(atomi,atomj);

	double stretch_critical=98.0/26.0;


	double sigma0=1/(material_distance*initial_distance);
	double epsilon0=1/(material_distance*initial_distance);

//	Point <dim> force_ij= interatomic_normal*4*epsilon0*
//			              (13*pow((sigma0*(stretch_critical)*interatomic_stretch_0),12)
//	                      -7*pow((sigma0*(stretch_critical)*interatomic_stretch_0),6));

	Point <dim> force_ij=interatomic_normal*(stretch_critical)*epsilon0;

	return (force_ij);

		}



template <int dim>
Point <dim> Force <dim>::SpatialBondForce(Atom <dim> atoma, Atom <dim> atomb)

		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	Bond <dim> bond;

	Point <dim> interatomic_normal=bond.SpatialBondNormal(atomi,atomj);
	double interatomic_stretch=bond.SpatialBondStretch(atomi,atomj);

	Point <dim> force_ij= interatomic_normal*(4*(epsilon/interatomic_stretch))*
			              (12*pow((sigma/interatomic_stretch),12)-6*pow((sigma/interatomic_stretch),6));

	return (force_ij);

		}

template <int dim>
Point <dim> Force <dim>:: ConfigurationalForce(Atom<dim> *atoma)

{

	this->atomi=*atoma;
	Point <dim> config_force;

	config_force.SetXCoord(0.0);
	config_force.SetYCoord(0.0);
	config_force.SetZCoord(0.0);

	vector < Atom<dim>* > neighbors=atomi.Neighbor(); //defined in atoms class
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();
	typedef typename vector < Atom<dim>* >::iterator Nghbr;

	for(Nghbr neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		Point <dim> config_interaction_force;

		config_interaction_force=MaterialBondForce(atomi,*neigh);
		config_force=config_force+config_interaction_force;

	}

	for(Nghbr neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		Point <dim> config_interaction_force;

		config_interaction_force=MaterialBondForce(atomi,*neigh);
		config_force=config_force+config_interaction_force;

	}

	return(config_force);

}


template <int dim>
double Force <dim>:: EnergyRelease(Atom<dim> *atoma, double delta_x)

{

	this->atomi=*atoma;
	double config_force;


	vector < Atom<dim>* > neighbors=atomi.Neighbor(); //defined in atoms class
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();
	typedef typename vector < Atom<dim>* >::iterator Nghbr;

	for(Nghbr neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		Point <dim> config_interaction_force(0.,0.,0.);

		config_interaction_force=MaterialBondForce(atomi,*neigh);

		config_force=config_force+(config_interaction_force.GetXCoord()*delta_x);

	}

	for(Nghbr neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		Point <dim> config_interaction_force(0.,0.,0.);

		config_interaction_force=MaterialBondForce(atomi,*neigh);

		config_force=config_force+(config_interaction_force.GetXCoord()*delta_x);

	}

	return(config_force);

}

//calculate the interaction force for each atom with neighbors.
//The neighbor list is in a way that atoms only interact with atoms that have larger id number.
//example: atom0 interacts with atom1,atom2,.. , atom1 interacts with atom2 and atom3

template <int dim>
Point <dim> Force <dim>:: ResultantForce(Atom <dim> *atoma)
{

    this -> atomi=*atoma;

    vector < Atom <dim>* > bond_neighbors;
    bond_neighbors=atomi.BondNeighbor();

    vector < Atom <dim>* > neighbors;
    neighbors=atomi.Neighbor();

    typedef typename vector < Atom <dim>* >::iterator At;

    Point <dim> resultant_force;

    resultant_force.SetXCoord(0.0);
    resultant_force.SetYCoord(0.0);
    resultant_force.SetZCoord(0.0);


    for(At neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
    {

    	 Atom <dim> *neigh;
    	 neigh=*neighbor;

         Point <dim> interatomic_force=BondForce(atomi,*neigh);

         resultant_force=resultant_force+interatomic_force;

    }

//    Point <dim> resultant_force= atoma -> GetForce();

    for(At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
    {

    	 Atom <dim> *neigh;
    	 neigh=*neighbor;

         Point <dim> interatomic_force=BondForce(atomi,*neigh);

//    	 Point <dim> neigh_force_old= neigh->GetForce() ;
//
//    	 neigh -> SetForce((neigh_force_old-interatomic_force)*-1);

         resultant_force=resultant_force+interatomic_force;

    }
//    atoma -> SetForce(resultant_force);

    return (resultant_force);


}

template <int dim>
Point <dim> Force<dim>:: ResultantConfigForce(Atom <dim>* atoma)
{

	this->atomi=*atoma;


	Bond <dim> bond;

	double material_distance=bond.MaterialBondDistance(atomi,atomj);
	double sigma0=1/material_distance;
	double epsilon0=1/material_distance;

	Energy <dim> energy(sigma0,epsilon0);

	vector < Atom<dim>* > neighbors=atomi.Neighbor();
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();
	Point <dim> config_force;

	config_force.SetXCoord(0.0);
	config_force.SetYCoord(0.0);
	config_force.SetZCoord(0.0);

	typedef typename vector < Atom<dim>* >::iterator Nghbr;

	for(Nghbr neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		double deriv_interatomic_config_eng=DerivMaterialEnergy(atomi,*neigh);

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*neigh);
		Point <dim> spatial_normal=bond.SpatialBondNormal(atomi,*neigh);
		double material_energy=energy.ConfigInteratomicEnergy(atomi,*neigh);

		Point <dim> config_interaction_force=spatial_normal*(deriv_interatomic_config_eng*material_bond_stretch
				                             +material_energy)*-1;

		config_force=config_force+config_interaction_force;

	}

	for(Nghbr neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
	{

   	    Atom <dim> *neigh;
   	    neigh=*neighbor;

		double deriv_interatomic_config_eng=DerivMaterialEnergy(atomi,*neigh);

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*neigh);
		Point <dim> spatial_normal=bond.SpatialBondNormal(atomi,*neigh);
		double material_energy=energy.ConfigInteratomicEnergy(atomi,*neigh);

		Point <dim> config_interaction_force=spatial_normal*(deriv_interatomic_config_eng*material_bond_stretch
				                             +material_energy)*-1;

		config_force=config_force+config_interaction_force;

	}

	return(config_force);


}


template <int dim>
double Force <dim>::DerivMaterialEnergy(Atom <dim> atoma, Atom <dim> atomb)

		{

	this -> atomi=atoma;
	this -> atomj=atomb;

	double material_distance=bond.MaterialBondDistance(atomi,atomj);
	double sigma0=1/material_distance;
	double epsilon0=1/material_distance;

	Bond <dim> bond;

	double interatomic_stretch=bond.MaterialBondStretch(atomi,atomj);

	double force_ij= 4*epsilon0*(13*pow((sigma0*interatomic_stretch),12)-7*pow((sigma0*interatomic_stretch),6));

	return (force_ij);

		}


double repulsionControlDerivative (double r_sn, double interatomic_distance)
{

	double repulsion_control_deriv= -1*4*pow ((r_sn-interatomic_distance),3) ;

	return (repulsion_control_deriv);

}



double Psi(double r_cut, double h, double interatomic_distance)
{

	double psi;

	if (interatomic_distance < r_cut)
	{

	double inv_h=1.0/h;
    double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );
    double denominator=1.0/(1+psi_param_fourth);

    psi=psi_param_fourth*denominator;

	}

	else
	{
		double psi=0.0;
	}



    return (psi);

}

double psiDerivative(double r_cut, double h, double interatomic_distance)
{

//	double inv_h=1.0/h;
//	double psi_param_third=pow( ((interatomic_distance-r_cut)*inv_h),3 );
//	double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );
//
//	double denominator=1.0/pow((1+ psi_param_fourth),2);
//
//	double psi_deriv= (4*psi_param_third)*denominator;

	double psi_deriv;

	if (interatomic_distance < r_cut)
	{

	double nominator=4*pow(h,4)*pow((interatomic_distance-r_cut),3);

	double denominator=1.0/(pow((pow((interatomic_distance-r_cut),4)+pow(h,4)),2));

	psi_deriv=nominator*denominator;
	}

	else
	{
		psi_deriv=0.0;

	}



	return (psi_deriv);

}

double calculateM(double r_0, double interatomic_distance, double alpha)
{

	double M=exp(-1*2*alpha*(interatomic_distance-r_0))-(2*exp(-1*alpha*((interatomic_distance-r_0))));
	return (M);

}

double derivM(double alpha, double r_0, double interatomic_distance)
{

	double deriv_M=-1*2*alpha*exp(-1*2*alpha* (interatomic_distance-r_0)) + 2*alpha*exp(-1*alpha*((interatomic_distance-r_0))) ;

	return (deriv_M);

}



template <int dim>
Point <dim> Force <dim> :: derivativeMorseFunction(Atom <dim>* atoma, double alpha1, double alpha2,
                                                   double r_01, double r_02, double E1, double E2, double delta,
                                                   double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                                   double s1, double s2, double s3)

{

	Point <dim> atoma_force;

	atoma_force.SetXCoord(0.0);
	atoma_force.SetYCoord(0.0);
	atoma_force.SetZCoord(0.0);

	this->atomi=*atoma;

	typedef typename vector < Atom <dim> * > ::iterator neigh;

	vector < Atom <dim> * > neighbors=atomi.Neighbor();
	vector < Atom <dim> * > bond_neighbors=atomi.BondNeighbor();

	for (neigh neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		double interatomic_bond=bond.SpatialBondDistance(atomi,*neigh);

		double part1=0.0;
		double part2_1=0;
		double part2_2=0.0;
		double part2_3=0.0;
		double morse_deriv_ij=0;
		Point <dim> morse_force_ij;
		morse_force_ij.SetXCoord(0);
		morse_force_ij.SetYCoord(0);
		morse_force_ij.SetZCoord(0);

	    if (interatomic_bond < r_cut)
	    {

		     part1=(E1*derivM(alpha1, r_01, interatomic_bond)+E2*derivM(alpha2, r_02, interatomic_bond))*Psi(r_cut, h, interatomic_bond)
		    	   + (E1*calculateM(r_01, interatomic_bond, alpha1)+E2*calculateM(r_02, interatomic_bond, alpha2)+delta)*psiDerivative(r_cut, h, interatomic_bond);

	    }



	    if ( r_s1 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_1= s1*repulsionControlDerivative (r_s1, interatomic_bond);


	    }


	    if ( r_s2 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_2=s2*repulsionControlDerivative (r_s2, interatomic_bond);

	    }


	    if ( r_s3 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_3=s3*repulsionControlDerivative (r_s3, interatomic_bond);

	    }


	    morse_deriv_ij=part1-part2_1-part2_2-part2_3;

	    Point <dim> bond_direction=bond.SpatialBondNormal(atomi,*neigh);

	    morse_force_ij=bond_direction*morse_deriv_ij;

	    atoma_force = atoma_force + morse_force_ij;

	}


	for (neigh bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	{

		Atom <dim> *bond_neigh;
		bond_neigh=*bond_neighbor;

		double interatomic_bond=bond.SpatialBondDistance(atomi,*bond_neigh);

		double part1=0.0;
		double part2_1=0;
		double part2_2=0.0;
		double part2_3=0.0;
		double morse_deriv_ij=0;

		Point <dim> morse_force_ij;
		morse_force_ij.SetXCoord(0);
		morse_force_ij.SetYCoord(0);
		morse_force_ij.SetZCoord(0);

	    if (interatomic_bond < r_cut)
	    {

		     part1=(E1*derivM(alpha1, r_01, interatomic_bond)+E2*derivM(alpha2, r_02, interatomic_bond))*Psi(r_cut, h, interatomic_bond)
		    	    + (E1*calculateM(r_01, interatomic_bond, alpha1)+E2*calculateM(r_02, interatomic_bond, alpha2)+delta)*psiDerivative(r_cut, h, interatomic_bond);

	    }



	    if ( r_s1 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_1= s1*repulsionControlDerivative (r_s1, interatomic_bond);


	    }


	    if ( r_s2 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_2=s2*repulsionControlDerivative (r_s2, interatomic_bond);

	    }


	    if ( r_s3 > interatomic_bond )
	    {
//	    	cout << "here" << endl;

	    	part2_3=s3*repulsionControlDerivative (r_s3, interatomic_bond);

	    }


	    morse_deriv_ij=part1-part2_1-part2_2-part2_3;

	    Point <dim> bond_direction=bond.SpatialBondNormal(atomi,*bond_neigh);

	    morse_force_ij=bond_direction*morse_deriv_ij;

	    atoma_force = atoma_force + morse_force_ij;

	}


	return(atoma_force);


}


template <int dim>
double Force <dim> ::electronDensity (Atom<dim>* atoma, double a, double beta1,
                                     double beta2, double r_03, double r_04,
                                     double h, double r_cut)
{

	Point <dim> atoma_force;
	this->atomi=*atoma;
	double atom_electron_density=0;

	double inv_h=1.0/h;


	typedef typename vector < Atom <dim> * > ::iterator neigh;

	vector < Atom <dim> * > neighbors=atomi.Neighbor();
	vector < Atom <dim> * > bond_neighbors=atomi.BondNeighbor();


	for (neigh neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;
		double psi=0;

		double interatomic_distance= bond.SpatialBondDistance(atomi, *neigh);

	    if (interatomic_distance < r_cut)
	    {

	       double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );

	       double denominator=1.0/(1+psi_param_fourth);

	       psi=psi_param_fourth*denominator;

	    }

	    else
	    {
	    	psi=0.0;

	    }

		double electron_density_ij= (a * exp (-1*beta1 * pow((interatomic_distance-r_03),2)) + exp (-1*beta2 *(interatomic_distance-r_04)) ) * psi;
		atom_electron_density+=electron_density_ij;

	}


	for (neigh bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	{

		Atom <dim> *bond_neigh;
		bond_neigh=*bond_neighbor;
		double psi=0;

		double interatomic_distance= bond.SpatialBondDistance(atomi, *bond_neigh);

	    if (interatomic_distance < r_cut)
	    {

	       double psi_param_fourth=pow( ((interatomic_distance-r_cut)*inv_h),4 );

	       double denominator=1.0/(1+psi_param_fourth);

	       psi=psi_param_fourth*denominator;

	    }

	    else
	    {
	    	psi=0.0;

	    }

		double electron_density_ij= (a * exp (-1*beta1 * pow((interatomic_distance-r_03),2)) + exp (-1*beta2 * (interatomic_distance-r_04)) ) * psi;
		atom_electron_density+=electron_density_ij;

	}

	return (atom_electron_density);


}

double derivativeElectronDensity(double interatomic_distance, double a, double beta1,
                                 double beta2, double r_03, double r_04,
                                 double h, double r_cut)
{

	double deriv_electron_density=((-2*a*beta1*(interatomic_distance-r_03))*exp(-1*beta1*pow ((interatomic_distance-r_03),2))
			                      -beta2*exp(-1*beta2*(interatomic_distance-r_04)) )*Psi(r_cut, h, interatomic_distance) +
			                      (a*exp(-1*beta1*pow ((interatomic_distance-r_03),2))+exp(-1*beta2*(interatomic_distance-r_04)))*psiDerivative(r_cut, h, interatomic_distance);


	return(deriv_electron_density);

}

//F'(rho)
double derivativeEmbeddingEnergy (double electron_density,
                                  double F0, double F2,
                                  double q1, double q2,
								  double q3, double q4,
								  double Q1, double Q2,
                                  double h, double r_cut )
{

	double derivative_embedding_energy=0.0;

	if (electron_density < 1)
	{

	derivative_embedding_energy=F2*(electron_density-1)
			                    + 3*q1*pow ((electron_density-1),2)
	                            + 4*q2*pow ((electron_density-1),3)
	                            + 5*q3*pow ((electron_density-1),4)
	                            + 6*q4*pow ((electron_density-1),5);

	}

	if (electron_density > 1)

	{

        double denominator= 1.0/(pow((1+Q2*pow((electron_density-1),3)),2));

        double nominator= (F2*(electron_density-1)+3*q1*pow((electron_density-1),2)+4*Q1*pow((electron_density-1),3)) * (1+Q2*pow((electron_density-1),3))
        	              - (F0+0.5*F2*pow((electron_density-1),2)+q1*pow((electron_density-1),3)+Q1*pow((electron_density-1),4))*(3*Q2*pow((electron_density-1),2));

		derivative_embedding_energy=nominator*denominator;

	}

	return (derivative_embedding_energy);

}




template <int dim>
Point <dim> Force <dim> :: embeddingForce(Atom <dim> * atoma, double a, double beta1,
                                          double beta2, double r_03, double r_04, double F0, double F2,
                                          double q1, double q2, double q3, double q4, double Q1, double Q2,
		                                  double h, double r_cut)

{

	Force <dim> force;

	Point <dim> atoma_force;

	atoma_force.SetXCoord(0.0);
	atoma_force.SetYCoord(0.0);
	atoma_force.SetZCoord(0.0);

	this->atomi=*atoma;
	double atom_electron_density=force.electronDensity(atoma, a, beta1,
                                                 beta2, r_03, r_04,
                                                 h, r_cut);


	double atom_deriv_embedding_eng=derivativeEmbeddingEnergy(atom_electron_density,
                                                     F0, F2,
                                                     q1, q2,
			                                         q3, q4,
			                                         Q1, Q2,
                                                     h, r_cut);


	double inv_h=1.0/h;
	double psi=0.0;

	typedef typename vector < Atom <dim> * > ::iterator neigh;

	vector < Atom <dim> * > neighbors=atomi.Neighbor();


	for (neigh neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> Embedding_Force;

		Embedding_Force.SetXCoord(0.0);
		Embedding_Force.SetYCoord(0.0);
		Embedding_Force.SetZCoord(0.0);
		double interatomic_distance =0;
		double neighbor_electron_density=0;
		double neighbor_deriv_embedding_eng=0;
		double deriv_electron_density=0;

		interatomic_distance = bond.SpatialBondDistance(atomi, *neigh) ;

		neighbor_electron_density=force.electronDensity (neigh, a, beta1,
	                                                            beta2, r_03, r_04,
	                                                            h, r_cut);


		neighbor_deriv_embedding_eng=derivativeEmbeddingEnergy(neighbor_electron_density,
		                                                              F0, F2,
		                                                              q1, q2,
					                                                  q3, q4,
					                                                  Q1, Q2,
		                                                              h, r_cut);


		deriv_electron_density= derivativeElectronDensity(interatomic_distance, a, beta1,
		                                                         beta2, r_03, r_04,
		                                                         h, r_cut);


		Point <dim> bond_direction=bond.SpatialBondNormal(atomi,*neigh);


		Embedding_Force=bond_direction*(neighbor_deriv_embedding_eng+atom_deriv_embedding_eng)*deriv_electron_density;

		atoma_force=atoma_force+Embedding_Force;

	}




	vector < Atom <dim> * > bond_neighbors=atomi.BondNeighbor();

	for (neigh neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> Embedding_Force;

		Embedding_Force.SetXCoord(0.0);
		Embedding_Force.SetYCoord(0.0);
		Embedding_Force.SetZCoord(0.0);
		double interatomic_distance =0;
		double neighbor_electron_density=0;
		double neighbor_deriv_embedding_eng=0;
		double deriv_electron_density=0;

		interatomic_distance = bond.SpatialBondDistance(atomi, *neigh) ;

		neighbor_electron_density=force.electronDensity (neigh, a, beta1,
	                                                     beta2, r_03, r_04,
	                                                     h, r_cut);



		neighbor_deriv_embedding_eng=derivativeEmbeddingEnergy(neighbor_electron_density,
		                                                       F0, F2,
		                                                       q1, q2,
					                                           q3, q4,
					                                           Q1, Q2,
		                                                       h, r_cut);


		deriv_electron_density= derivativeElectronDensity(interatomic_distance, a, beta1,
		                                                         beta2, r_03, r_04,
		                                                         h, r_cut);


		Point <dim> bond_direction=bond.SpatialBondNormal(atomi,*neigh);


		Embedding_Force=bond_direction*(neighbor_deriv_embedding_eng+atom_deriv_embedding_eng)*deriv_electron_density;


		atoma_force=atoma_force+Embedding_Force;

	}

	return(atoma_force);

}


template <int dim>
Point <dim> Force <dim>:: configForceEamPair(Atom <dim>* atoma, double alpha1, double alpha2,
                                             double r_01, double r_02, double E1, double E2, double delta,
                                             double r_cut, double h, double r_s1, double r_s2, double r_s3,
                                             double s1, double s2, double s3)
{

	Point <dim> atoma_force;

	atoma_force.SetXCoord(0.0);
	atoma_force.SetYCoord(0.0);
	atoma_force.SetZCoord(0.0);

	this->atomi=*atoma;

	typedef typename vector < Atom <dim> * > ::iterator neigh;

	vector < Atom <dim> * > neighbors=atomi.Neighbor();
	vector < Atom <dim> * > bond_neighbors=atomi.BondNeighbor();

	for (neigh neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*neigh);

		double inv_material_bond_stretch=1.0/material_bond_stretch;

		double interatomic_bond=bond.SpatialBondDistance(atomi,*neigh);

		double material_bond_distance=bond.MaterialBondDistance(atomi,*neigh);

		double alpha1_mod=alpha1*material_bond_distance;
		double alpha2_mod=material_bond_distance*alpha2;
		double r_01_mod=1.0/material_bond_distance*r_01;
		double r_02_mod=1.0/material_bond_distance*r_02;
		double r_cut_mod=1.0/material_bond_distance*r_cut;
		double h_mod=1.0/material_bond_distance*h;
		double r_s1_mod=1.0/material_bond_distance*r_s1;
		double r_s2_mod=1.0/material_bond_distance*r_s2;
		double r_s3_mod=1.0/material_bond_distance*r_s3;
		double s1_mod=pow(material_bond_distance,4)*s1;
		double s2_mod=pow(material_bond_distance,4)*s2;
		double s3_mod=pow(material_bond_distance,4)*s3;

		double part1=0.0;
		double part2_1=0;
		double part2_2=0.0;
		double part2_3=0.0;
		double morse_deriv_ij=0;
		Point <dim> morse_force_ij;
		morse_force_ij.SetXCoord(0);
		morse_force_ij.SetYCoord(0);
		morse_force_ij.SetZCoord(0);

	    if (interatomic_bond < r_cut_mod)
	    {

		     part1=(E1*derivM(alpha1_mod, r_01_mod, inv_material_bond_stretch)+E2*derivM(alpha2_mod, r_02_mod, inv_material_bond_stretch))*Psi(r_cut_mod, h_mod, inv_material_bond_stretch)
		    	   + (E1*calculateM(r_01_mod, inv_material_bond_stretch, alpha1_mod)+
		    	   E2*calculateM(r_02_mod, inv_material_bond_stretch, alpha2_mod)+delta)*psiDerivative(r_cut_mod, h_mod, inv_material_bond_stretch);

	    }



	    if ( r_s1_mod > interatomic_bond )
	    {

	    	part2_1= s1_mod*repulsionControlDerivative (r_s1_mod, inv_material_bond_stretch);


	    }


	    if ( r_s2_mod > interatomic_bond )
	    {

	    	part2_2=s2_mod*repulsionControlDerivative (r_s2_mod, inv_material_bond_stretch);

	    }


	    if ( r_s3_mod > interatomic_bond )
	    {

	    	part2_3=s3_mod*repulsionControlDerivative (r_s3_mod, inv_material_bond_stretch);

	    }


	    morse_deriv_ij=part1-material_bond_distance*(part2_1+part2_2+part2_3);

	    Point <dim> bond_direction=bond.SpatialBondNormal(atomi,*neigh);

	    morse_force_ij=bond_direction*morse_deriv_ij;

	    atoma_force = atoma_force + morse_force_ij;

	}


	for (neigh bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	{

		Atom <dim> *bond_neigh;
		bond_neigh=*bond_neighbor;

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*bond_neigh);

		double inv_material_bond_stretch=1.0/material_bond_stretch;

		double interatomic_bond=bond.SpatialBondDistance(atomi,*bond_neigh);

		double material_bond_distance=bond.MaterialBondDistance(atomi,*bond_neigh) ;

		double alpha1_mod=alpha1*material_bond_distance;
		double alpha2_mod=material_bond_distance*alpha2;
		double r_01_mod=1.0/material_bond_distance*r_01;
		double r_02_mod=1.0/material_bond_distance*r_02;
		double r_cut_mod=1.0/material_bond_distance*r_cut;
		double h_mod=1.0/material_bond_distance*h;
		double r_s1_mod=1.0/material_bond_distance*r_s1;
		double r_s2_mod=1.0/material_bond_distance*r_s2;
		double r_s3_mod=1.0/material_bond_distance*r_s3;
		double s1_mod=pow(material_bond_distance,4)*s1;
		double s2_mod=pow(material_bond_distance,4)*s2;
		double s3_mod=pow(material_bond_distance,4)*s3;

		double part1=0.0;
		double part2_1=0;
		double part2_2=0.0;
		double part2_3=0.0;
		double morse_deriv_ij=0;

		Point <dim> morse_force_ij;
		morse_force_ij.SetXCoord(0);
		morse_force_ij.SetYCoord(0);
		morse_force_ij.SetZCoord(0);

	    if (interatomic_bond < r_cut_mod)
	    {

		     part1=(E1*derivM(alpha1_mod, r_01_mod, inv_material_bond_stretch)+E2*derivM(alpha2_mod, r_02_mod, inv_material_bond_stretch))*Psi(r_cut_mod, h_mod, inv_material_bond_stretch)
		    	    + (E1*calculateM(r_01_mod, inv_material_bond_stretch, alpha1_mod)+
		    	    E2*calculateM(r_02_mod, inv_material_bond_stretch, alpha2_mod)+delta)*psiDerivative(r_cut_mod, h_mod, inv_material_bond_stretch);

	    }



	    if ( r_s1_mod > interatomic_bond )
	    {

	    	part2_1= s1_mod*repulsionControlDerivative (r_s1_mod, inv_material_bond_stretch);


	    }


	    if ( r_s2_mod > interatomic_bond )
	    {

	    	part2_2=s2_mod*repulsionControlDerivative (r_s2_mod, inv_material_bond_stretch);

	    }


	    if ( r_s3_mod > interatomic_bond )
	    {

	    	part2_3=s3_mod*repulsionControlDerivative (r_s3_mod, inv_material_bond_stretch);

	    }


	    morse_deriv_ij=part1-material_bond_distance*(part2_1+part2_2+part2_3);

	    Point <dim> bond_direction=bond.MaterialBondNormal(atomi,*bond_neigh);

	    morse_force_ij=bond_direction*morse_deriv_ij;

	    atoma_force = atoma_force + morse_force_ij;

	}


	return(atoma_force);


}

template <int dim>
Point <dim> Force <dim> :: configForceEamEmbedding(Atom <dim> * atoma, double a, double beta1,
                                                   double beta2, double r_03, double r_04, double F0, double F2,
                                                   double q1, double q2, double q3, double q4, double Q1, double Q2,
                                                   double h, double r_cut)
{

	Force <dim> force;

	Point <dim> atoma_force;

	atoma_force.SetXCoord(0.0);
	atoma_force.SetYCoord(0.0);
	atoma_force.SetZCoord(0.0);

	this->atomi=*atoma;
	double atom_electron_density=force.electronDensity(atoma, a, beta1,
                                                       beta2, r_03, r_04,
                                                       h, r_cut);

	double atom_deriv_embedding_eng=derivativeEmbeddingEnergy(atom_electron_density,
                                                              F0, F2,
                                                              q1, q2,
			                                                  q3, q4,
			                                                  Q1, Q2,
                                                              h, r_cut);


	double inv_h=1.0/h;
	double psi=0.0;

	typedef typename vector < Atom <dim> * > ::iterator neigh;

	vector < Atom <dim> * > neighbors=atomi.Neighbor();


	for (neigh neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> Embedding_Force;

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*neigh);

		double inv_material_bond_stretch=1.0/material_bond_stretch;

		Embedding_Force.SetXCoord(0.0);
		Embedding_Force.SetYCoord(0.0);
		Embedding_Force.SetZCoord(0.0);
		double interatomic_distance =0;
		double neighbor_electron_density=0;
		double neighbor_deriv_embedding_eng=0;
		double deriv_electron_density=0;

		double material_bond_distance=bond.MaterialBondDistance(atomi,*neigh) ;

		double beta1_mod=pow(material_bond_distance,2)*beta1;
		double beta2_mod=material_bond_distance*beta2;

		double r_03_mod=(1.0/material_bond_distance)*r_03;
		double r_04_mod=(1.0/material_bond_distance)*r_04;

		double h_mod=(1.0/material_bond_distance)*h;
		double r_cut_mod=(1.0/material_bond_distance)*r_cut;

		interatomic_distance = bond.SpatialBondDistance(atomi, *neigh) ;

		neighbor_electron_density=force.electronDensity (neigh, a, beta1_mod,
	                                                            beta2_mod, r_03_mod, r_04_mod,
	                                                            h_mod, r_cut_mod);


		neighbor_deriv_embedding_eng=derivativeEmbeddingEnergy(neighbor_electron_density,
		                                                       F0, F2,
		                                                       q1, q2,
					                                           q3, q4,
					                                           Q1, Q2,
		                                                       h_mod, r_cut_mod);


		deriv_electron_density= derivativeElectronDensity(inv_material_bond_stretch, a, beta1_mod,
		                                                  beta2_mod, r_03_mod, r_04_mod,
		                                                  h_mod, r_cut_mod);


		Point <dim> bond_direction=bond.MaterialBondNormal(atomi,*neigh);


		Embedding_Force=bond_direction*(neighbor_deriv_embedding_eng+atom_deriv_embedding_eng)*deriv_electron_density;

		atoma_force=atoma_force+Embedding_Force;

	}




	vector < Atom <dim> * > bond_neighbors=atomi.BondNeighbor();

	for (neigh neighbor=bond_neighbors.begin(); neighbor!=bond_neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> Embedding_Force;

		double material_bond_stretch=bond.MaterialBondStretch(atomi,*neigh);

		double inv_material_bond_stretch=1.0/material_bond_stretch;

		Embedding_Force.SetXCoord(0.0);
		Embedding_Force.SetYCoord(0.0);
		Embedding_Force.SetZCoord(0.0);
		double interatomic_distance =0;
		double neighbor_electron_density=0;
		double neighbor_deriv_embedding_eng=0;
		double deriv_electron_density=0;

		double material_bond_distance=bond.MaterialBondDistance(atomi,*neigh) ;

		double beta1_mod=pow(material_bond_distance,2)*beta1;
		double beta2_mod=material_bond_distance*beta2;

		double r_03_mod=(1.0/material_bond_distance)*r_03;
		double r_04_mod=(1.0/material_bond_distance)*r_04;

		double h_mod=(1.0/material_bond_distance)*h;
		double r_cut_mod=(1.0/material_bond_distance)*r_cut;

		interatomic_distance = bond.SpatialBondDistance(atomi, *neigh) ;

		neighbor_electron_density=force.electronDensity (neigh, a, beta1_mod,
	                                                     beta2_mod, r_03_mod, r_04_mod,
	                                                     h_mod, r_cut_mod);



		neighbor_deriv_embedding_eng=derivativeEmbeddingEnergy(neighbor_electron_density,
		                                                       F0, F2,
		                                                       q1, q2,
					                                           q3, q4,
					                                           Q1, Q2,
		                                                       h_mod, r_cut_mod);


		deriv_electron_density= derivativeElectronDensity(inv_material_bond_stretch, a, beta1_mod,
		                                                         beta2_mod, r_03_mod, r_04_mod,
		                                                         h_mod, r_cut_mod);


		Point <dim> bond_direction=bond.MaterialBondNormal(atomi,*neigh);

		Embedding_Force=bond_direction*(neighbor_deriv_embedding_eng+atom_deriv_embedding_eng)*deriv_electron_density;

		atoma_force=atoma_force+Embedding_Force;

	}

	return(atoma_force);

}

template <int dim>
Point <dim> Force <dim>::StillingerWeberThreeBodyForceH2H3(Atom <dim> &atom_Beta, Atom <dim> &atom_Alpha,
		                                               Atom <dim> &atom_Gamma, double sigma_AlphaBeta,
                                                       double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
													   double gamma, double lambda_AlphaBetaGamma
                                                       ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{

	this -> atomi=atom_Alpha;
	this -> atomj=atom_Beta;
	this -> atomk=atom_Gamma;

	Point <dim> spatialPositionAlpha(0.,0.,0.);
	Point <dim> spatialPositionBeta(0.,0.,0.);
	Point <dim> spatialPositionGamma(0.,0.,0.);

	spatialPositionAlpha=atomi.GetSpatialPosition();
	spatialPositionBeta=atomj.GetSpatialPosition();
	spatialPositionGamma=atomk.GetSpatialPosition();

	Point <dim> interatomic_distance_AlphaBeta(0.,0.,0.);
	Point <dim> interatomic_distance_AlphaGamma(0.,0.,0.);

	interatomic_distance_AlphaBeta=spatialPositionBeta-spatialPositionAlpha;
	interatomic_distance_AlphaGamma=spatialPositionGamma-spatialPositionAlpha;
	
	double distance_AlphaBeta=sqrt((pow(interatomic_distance_AlphaBeta.GetXCoord(),2) +
			                   pow(interatomic_distance_AlphaBeta.GetYCoord(),2) +
			                   pow(interatomic_distance_AlphaBeta.GetZCoord(),2)));

	double distance_AlphaGamma=sqrt((pow(interatomic_distance_AlphaGamma.GetXCoord(),2) +
			                    pow(interatomic_distance_AlphaGamma.GetYCoord(),2) +
			                    pow(interatomic_distance_AlphaGamma.GetZCoord(),2)));

	double distance_BetaGamma=bond.SpatialBondDistance(atomj,atomk);

	Point <dim> SWThreeBodyForce;

	SWThreeBodyForce.SetXCoord(0.0);
	SWThreeBodyForce.SetYCoord(0.0);
	SWThreeBodyForce.SetZCoord(0.0);

	cout << distance_AlphaBeta<< "\n";

	Point <dim> normalVectorAlphaBeta(0.,0.,0.);
	Point <dim> normalVectorAlphaGamma(0.,0.,0.);

	normalVectorAlphaBeta=bond.SpatialBondNormal(atomj, atomi);
	normalVectorAlphaGamma=bond.SpatialBondNormal(atomk, atomi);
	
	double cosine_teta_BetaAlphaGamma=(pow(distance_AlphaBeta,2)+pow(distance_AlphaGamma,2)-pow(distance_BetaGamma,2))
									   /(2*distance_AlphaBeta*distance_AlphaGamma);
	

	double term1=(gamma*sigma_AlphaBeta)/(distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta));
	double term2=(gamma*sigma_AlphaGamma)/(distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma));
	double terms_sum=term1+term2;

	double exponent=exp(terms_sum);
	
	double difference_cosines=(cosine_teta_BetaAlphaGamma-cosine_teta0)/100;
	
	double difference_cosines_squared=difference_cosines*difference_cosines;
	
	Point <dim> term3(0.,0.,0.);

	term3.SetXCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetXCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta));

	term3.SetYCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta));

	term3.SetZCoord(cosine_teta_BetaAlphaGamma*(normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta)
					-(normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta));

	Point <dim> term4(0.,0.,0.);

	term4.SetXCoord(term1*normalVectorAlphaBeta.GetXCoord());
	term4.SetYCoord(term1*normalVectorAlphaBeta.GetYCoord());
	term4.SetZCoord(term1*normalVectorAlphaBeta.GetZCoord());

	SWThreeBodyForce.SetXCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetXCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetXCoord());

	SWThreeBodyForce.SetYCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetYCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetYCoord());

	SWThreeBodyForce.SetZCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetZCoord()+
								lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines_squared*term4.GetZCoord());


	return(SWThreeBodyForce);


}

template <int dim>
Point <dim> Force <dim>::StillingerWeberThreeBodyForceH1(Atom <dim> *atom_Beta, Atom <dim> *atom_Alpha,
		                                                 Atom <dim> *atom_Gamma, double sigma_AlphaBeta,
                                                         double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
													     double gamma, double lambda_AlphaBetaGamma
                                                        ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)

{

	this -> atomi=*atom_Alpha;
	this -> atomj=*atom_Beta;
	this -> atomk=*atom_Gamma;

	Point <dim> spatialPositionAlpha=(*atom_Alpha).GetSpatialPosition();
	Point <dim> spatialPositionBeta=(*atom_Beta).GetSpatialPosition();
	Point <dim> spatialPositionGamma=(*atom_Gamma).GetSpatialPosition();

	double distance_AlphaBeta=bond.SpatialBondDistance(atomj,atomi);

	double distance_AlphaGamma=bond.SpatialBondDistance(atomk,atomi);
	
	double distance_BetaGamma=bond.SpatialBondDistance(atomk,atomj);

	Point <dim> SWThreeBodyForce;

	SWThreeBodyForce.SetXCoord(0.0);
	SWThreeBodyForce.SetYCoord(0.0);
	SWThreeBodyForce.SetZCoord(0.0);

	if (distance_AlphaBeta < (sigma_AlphaBeta*a_AlphaBeta)
	    && distance_AlphaGamma < (sigma_AlphaGamma*a_AlphaGamma))
	{

		Point <dim> normalVectorAlphaBeta=bond.SpatialBondNormal(atomj, atomi);
		Point <dim> normalVectorAlphaGamma=bond.SpatialBondNormal(atomk, atomi);
		Point <dim> normalVectorBetaGamma=bond.SpatialBondNormal(atomj, atomk);
		
		
		double cosine_teta_BetaAlphaGamma=(pow(distance_AlphaBeta,2)+pow(distance_AlphaGamma,2)-pow(distance_BetaGamma,2))
						                   /(2*distance_AlphaBeta*distance_AlphaGamma);
	
		double term1=(gamma*sigma_AlphaBeta)/(distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta));
		double term2=(gamma*sigma_AlphaGamma)/(distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma));
	
		double exponent=exp(term1)*exp(term2);
		
		float difference_cosines=cosine_teta_BetaAlphaGamma-cosine_teta0;
			
		double difference_cosines_squared=difference_cosines*difference_cosines;
	
		Point <dim> term3(0.,0.,0.);
		
		term3.SetXCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetXCoord()/distance_AlphaBeta
				         -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetXCoord()/distance_AlphaGamma)
				
						+(normalVectorAlphaBeta.GetXCoord()/distance_AlphaGamma+
						  normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta));
	
		term3.SetYCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta
				        -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetYCoord()/distance_AlphaGamma)
				
						+(normalVectorAlphaBeta.GetYCoord()/distance_AlphaGamma+
						normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta));
	
	
		term3.SetZCoord((-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta
				        -cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetZCoord()/distance_AlphaGamma)
				
						+(normalVectorAlphaBeta.GetZCoord()/distance_AlphaGamma+
						normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta));
	
	
		Point <dim> term4(0.,0.,0.);
	
		double prefactor1=0.;
		double prefactor2=0.;
	
		prefactor1=sigma_AlphaBeta/((distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta))*
				                      (distance_AlphaBeta-(a_AlphaBeta*sigma_AlphaBeta)));
	
		prefactor2=sigma_AlphaGamma/((distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma))*
				                       (distance_AlphaGamma-(a_AlphaGamma*sigma_AlphaGamma)));
		
	
		term4.SetXCoord(prefactor1*normalVectorAlphaBeta.GetXCoord()+prefactor2*normalVectorAlphaGamma.GetXCoord());
		term4.SetYCoord(prefactor1*normalVectorAlphaBeta.GetYCoord()+prefactor2*normalVectorAlphaGamma.GetYCoord());
		term4.SetZCoord(prefactor1*normalVectorAlphaBeta.GetZCoord()+prefactor2*normalVectorAlphaGamma.GetZCoord());
		
		
		Energy <dim> energy;
		
		double threebodyEng=energy.stillingerWeberThreeBody(atomj, atomi,
                                                            atomk, sigma_AlphaBeta,
                                                            sigma_AlphaGamma, epsilon_AlphaBetaGamma,
			                                                gamma, lambda_AlphaBetaGamma
                                                            , cosine_teta0, a_AlphaBeta, a_AlphaGamma);
		
	
		 
			
		SWThreeBodyForce.SetXCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetXCoord()+
				                   gamma*threebodyEng*term4.GetXCoord());
	
		SWThreeBodyForce.SetYCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetYCoord()+
				                    gamma*threebodyEng*term4.GetYCoord());
	
		SWThreeBodyForce.SetZCoord(-2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*term3.GetZCoord()+
								   gamma*threebodyEng*term4.GetZCoord());
		
		Point <dim> force_atom_beta (0.,0.,0.);
		
		
		force_atom_beta.SetXCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
                                  *exponent*difference_cosines*
								  normalVectorAlphaBeta.GetXCoord()*(1/distance_AlphaGamma-(cosine_teta_BetaAlphaGamma/distance_AlphaBeta))
	                              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetXCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								  normalVectorBetaGamma.GetXCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
				
		force_atom_beta.SetYCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetYCoord()/distance_AlphaBeta
					              +normalVectorAlphaBeta.GetYCoord()/distance_AlphaGamma)
					              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetYCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetYCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		
		force_atom_beta.SetZCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaBeta.GetZCoord()/distance_AlphaBeta
					              +normalVectorAlphaBeta.GetZCoord()/distance_AlphaGamma)
					              -gamma*threebodyEng*(prefactor1*normalVectorAlphaBeta.GetZCoord())
								  -2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetZCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		

		
		atom_Beta->SetForce(atom_Beta->GetForce()+force_atom_beta);
		
		Point <dim> force_atom_gamma (0.,0.,0.);
		
		force_atom_gamma.SetXCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetXCoord()/distance_AlphaGamma
					              +normalVectorAlphaGamma.GetXCoord()/distance_AlphaBeta)
					              -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetXCoord())
								  +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetXCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		
		force_atom_gamma.SetYCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                               *exponent*difference_cosines*
				                   (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetYCoord()/distance_AlphaGamma
					               +normalVectorAlphaGamma.GetYCoord()/distance_AlphaBeta)
					               -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetYCoord())
								   +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   	normalVectorBetaGamma.GetYCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		
		force_atom_gamma.SetZCoord(2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma
	                              *exponent*difference_cosines*
				                  (-cosine_teta_BetaAlphaGamma*normalVectorAlphaGamma.GetZCoord()/distance_AlphaGamma
					              +normalVectorAlphaGamma.GetZCoord()/distance_AlphaBeta)
					              -gamma*threebodyEng*(prefactor2*normalVectorAlphaGamma.GetZCoord())
								  +2*lambda_AlphaBetaGamma*epsilon_AlphaBetaGamma*exponent*difference_cosines*
								   normalVectorBetaGamma.GetZCoord()*distance_BetaGamma/(distance_AlphaBeta*distance_AlphaGamma));
		
		atom_Gamma->SetForce(atom_Gamma->GetForce()+force_atom_gamma);

		
	}

	return(SWThreeBodyForce);

}

template <int dim>
Point <dim> Force <dim>:: StillingerWeberTwoBodyForce(Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta,
                                                      double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                      double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                                      double q_AlphaBeta, double a_AlphaBeta)
{

	this -> atomi=atom_Alpha;
	this -> atomj=atom_Beta;
	
	double distance=bond.SpatialBondDistance(atomi,atomj);

	double term1=sigma_AlphaBeta/pow((distance-(a_AlphaBeta*sigma_AlphaBeta)),2);
	
	Energy <dim> energy;

	double energy_ij=energy.stillingerWeberTwoBody(atomi, atomj, sigma_AlphaBeta, epsilon_AlphaBeta,
                                                   A_AlphaBeta, B_AlphaBeta, p_AlphaBeta, 
												    q_AlphaBeta, a_AlphaBeta);

    Point <dim> normalVectorAlphaBeta=bond.SpatialBondNormal(atomj, atomi);

    double term2_nominator=0.;
    double term2_denominator=0;

    term2_nominator=(p_AlphaBeta*B_AlphaBeta*pow((sigma_AlphaBeta/distance),(p_AlphaBeta)))/distance-
    		        (q_AlphaBeta*pow((sigma_AlphaBeta/distance),(q_AlphaBeta)))/distance;

    term2_denominator=B_AlphaBeta*pow((sigma_AlphaBeta/distance),p_AlphaBeta)
                      -pow((sigma_AlphaBeta/distance),q_AlphaBeta);

    double term2=term2_nominator/term2_denominator;

    Point <dim> SWTwoBodyForce;

    SWTwoBodyForce.SetXCoord(0.);
    SWTwoBodyForce.SetYCoord(0.);
    SWTwoBodyForce.SetZCoord(0.);
    
	if (distance < (sigma_AlphaBeta*a_AlphaBeta))
	{

		SWTwoBodyForce.SetXCoord(energy_ij*normalVectorAlphaBeta.GetXCoord()*(term1+term2));
		SWTwoBodyForce.SetYCoord(energy_ij*normalVectorAlphaBeta.GetYCoord()*(term1+term2));
		SWTwoBodyForce.SetZCoord(energy_ij*normalVectorAlphaBeta.GetZCoord()*(term1+term2));
	}

	return(SWTwoBodyForce);

}

template <int dim>
Point <dim> Force <dim>::ResultantSWTwoBodyForce(Atom <dim> *atoma, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                 double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                                 double q_AlphaBeta, double a_AlphaBeta)
{

	this->atomi=*atoma;

	Point <dim> atoma_force;
	
	atoma_force.SetXCoord(0.);
	atoma_force.SetYCoord(0.);
	atoma_force.SetZCoord(0.);

	vector < Atom <dim>* > neighbors=atomi.Neighbor();
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();

	typedef typename vector < Atom <dim>* >::const_iterator Neigh;

	for (Neigh neighbor=neighbors.begin();  neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> SW2BodyForce(0.0,0.0,0.0);

		SW2BodyForce=StillingerWeberTwoBodyForce((*atoma), (*neigh),
		                                         sigma_AlphaBeta, epsilon_AlphaBeta,
		                                         A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                         q_AlphaBeta, a_AlphaBeta);
		
		atoma_force=atoma_force+SW2BodyForce;



	}

	for (Neigh bond_neighbor=bond_neighbors.begin();  bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	{

		Atom <dim> *bond_neigh;
		bond_neigh=*bond_neighbor;

		Point <dim> SW2BodyForce(0.0,0.0,0.0);

		SW2BodyForce=StillingerWeberTwoBodyForce((*atoma), (*bond_neigh),
		                                         sigma_AlphaBeta, epsilon_AlphaBeta,
		                                         A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                         q_AlphaBeta, a_AlphaBeta);
		
		atoma_force=atoma_force+SW2BodyForce;



	}

	return(atoma_force);
}

template <int dim>
Point <dim> Force <dim>::ResultantSWThreeBodyForce( Atom <dim> *atoma, double sigma_AlphaBeta,
                                                    double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
		                                            double mu, double lambda_AlphaBetaGamma
                                                   ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{

	typedef typename vector < Atom <dim>* >::const_iterator Neigh;

	this->atomi=*atoma;

	vector < Atom <dim>* > neighbors=atomi.Neighbor();
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();
	vector < Atom <dim>* > all_neighbors;
	
	Point <dim> atomaThreeBodyForce;

	atomaThreeBodyForce.SetXCoord(0.);
	atomaThreeBodyForce.SetYCoord(0.);
	atomaThreeBodyForce.SetZCoord(0.);

	for (Neigh beta=neighbors.begin(); beta!=neighbors.end(); ++beta)
	{
		Atom <dim> *Beta;
		Beta=*beta;
		
		int Beta_id=Beta->GetID();
		
		for (Neigh gamma=neighbors.begin(); gamma!=neighbors.end(); ++gamma)
		{
			
			Atom <dim> *Gamma;
			Gamma=*gamma;
			int Gamma_id=Gamma -> GetID();
			
			if (Beta_id < Gamma_id)
				
				{
				
					Point <dim> ThreeBodyForceH1(0.,0.,0.);
		
					ThreeBodyForceH1=StillingerWeberThreeBodyForceH1( Gamma,atoma,Beta,
																	  sigma_AlphaBeta,
																	  sigma_AlphaGamma,epsilon_AlphaBetaGamma,
																	  mu, lambda_AlphaBetaGamma
																	  ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
					
					atomaThreeBodyForce=atomaThreeBodyForce+ThreeBodyForceH1;
				}
			
			
		}
		
		for (Neigh bond_gamma=bond_neighbors.begin(); bond_gamma!=bond_neighbors.end(); ++bond_gamma)
		{
			
			
			Atom <dim> *Bond_Gamma;
			Bond_Gamma=*bond_gamma;
			
			int Bond_Gamma_id=Bond_Gamma -> GetID();
				
			Point <dim> ThreeBodyForceH1(0.,0.,0.);

			ThreeBodyForceH1=StillingerWeberThreeBodyForceH1( Bond_Gamma,atoma,Beta,
															  sigma_AlphaBeta,
															  sigma_AlphaGamma,epsilon_AlphaBetaGamma,
															  mu, lambda_AlphaBetaGamma
															  ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
			
			atomaThreeBodyForce=atomaThreeBodyForce+ThreeBodyForceH1;
			
			
		}
		

	}

	for (Neigh bond_beta=bond_neighbors.begin(); bond_beta!=bond_neighbors.end(); ++bond_beta)
	{
		
		Atom <dim> *Bond_Beta;
		Bond_Beta=*bond_beta;
		
		int Bond_Beta_id=Bond_Beta -> GetID();
		
		
		for (Neigh bond_gamma=bond_neighbors.begin(); bond_gamma!=bond_neighbors.end(); ++bond_gamma)
		{
			
			Atom <dim> *Bond_Gamma;
			Bond_Gamma=*bond_gamma;
			
			int Bond_Gamma_id=Bond_Gamma -> GetID();
			
			if (Bond_Beta_id < Bond_Gamma_id)
				{
			
					Point <dim> ThreeBodyForceH1(0.,0.,0.);
		
					ThreeBodyForceH1=StillingerWeberThreeBodyForceH1( Bond_Gamma,atoma,Bond_Beta,
																	  sigma_AlphaBeta,
																	  sigma_AlphaGamma,epsilon_AlphaBetaGamma,
																	  mu, lambda_AlphaBetaGamma
																	  ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
					
					atomaThreeBodyForce=atomaThreeBodyForce+ThreeBodyForceH1;
					
				}
			
			
		}


	}
	
	Point <dim> initial_force(0.,0.,0.);
	
	initial_force=(*atoma).GetForce();
	
	(*atoma).SetForce(initial_force+atomaThreeBodyForce);

	return(initial_force+atomaThreeBodyForce);

}

template <int dim>
Point <dim> Force <dim>:: deformationalBondSWForce(Atom <dim> &atom_Alpha, Atom <dim> &atom_Beta,
                                                   double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                   double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
                                                   double J_AlphaBeta, double a_AlphaBeta)
{

	this -> atomi=atom_Alpha;
	this -> atomj=atom_Beta;


	double material_bond_distance=bond.MaterialBondDistance(atomi,atomj);

	double sigma0_AlphaBeta=sigma_AlphaBeta/material_bond_distance;
	double epsilon0_AlphaBeta=epsilon_AlphaBeta/material_bond_distance;

	Point <dim> spatial_bond_normal=bond.SpatialBondNormal(atomi,atomj);

	double spatial_stretch=bond.SpatialBondStretch(atomi,atomj);

	Energy <dim> energy;

	double energy_AlphaBeta=energy.deformationalSWtwoBody(atomi, atomj,
                                                          sigma_AlphaBeta, epsilon_AlphaBeta,
                                                          A_AlphaBeta, B_AlphaBeta, I_AlphaBeta,
                                                          J_AlphaBeta, a_AlphaBeta);
	

	double term1=(I_AlphaBeta*B_AlphaBeta/spatial_stretch)*pow((sigma0_AlphaBeta/spatial_stretch),I_AlphaBeta)
			     -(J_AlphaBeta/spatial_stretch)*pow((sigma0_AlphaBeta/spatial_stretch),J_AlphaBeta);

	double term2=B_AlphaBeta*pow((sigma0_AlphaBeta/spatial_stretch),I_AlphaBeta)
	             -pow((sigma0_AlphaBeta/spatial_stretch),J_AlphaBeta);

	double term3=term1/term2;

	double term4=sigma0_AlphaBeta/pow((spatial_stretch-a_AlphaBeta*sigma0_AlphaBeta),2);

	Point <dim> BondSWForce(0.,0.,0.);

	BondSWForce.SetXCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetXCoord());
	BondSWForce.SetYCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetYCoord());
	BondSWForce.SetZCoord(energy_AlphaBeta*(term3+term4)*spatial_bond_normal.GetZCoord());
	
	return(BondSWForce);

}


template <int dim>
Point <dim> Force <dim>::deformationalprSW3BodyForce(Atom <dim> &atom_alpha, Atom <dim> &atom_beta, Atom <dim> &atom_gamma,
                                                     double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
	                                                 double mu, double eta_AlphaBetaGamma,
                                                     double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{

	this -> atomi = atom_alpha;
	this -> atomj = atom_beta;
	this -> atomk = atom_gamma;
	
	Point <dim> parallelForce (0.,0.,0.);
	
	Energy <dim> energy;
	
	double deformational_energy=energy.deformationalSWthreeBody(atomi, atomj, atomk,
                                                                sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                                mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
	
	double mat_dist_alphabeta=bond.MaterialBondDistance(atomj, atomi);
	double mat_dist_alphagamma=bond.MaterialBondDistance(atomk, atomi);
	double mat_dist_betagamma=bond.MaterialBondDistance(atomj, atomk);

	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;

	double area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));
	
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;

	double spatial_stretch_AlphaBeta=bond.SpatialBondStretch(atomj,atomi);
	double spatial_stretch_AlphaGamma=bond.SpatialBondStretch(atomk,atomi);
	
	if(spatial_stretch_AlphaBeta<(a_AlphaBeta*sigma_0_AlphaBeta)
	   && spatial_stretch_AlphaGamma<(a_AlphaGamma*sigma_0_AlphaGamma))
	{
	
	double term1=(mu*sigma_0_AlphaBeta/pow((spatial_stretch_AlphaBeta-(a_AlphaBeta*sigma_0_AlphaBeta)),2))*(1/mat_dist_alphabeta);
	double term2=(mu*sigma_0_AlphaGamma/pow((spatial_stretch_AlphaGamma-(a_AlphaGamma*sigma_0_AlphaGamma)),2))*(1/mat_dist_alphagamma);

	Point <dim> spatial_norm_AlphaBeta=bond.SpatialBondNormal(atomj,atomi);
	Point <dim> spatial_norm_AlphaGamma=bond.SpatialBondNormal(atomk,atomi);
	
	parallelForce.SetXCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetXCoord()
			                 +term2*spatial_norm_AlphaGamma.GetXCoord()));
	
	parallelForce.SetYCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetYCoord()
                            +term2*spatial_norm_AlphaGamma.GetYCoord()));
	
	parallelForce.SetZCoord(area*deformational_energy*(term1*spatial_norm_AlphaBeta.GetZCoord()
                            +term2*spatial_norm_AlphaGamma.GetZCoord()));
	
	}
	
	return parallelForce;
	
}

template <int dim>
Point <dim> Force <dim>::deformationalcpSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
                                                     double sigma_AlphaBeta, double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
	                                                 double mu, double eta_AlphaBetaGamma,
                                                     double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{

	this -> atomi = *atom_alpha;
	this -> atomj = *atom_beta;
	this -> atomk = *atom_gamma;
	Point <dim> Result(0.,0.,0.);

	double mat_dist_alphabeta=bond.MaterialBondDistance(atomi, atomj);
	double mat_dist_alphagamma=bond.MaterialBondDistance(atomi, atomk);
	double mat_dist_betagamma=bond.MaterialBondDistance(atomj, atomk);

	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;

	double area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));

	double epslion_0_AlphaBetaGamma=epsilon_AlphaBetaGamma/area;

	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;

	double spatial_stretch_AlphaBeta=bond.SpatialBondStretch(atomi,atomj);
	double spatial_stretch_AlphaGamma=bond.SpatialBondStretch(atomi,atomk);
	
	Point <dim> normalVectorAlphaBeta=bond.SpatialBondNormal(atomj, atomi);
	double spa_AlphaBeta_x=normalVectorAlphaBeta.GetXCoord();
	double spa_AlphaBeta_y=normalVectorAlphaBeta.GetYCoord();
	double spa_AlphaBeta_z=normalVectorAlphaBeta.GetZCoord();
	
	Point <dim> normalVectorAlphaGamma=bond.SpatialBondNormal(atomk, atomi);
	double spa_AlphaGamma_x=normalVectorAlphaGamma.GetXCoord();
	double spa_AlphaGamma_y=normalVectorAlphaGamma.GetYCoord();
	double spa_AlphaGamma_z=normalVectorAlphaGamma.GetZCoord();
	
	Point <dim> normalVectorBetaGamma=bond.SpatialBondNormal(atomk, atomj);
	double spa_BetaGamma_x=normalVectorBetaGamma.GetXCoord();
	double spa_BetaGamma_y=normalVectorBetaGamma.GetYCoord();
	double spa_BetaGamma_z=normalVectorBetaGamma.GetZCoord();
	
	Point <dim> normalVectorAlphaBetaMat=bond.MaterialBondNormal(atomj, atomi);
	double trace_matAlphaBeta=0;
	trace_matAlphaBeta=normalVectorAlphaBetaMat.GetXCoord()+normalVectorAlphaBetaMat.GetYCoord()+normalVectorAlphaBetaMat.GetZCoord();
	
	Point <dim> normalVectorAlphaGammaMat=bond.MaterialBondNormal(atomk, atomi);
	double trace_matAlphaGamma=0;
	trace_matAlphaGamma=normalVectorAlphaGammaMat.GetXCoord()+normalVectorAlphaGammaMat.GetYCoord()+normalVectorAlphaGammaMat.GetZCoord();

	double cosine_teta_BetaAlphaGamma=normalVectorAlphaBeta.GetXCoord()*normalVectorAlphaGamma.GetXCoord()
			                          +normalVectorAlphaBeta.GetYCoord()*normalVectorAlphaGamma.GetYCoord()
									  +normalVectorAlphaBeta.GetZCoord()*normalVectorAlphaGamma.GetZCoord();
	
	double exponet_AlphaBeta=0.;
	exponet_AlphaBeta=exp(mu*sigma_0_AlphaBeta/(spatial_stretch_AlphaBeta-(a_AlphaBeta*sigma_0_AlphaBeta)));
	
	double exponet_AlphaGamma=0.;
	exponet_AlphaGamma=exp(mu*sigma_0_AlphaGamma/(spatial_stretch_AlphaGamma-(a_AlphaGamma*sigma_0_AlphaGamma)));
	
	Point <dim> spatial_bond_vec_AlphaGamma=bond.SpatialBondVec(atomk, atomi);
	double spatial_bond_vec_AlphaGamma_x=spatial_bond_vec_AlphaGamma.GetXCoord();
	double spatial_bond_vec_AlphaGamma_y=spatial_bond_vec_AlphaGamma.GetYCoord();
	double spatial_bond_vec_AlphaGamma_z=spatial_bond_vec_AlphaGamma.GetZCoord();
	
	Point <dim> spatial_bond_vec_AlphaBeta=bond.SpatialBondVec(atomj, atomi);
	double spatial_bond_vec_AlphaBeta_x=spatial_bond_vec_AlphaBeta.GetXCoord();
	double spatial_bond_vec_AlphaBeta_y=spatial_bond_vec_AlphaBeta.GetYCoord();
	double spatial_bond_vec_AlphaBeta_z=spatial_bond_vec_AlphaBeta.GetZCoord();
	
	Point <dim> material_bond_vec_AlphaBeta=bond.MaterialBondVec(atomj, atomi);
	double material_bond_vec_AlphaBeta_x=material_bond_vec_AlphaBeta.GetXCoord();
	double material_bond_vec_AlphaBeta_y=material_bond_vec_AlphaBeta.GetYCoord();
	double material_bond_vec_AlphaBeta_z=material_bond_vec_AlphaBeta.GetZCoord();
	
	Point <dim> material_bond_vec_AlphaGamma=bond.MaterialBondVec(atomk, atomi);
	double material_bond_vec_AlphaGamma_x=material_bond_vec_AlphaGamma.GetXCoord();
	double material_bond_vec_AlphaGamma_y=material_bond_vec_AlphaGamma.GetYCoord();
	double material_bond_vec_AlphaGamma_z=material_bond_vec_AlphaGamma.GetZCoord();
	
	double spatial_bond_distance_alphabeta=bond.SpatialBondDistance(atomj,atomi);
	
	double spatial_bond_distance_alphagamma=bond.SpatialBondDistance(atomk,atomi);
	
	double spatial_bond_distance_betagamma=bond.SpatialBondDistance(atomk,atomj);
	
	double coefficient=2*eta_AlphaBetaGamma*epslion_0_AlphaBetaGamma*area*
			            exponet_AlphaBeta*exponet_AlphaGamma
						*(cosine_teta_BetaAlphaGamma-cosine_teta0);

	if (spatial_stretch_AlphaBeta < a_AlphaBeta*sigma_0_AlphaBeta
		&& spatial_stretch_AlphaGamma < a_AlphaGamma*sigma_0_AlphaGamma)
		
	{
		
		
		Result.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_x/spatial_bond_distance_alphabeta+
				                      spa_AlphaGamma_x/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_x/spatial_bond_distance_alphabeta)
				                      +(spa_AlphaBeta_x/spatial_bond_distance_alphagamma)));
				
		
		Result.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_y/spatial_bond_distance_alphabeta+
                                      spa_AlphaGamma_y/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_y/spatial_bond_distance_alphabeta)
                                      +(spa_AlphaBeta_y/spatial_bond_distance_alphagamma)));
		
		Result.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_z/spatial_bond_distance_alphabeta+
                                      spa_AlphaGamma_z/spatial_bond_distance_alphagamma)+(spa_AlphaGamma_z/spatial_bond_distance_alphabeta)
                                      +(spa_AlphaBeta_z/spatial_bond_distance_alphagamma)));
		
		Point <dim> deform_force_atom_beta(0.,0.,0.);
		
		deform_force_atom_beta.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_x/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_x/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_x);
		
		deform_force_atom_beta.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_y/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_y/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_y);
		
		deform_force_atom_beta.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaBeta_z/spatial_bond_distance_alphabeta)
                                         +(spa_AlphaBeta_z/spatial_bond_distance_alphagamma))
				                         +(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
										 spa_BetaGamma_z);
		
		atom_beta->SetDeformForce(atom_beta->GetDeformForce()+deform_force_atom_beta);
		
		Point <dim> deform_force_atom_gamma(0.,0.,0.);
		
		deform_force_atom_gamma.SetXCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_x/spatial_bond_distance_alphagamma)
				                          +(spa_AlphaGamma_x/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);
		
		deform_force_atom_gamma.SetYCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_y/spatial_bond_distance_alphagamma)
                                          +(spa_AlphaGamma_y/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);
		
		deform_force_atom_gamma.SetZCoord(coefficient*(-1*cosine_teta_BetaAlphaGamma*(spa_AlphaGamma_z/spatial_bond_distance_alphagamma)
                                          +(spa_AlphaGamma_z/spatial_bond_distance_alphabeta))
				                          -(spatial_bond_distance_betagamma/(spatial_bond_distance_alphabeta*spatial_bond_distance_alphagamma))*
														 spa_BetaGamma_x);
		
		atom_gamma->SetDeformForce(atom_gamma->GetDeformForce()+deform_force_atom_gamma);
		
	}
	
	return Result;

}


template <int dim>
Point <dim> Force <dim>::configurationalBondSWForce(Atom <dim> &atom_alpha,Atom <dim> &atom_beta,
		                                            double sigma_AlphaBeta, double epsilon_AlphaBeta,
		                                            double A_AlphaBeta, double B_AlphaBeta, double I_AlphaBeta,
					                                double J_AlphaBeta, double a_AlphaBeta)
{
	this -> atomi = atom_alpha;
	this -> atomj = atom_beta;
	
	double mat_dist_alphabeta=bond.MaterialBondDistance(atomi, atomj);

	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double epsilon_0_AlphaBeta=epsilon_AlphaBeta/mat_dist_alphabeta;
	
	double material_AlphaBeta=bond.MaterialBondStretch(atomj,atomi);
	
	Point <dim> material_bond_normal=bond.MaterialBondNormal(atomj,atomi);
	
	Energy <dim> energy;
	
	double Config2BodyEng=energy.configurationalSWtwoBody(atomi, atomj,
                                                          sigma_AlphaBeta, epsilon_AlphaBeta,
                                                          A_AlphaBeta,  B_AlphaBeta, I_AlphaBeta,
                                                          J_AlphaBeta, a_AlphaBeta);
	
	double term0=1/material_AlphaBeta;
	
	double term1=(B_AlphaBeta*I_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),I_AlphaBeta)
	             -J_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),J_AlphaBeta))
	            		 
	              /(B_AlphaBeta*material_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),I_AlphaBeta)
                    -material_AlphaBeta*pow((sigma_0_AlphaBeta*material_AlphaBeta),J_AlphaBeta));
	
	double term2=sigma_0_AlphaBeta/pow((1-(material_AlphaBeta*a_AlphaBeta*sigma_0_AlphaBeta)),2);
	
	double term3=-Config2BodyEng*(term0+term1+term2);
	
	Point <dim> config_BondSWForce(0.,0.,0.);
	
	config_BondSWForce.SetXCoord(term3*material_bond_normal.GetXCoord());
	config_BondSWForce.SetYCoord(term3*material_bond_normal.GetYCoord());
	config_BondSWForce.SetZCoord(term3*material_bond_normal.GetZCoord());
	
	return config_BondSWForce;
	
}


template <int dim>
Point <dim> Force <dim>::configurationalcpSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
		                                               double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
									                   double mu, double eta_AlphaBetaGamma,
                                                       double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{
	
	Point <dim> config_cpSW3BodyForce(0.,0.,0.);
	
	this -> atomi = *atom_alpha;
	this -> atomj = *atom_beta;
	this -> atomk = *atom_gamma;
		
	double spa_dist_alphabeta=bond.SpatialBondDistance(atomj, atomi);
	double spa_dist_alphagamma=bond.SpatialBondDistance(atomk, atomi);
	double spa_dist_betagamma=bond.SpatialBondDistance(atomj, atomk);
	
	double mat_dist_alphabeta=bond.MaterialBondDistance(atomj, atomi);
	double mat_dist_alphagamma=bond.MaterialBondDistance(atomk, atomi);
	double mat_dist_betagamma=bond.MaterialBondDistance(atomj, atomk);

	double m=(spa_dist_alphabeta+spa_dist_alphagamma+spa_dist_betagamma)/2;
	double M=(mat_dist_alphabeta+mat_dist_alphagamma+mat_dist_betagamma)/2;

	double area=sqrt(m*abs(m-spa_dist_alphabeta)*abs(m-spa_dist_alphagamma)*abs(m-spa_dist_betagamma));
	double Area=sqrt(M*abs(M-mat_dist_alphabeta)*abs(M-mat_dist_alphagamma)*abs(M-mat_dist_betagamma));
	
	double ratio_areas=Area/area;
	
	Energy <dim> energy;
	
	double config_energy=energy.configurationalSWthreeBody(atomi,atomj,atomk, sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                           mu, eta_AlphaBetaGamma,
                                                           cosine_teta0, a_AlphaBeta, a_AlphaGamma);
	
	double epsilon_AlphaBetaGamma_t=epsilon_AlphaBetaGamma/area;
	double epsilon_AlphaBetaGamma_0=epsilon_AlphaBetaGamma/Area;
	
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;

	double material_stretch_AlphaBeta=bond.MaterialBondStretch(atomj,atomi);
	double material_stretch_AlphaGamma=bond.MaterialBondStretch(atomk,atomi);
	
	SecondTensor <double> rotation_mat_alphabeta(3,3,0.);
	rotation_mat_alphabeta=bond.RotationMatrix(atomj,atomi);
	SecondTensor <double> rotation_mat_alphabeta_transpose(3,3,0.);
	rotation_mat_alphabeta_transpose=rotation_mat_alphabeta.transpose();
	
	SecondTensor <double> rotation_mat_alphagamma(3,3,0.);
	rotation_mat_alphagamma=bond.RotationMatrix(atomj,atomi);
	SecondTensor <double> rotation_mat_alphagamma_transpose(3,3,0.);
	rotation_mat_alphagamma_transpose=rotation_mat_alphagamma.transpose();
	
	Point <dim> spa_alphabeta=bond.SpatialBondVec(atomj,atomi);
	std::vector <double> spa_alphabeta_vec;
	spa_alphabeta_vec.push_back(spa_alphabeta.GetXCoord());
	spa_alphabeta_vec.push_back(spa_alphabeta.GetYCoord());
	spa_alphabeta_vec.push_back(spa_alphabeta.GetZCoord());
	
	Point <dim> spa_alphagamma=bond.SpatialBondVec(atomk,atomi);
	std::vector <double> spa_alphagamma_vec;
	spa_alphagamma_vec.push_back(spa_alphagamma.GetXCoord());
	spa_alphagamma_vec.push_back(spa_alphagamma.GetYCoord());
	spa_alphagamma_vec.push_back(spa_alphagamma.GetZCoord());
	
	std::vector <double> mat_alphabeta_vec=rotation_mat_alphabeta_transpose*spa_alphabeta_vec;
	std::vector <double> mat_alphagamma_vec=rotation_mat_alphagamma_transpose*spa_alphagamma_vec;
	
	double term1= mat_alphabeta_vec[0]*mat_alphagamma_vec[0]
                 +mat_alphabeta_vec[1]*mat_alphagamma_vec[1]
				 +mat_alphabeta_vec[2]*mat_alphagamma_vec[2];
	
	double term2=sqrt( pow(mat_alphabeta_vec[0],2)
			          +pow(mat_alphabeta_vec[1],2)
					  +pow(mat_alphabeta_vec[2],2))*
				 sqrt( pow(mat_alphagamma_vec[0],2)
				       +pow(mat_alphagamma_vec[1],2)
				       +pow(mat_alphagamma_vec[2],2));
							  
	
	double cosine_theta=term1/term2;
	double diff_cosines=cosine_theta-cosine_teta0;
			
	double term3=material_stretch_AlphaBeta*mu*sigma_0_AlphaBeta
				 /(1-(material_stretch_AlphaBeta*a_AlphaBeta*sigma_0_AlphaBeta));

	double term4=material_stretch_AlphaGamma*mu*sigma_0_AlphaGamma
				 /(1-(material_stretch_AlphaGamma*a_AlphaGamma*sigma_0_AlphaGamma));
	
	Point <dim> mat_bond_normal_alphabeta=bond.MaterialBondNormal(atomj,atomi);
	Point <dim> mat_bond_normal_alphagamma=bond.MaterialBondNormal(atomk,atomi);
	Point <dim> mat_bond_normal_betagamma=bond.MaterialBondNormal(atomk,atomj);
	
	if(material_stretch_AlphaBeta < a_AlphaBeta*sigma_0_AlphaBeta
	   && material_stretch_AlphaGamma < a_AlphaGamma*sigma_0_AlphaGamma)
		
	{
		
		config_cpSW3BodyForce.SetXCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphagamma)))
														
										-(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetXCoord()
										+((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetXCoord()));
		
		config_cpSW3BodyForce.SetYCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphagamma)))
										 -(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetYCoord()
										 +((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetYCoord()));
		
		config_cpSW3BodyForce.SetZCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
				                        *diff_cosines*exp(term3)*exp(term4)*
										((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphabeta)+
												        (mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphagamma)))
										-(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphabeta.GetZCoord()
										+((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))/(2*Area))*mat_bond_normal_alphagamma.GetZCoord()));
		
		Point <dim> config_force_atom_beta(0.,0.,0.);
		
		config_force_atom_beta.SetXCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetXCoord())
										  
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetXCoord())
										 );	
		
		config_force_atom_beta.SetYCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetYCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetXCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetYCoord())
										  
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetYCoord())
										 );
		
		config_force_atom_beta.SetZCoord(-2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
										 *diff_cosines*exp(term3)*exp(term4)*
										 ((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
										 +(mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphagamma)
										 -cosine_theta*((mat_bond_normal_alphabeta.GetZCoord()/mat_dist_alphabeta))
										 +(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										 *mat_bond_normal_betagamma.GetXCoord())
										  
										 +(config_energy/Area)*(((M*(M-mat_dist_alphagamma)*(M-mat_dist_betagamma))
										 /(2*Area))*mat_bond_normal_alphabeta.GetZCoord())
										 );
		
		atom_beta->SetForce(atom_beta->GetForce()+config_force_atom_beta);
		
		
		Point <dim> config_force_atom_gamma(0.,0.,0.);
		
		config_force_atom_gamma.SetXCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
				                          ((mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphabeta)
				                          -cosine_theta*(mat_bond_normal_alphagamma.GetXCoord()/mat_dist_alphagamma))
										  
				                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
				                          /(2*Area))*mat_bond_normal_alphagamma.GetXCoord()));	
		
		config_force_atom_gamma.SetYCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
                                          ((mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphabeta)
                                          -cosine_theta*(mat_bond_normal_alphagamma.GetYCoord()/mat_dist_alphagamma)
										  -(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										  *mat_bond_normal_betagamma.GetYCoord())
										  
                                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
                                          /(2*Area))*mat_bond_normal_alphagamma.GetYCoord()));
		
		config_force_atom_gamma.SetZCoord(2*epsilon_AlphaBetaGamma_0*eta_AlphaBetaGamma*ratio_areas
                                          *diff_cosines*exp(term3)*exp(term4)*
                                          ((mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphabeta)
                                          -cosine_theta*(mat_bond_normal_alphagamma.GetZCoord()/mat_dist_alphagamma))
										  -(mat_dist_betagamma/(mat_dist_alphagamma*mat_dist_alphabeta))
										  *mat_bond_normal_betagamma.GetZCoord()
										  
                                          -(config_energy/Area)*(((M*(M-mat_dist_alphabeta)*(M-mat_dist_betagamma))
                                          /(2*Area))*mat_bond_normal_alphagamma.GetZCoord()));
		
		atom_gamma->SetForce(atom_gamma->GetForce()+config_force_atom_gamma);
		
		
	}
	
	
	return config_cpSW3BodyForce;
}
																   
																


template <int dim>
Point <dim> Force <dim>::configurationalprSW3BodyForce(Atom <dim> *atom_alpha, Atom <dim> *atom_beta, Atom <dim> *atom_gamma,
		                                                           double sigma_AlphaBeta,double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
									                               double mu, double eta_AlphaBetaGamma,
                                                                   double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
																   {
	Point <dim> config_prSW3BodyForce(0.,0.,0.);
	
	this -> atomi = *atom_alpha;
	this -> atomj = *atom_beta;
	this -> atomk = *atom_gamma;
	
	Energy <dim> energy;
	
	double configurational_energy=energy.configurationalSWthreeBody(atomi, atomj, atomk,
                                                                    sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                                    mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
	
	double spa_dist_alphabeta=bond.SpatialBondDistance(atomj, atomi);
	double spa_dist_alphagamma=bond.SpatialBondDistance(atomk, atomi);
	double spa_dist_betagamma=bond.SpatialBondDistance(atomj, atomk);
	
	double mat_dist_alphabeta=bond.MaterialBondDistance(atomj, atomi);
	double mat_dist_alphagamma=bond.MaterialBondDistance(atomk, atomi);

	double m=(spa_dist_alphabeta+spa_dist_alphagamma+spa_dist_betagamma)/2;

	double area=sqrt(m*abs(m-spa_dist_alphabeta)*abs(m-spa_dist_alphagamma)*abs(m-spa_dist_betagamma));
	
	double sigma_0_AlphaBeta=sigma_AlphaBeta/mat_dist_alphabeta;
	double sigma_0_AlphaGamma=sigma_AlphaGamma/mat_dist_alphagamma;

	double material_stretch_AlphaBeta=bond.MaterialBondStretch(atomj,atomi);
	double material_stretch_AlphaGamma=bond.MaterialBondStretch(atomk,atomi);
	
	Point <dim> material_bond_normal_AlphaBeta=bond.MaterialBondNormal(atomj,atomi);
	Point <dim> material_bond_normal_AlphaGamma=bond.MaterialBondNormal(atomk,atomi);
	
	double term1=(mu*sigma_0_AlphaBeta)/((1-(material_stretch_AlphaBeta)*a_AlphaBeta*sigma_0_AlphaBeta));
	double term2=(mu*sigma_0_AlphaGamma)/((1-(material_stretch_AlphaGamma)*a_AlphaGamma*sigma_0_AlphaGamma));
	
	double config_energy=energy.configurationalSWthreeBody(atomi, atomj, atomk,
                                                           sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
                                                           mu, eta_AlphaBetaGamma,cosine_teta0, a_AlphaBeta, a_AlphaGamma);
	
	config_prSW3BodyForce.SetXCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetXCoord()
			                                            +term1*material_bond_normal_AlphaGamma.GetXCoord()));
	
	config_prSW3BodyForce.SetYCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetYCoord()
                                                         +term2*material_bond_normal_AlphaGamma.GetYCoord()));
	
	config_prSW3BodyForce.SetZCoord(-config_energy*(term1*material_bond_normal_AlphaBeta.GetZCoord()
                                                        +term2*material_bond_normal_AlphaGamma.GetZCoord()));
	
	Point <dim> config_forcepr_atom_beta(0.,0.,0.);
	
	config_forcepr_atom_beta.SetXCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetXCoord()));
	config_forcepr_atom_beta.SetYCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetYCoord()));
	config_forcepr_atom_beta.SetZCoord(config_energy*(term1*material_bond_normal_AlphaBeta.GetZCoord()));
	
	atom_beta->SetForce(atom_beta->GetForce()+config_forcepr_atom_beta);
	
	Point <dim> config_forcepr_atom_gamma(0.,0.,0.);
	
	config_forcepr_atom_gamma.SetXCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetXCoord()));
	config_forcepr_atom_gamma.SetYCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetYCoord()));
	config_forcepr_atom_gamma.SetZCoord(config_energy*(term1*material_bond_normal_AlphaGamma.GetZCoord()));
	
	atom_gamma->SetForce(atom_gamma->GetForce()+config_forcepr_atom_gamma);
	
	return config_prSW3BodyForce;
	
																   }




template <int dim>
Point <dim> Force <dim>::ResultantConfigSWTwoBodyForce(Atom <dim> *atoma, double sigma_AlphaBeta, double epsilon_AlphaBeta,
                                                       double A_AlphaBeta, double B_AlphaBeta, double p_AlphaBeta,
                                                       double q_AlphaBeta, double a_AlphaBeta)
{

	this->atomi=*atoma;

	Point <dim> atoma_force(0.,0.,0.);

	vector < Atom <dim>* > neighbors=atomi.Neighbor();
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();

	typedef typename vector < Atom <dim>* >::const_iterator Neigh;

	for (Neigh neighbor=neighbors.begin();  neighbor!=neighbors.end(); ++neighbor)
	{

		Atom <dim> *neigh;
		neigh=*neighbor;

		Point <dim> SW2BodyForce(0.0,0.0,0.0);

		SW2BodyForce=configurationalBondSWForce((*atoma), (*neigh),
		                                         sigma_AlphaBeta, epsilon_AlphaBeta,
		                                         A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                         q_AlphaBeta, a_AlphaBeta);

		atoma_force.SetXCoord(atoma_force.GetXCoord()+SW2BodyForce.GetXCoord());
		atoma_force.SetYCoord(atoma_force.GetYCoord()+SW2BodyForce.GetYCoord());
		atoma_force.SetZCoord(atoma_force.GetZCoord()+SW2BodyForce.GetZCoord());


	}

	for (Neigh bond_neighbor=bond_neighbors.begin();  bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
	{

		Atom <dim> *bond_neigh;
		bond_neigh=*bond_neighbor;

		Point <dim> SW2BodyForce(0.0,0.0,0.0);

		SW2BodyForce=configurationalBondSWForce((*atoma), (*bond_neigh),
		                                         sigma_AlphaBeta, epsilon_AlphaBeta,
		                                         A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
		                                         q_AlphaBeta, a_AlphaBeta);

		atoma_force.SetXCoord(atoma_force.GetXCoord()+SW2BodyForce.GetXCoord());
		atoma_force.SetYCoord(atoma_force.GetYCoord()+SW2BodyForce.GetYCoord());
		atoma_force.SetZCoord(atoma_force.GetZCoord()+SW2BodyForce.GetZCoord());


	}

	return(atoma_force);
}

template <int dim>
Point <dim> Force <dim>::ResultantConfigSWThreeBodyForce( Atom <dim> *atoma, double sigma_AlphaBeta,
                                                    double sigma_AlphaGamma,double epsilon_AlphaBetaGamma,
		                                            double Gamma, double lambda_AlphaBetaGamma
                                                   ,double cosine_teta0, double a_AlphaBeta, double a_AlphaGamma)
{

	typedef typename vector < Atom <dim>* >::const_iterator Neigh;

	this->atomi=*atoma;

	vector < Atom <dim>* > neighbors=atomi.Neighbor();
	vector < Atom <dim>* > bond_neighbors=atomi.BondNeighbor();

	Point <dim> atomaThreeBodyForce(0.,0.,0.);

	for (Neigh beta=neighbors.begin(); beta!=neighbors.end(); ++beta)
	{

		int beta_id=(*beta) -> GetID();

		Atom <dim> *neigh_beta;
		neigh_beta=*beta;

		for (Neigh gamma=neighbors.begin(); gamma!=neighbors.end(); ++gamma)
		{

			int gamma_id=0.;

			gamma_id=(*gamma) -> GetID();

			Atom <dim> *neigh_gamma;
			neigh_gamma=*gamma;

			if (gamma_id > beta_id)
			{

				Point <dim> ThreeBodyForceH1(0.,0.,0.);

				ThreeBodyForceH1=configurationalcpSW3BodyForce(atoma,neigh_beta,neigh_gamma,
						                                       sigma_AlphaBeta,
						                                       sigma_AlphaGamma,epsilon_AlphaBetaGamma,
								                               Gamma, lambda_AlphaBetaGamma
						                                       ,cosine_teta0, a_AlphaBeta, a_AlphaGamma)
						                                        		  
						           +configurationalprSW3BodyForce(atoma,neigh_beta,neigh_gamma,
																  sigma_AlphaBeta,
																  sigma_AlphaGamma,epsilon_AlphaBetaGamma,
																  Gamma, lambda_AlphaBetaGamma
																  ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);

				atomaThreeBodyForce.SetXCoord(atomaThreeBodyForce.GetXCoord()+ThreeBodyForceH1.GetXCoord());

				atomaThreeBodyForce.SetYCoord(atomaThreeBodyForce.GetYCoord()+ThreeBodyForceH1.GetYCoord());

				atomaThreeBodyForce.SetZCoord(atomaThreeBodyForce.GetZCoord()+ThreeBodyForceH1.GetZCoord());


			}


		}
		
		
		for (Neigh bond_gamma=bond_neighbors.begin(); bond_gamma!=bond_neighbors.end(); ++bond_gamma)
		{

			int bond_gamma_id=0.;
			bond_gamma_id=(*bond_gamma) -> GetID();

			Atom <dim> *bond_neigh_gamma;
			bond_neigh_gamma=*bond_gamma;

			if (bond_gamma_id > beta_id)
			{
				
				Point <dim> ThreeBodyForceH1(0.,0.,0.);

				ThreeBodyForceH1=configurationalcpSW3BodyForce(atoma,neigh_beta,bond_neigh_gamma,
						                                        sigma_AlphaBeta,
						                                        sigma_AlphaGamma,epsilon_AlphaBetaGamma,
								                                Gamma, lambda_AlphaBetaGamma
						                                        ,cosine_teta0, a_AlphaBeta, a_AlphaGamma)
						                                        		  
						           +configurationalprSW3BodyForce(atoma,neigh_beta,bond_neigh_gamma,
	                                                               sigma_AlphaBeta,
	                                                               sigma_AlphaGamma,epsilon_AlphaBetaGamma,
			                                                       Gamma, lambda_AlphaBetaGamma
	                                                               ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);


				atomaThreeBodyForce.SetXCoord(atomaThreeBodyForce.GetXCoord()+ThreeBodyForceH1.GetXCoord());

				atomaThreeBodyForce.SetYCoord(atomaThreeBodyForce.GetYCoord()+ThreeBodyForceH1.GetYCoord());

				atomaThreeBodyForce.SetZCoord(atomaThreeBodyForce.GetZCoord()+ThreeBodyForceH1.GetZCoord());


			}


		}


	}

	for (Neigh bond_beta=bond_neighbors.begin(); bond_beta!=bond_neighbors.end(); ++bond_beta)
	{

		int beta_id=0.;

		beta_id=(*bond_beta) -> GetID();

		Atom <dim> *bond_neigh_beta;
		bond_neigh_beta=*bond_beta;

		for (Neigh bond_gamma=bond_neighbors.begin(); bond_gamma!=bond_neighbors.end(); ++bond_gamma)
		{

			int gamma_id=0.;
			gamma_id=(*bond_gamma) -> GetID();

			Atom <dim> *bond_neigh_gamma;
			bond_neigh_gamma=*bond_gamma;

			if (gamma_id > beta_id)
			{
				
				Point <dim> ThreeBodyForceH1(0.,0.,0.);

				ThreeBodyForceH1=configurationalcpSW3BodyForce(atoma,bond_neigh_beta,bond_neigh_gamma,
						                                          sigma_AlphaBeta,
						                                          sigma_AlphaGamma,epsilon_AlphaBetaGamma,
								                                  Gamma, lambda_AlphaBetaGamma
						                                          ,cosine_teta0, a_AlphaBeta, a_AlphaGamma)
						                                        		  
						           +configurationalprSW3BodyForce( atoma,bond_neigh_beta,bond_neigh_gamma,
	                                                               sigma_AlphaBeta,
	                                                               sigma_AlphaGamma,epsilon_AlphaBetaGamma,
			                                                       Gamma, lambda_AlphaBetaGamma
	                                                               ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);


				atomaThreeBodyForce.SetXCoord(atomaThreeBodyForce.GetXCoord()+ThreeBodyForceH1.GetXCoord());

				atomaThreeBodyForce.SetYCoord(atomaThreeBodyForce.GetYCoord()+ThreeBodyForceH1.GetYCoord());

				atomaThreeBodyForce.SetZCoord(atomaThreeBodyForce.GetZCoord()+ThreeBodyForceH1.GetZCoord());


			}


		}


	}
	
	Point <dim> initial_force(0.,0.,0.);
	
	initial_force=(*atoma).GetConfigForce();
	
	(*atoma).SetConfigForce(initial_force+atomaThreeBodyForce);

	return(initial_force+atomaThreeBodyForce);

}

template class Force<2>;
template class Force<3>;
