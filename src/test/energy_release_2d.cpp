/**
 * energy_release_2d.cpp
 *
 *  Created on: Aug 26, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//---
#include "assign_region.h"
#include "atom.h"
#include "config_force_criterion.h"
#include "point.h"
#include "energy.h"
#include "force.h"
#include "boundary.h"
#include "bond.h"
#include "matrix.h"
#include "neighbors.h"
#include "unrelaxed_config.h"
#include "write_data_file.h"

//#include <eigen3/Eigen/Core>
#include "ceres/ceres.h"
#include "glog/logging.h"
//#include "gflags/gflags.h"


class MolecularStatic : public ceres::FirstOrderFunction
{

private:

	const vector < Atom <2>* > Atoms;
	double sigma;
	double epsilon;
	int dof;

public:


    MolecularStatic(const vector < Atom <2>* > atoms, double sigma0, double epsilon0, int dof0):
    	           Atoms(atoms), sigma(sigma0), epsilon(epsilon0), dof(dof0){}

    typedef typename vector < Atom <2>* >::const_iterator At;


	virtual bool Evaluate (const double* parameters, double* cost, double* gradient) const
	{

	    Energy <2> energy(1.0, 1.0);
	    Force <2> force(1.0, 1.0);

	    int iter=0;

	    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
	    {
	    	//only for atoms that are labeled

	    	int atom_region=(*atom)->GetAtomRegion();
	    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

	    	{

        		Point <2> spatial_position_2d;

        		spatial_position_2d.SetXCoord(parameters[iter]);
        		spatial_position_2d.SetYCoord(parameters[iter+1]);

	        	(*atom)->SetSpatialPosition(spatial_position_2d);

	        	iter+=2;

	    	}

	    }


        int iter_force=0;
        if (gradient!=NULL)
        {

        	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
        	{

        		int atom_region=(*atom)->GetAtomRegion();
        		if( atom_region!=1 && atom_region!=2 && atom_region!=3 )
        		{

            		Point <2> resultant_force=force.ResultantForce(*atom);

            		gradient[iter_force]=-1*resultant_force.GetXCoord();
            		gradient[iter_force+1]=-1*resultant_force.GetYCoord();

            		iter_force+=2;

        		}


        	}

        }

    	cost[0]=energy.TotPotentialEnergy(Atoms);

        return true;

	}

	virtual int NumParameters() const
	{
		return{dof};

	}

};



int main(){

    vector < Atom <2>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <2> (53, 49, 0, 1.1, 2, "model");

    Boundary <2> zone;

    double zone_x_min=2.08;
    double zone_y_min=1.9;

    double zone_x_max=55.12;
    double zone_y_max=43.9;
    int zone_id=0;

    zone.SetBoundaryRegion(zone_x_min, zone_y_min, zone_x_max,
    		               zone_y_max, zone_id);

    Boundary <2> boundary_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;

    double bottom_BC_x_max=58.;
    double bottom_BC_y_max=0.1;
    int bottom_BC_id=1;

    boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min, bottom_BC_x_max,
    		                          bottom_BC_y_max, bottom_BC_id);

    Boundary <2> boundary_top;

    double top_BC_x_min=0.0;
    double top_BC_y_min=45.2;

    double top_BC_x_max=58.;
    double top_BC_y_max=45.8;
    int top_BC_id=2;

    boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min,top_BC_x_max,
    		                       top_BC_y_max, top_BC_id);

    Boundary <2> crack_top;

    double crack_top_x_min=0;
    double crack_top_y_min=12.0*1.1*sqrt(3)+0.1;

    double crack_top_x_max=13.5*1.1+0.1;
    double crack_top_y_max=13.0*1.1*sqrt(3)+0.1;
    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min,
        		               crack_top_x_max, crack_top_y_max,crack_top_id);

    Boundary <2> crack_bottom;

    double crack_bottom_x_min=0;
    double crack_bottom_y_min=11.0*1.1*sqrt(3)+0.1;

    double crack_bottom_x_max=13.5*1.1+0.1;
    double crack_bottom_y_max=12.0*1.1*sqrt(3)+0.1;
    int crack_bottom_id=4;

    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min,
                                   crack_bottom_x_max, crack_bottom_y_max, crack_bottom_id);

//    Boundary <2> crack;
//
//    double crack_x_min=0;
//    double crack_y_min=18;
//
//    double crack_x_max=16.;
//    double crack_y_max=18.2;
//    int crack_id=3;
//
//    crack.SetBoundaryRegion(crack_x_min, crack_y_min, crack_x_max, crack_y_max,crack_id);


    vector < Boundary <2> > boundaries;
    boundaries.push_back(zone);
    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);

    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);

//    boundaries.push_back(crack);

    typedef typename vector < Atom <2>* >::iterator At;
    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,1.5);

	vector < Atom <2>* > optimized_atoms;
	int i=0;
    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ( (*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5 )
    	{
            optimized_atoms.push_back(*atom);
    	}

    }

    int num_optimized_atoms=optimized_atoms.size()*2;

	int iterator=0;
	int atom_num=0;
	double parameters [num_optimized_atoms];
	while (iterator<num_optimized_atoms)
	{

		Point <2> optimized_atom=optimized_atoms[atom_num]->GetMaterialPosition();
		double x=optimized_atom.GetXCoord();
		double y=optimized_atom.GetYCoord();

		parameters[iterator]=x;
		parameters[iterator+1]=y;

		iterator+=2;
		atom_num+=1;

	}

	string material_configuration="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/material_positions";

	ofstream mfile;
	mfile.open(material_configuration.c_str());


	if (mfile.is_open()){

		mfile <<"ITEM: TIMESTEP" <<endl;
		mfile << "1" <<endl;
		mfile <<"ITEM: NUMBER OF ATOMS" << endl;
		mfile << atoms.size() << endl;
		mfile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		mfile << "0" << " " << "1.6" << endl;
		mfile << "0" << " " << "1.6" << endl;
		mfile << "0" << " " << "1.6" << endl;
		mfile <<"ITEM: ATOMS id type x y z" <<endl;

		for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
		{
			Point <2> matposition = (*atom)-> GetMaterialPosition();
			int atom_id=(*atom)->GetID();

			int atom_region=(*atom)->GetAtomRegion();

			double materialp_x=matposition.GetXCoord();
			double materialp_y=matposition.GetYCoord();


			switch (atom_region) {

			case 0:  mfile << atom_id<< " " << "0" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << "0.0"<<endl;
			break;

			case 1: mfile << atom_id<< " " << "1" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << "0.0"<<endl;
			break;

			case 2: mfile << atom_id<< " " << "2" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << "0.0"<<endl;
			break;

			case 4: mfile << atom_id<< " " << "4" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << "0.0"<<endl;

			break;

			case 5: mfile << atom_id<< " " << "5" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << "0.0"<<endl;

			break;

			}
	    }
	}




//	int id=0;
//
//    double load=0;
//    double unload=0;
//    double reload=0;
//
//    int load_step=0;
//    int max_load_step=1680;
//    double totEngOld=0.;
//
//    Point <2> config_force_crack_tip(0.,0.);
//
//    string path2="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/total_energy";
//    ofstream efile;
//    efile.open(path2.c_str());
//
//    string path3="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/external_work";
//    ofstream wfile;
//    wfile.open(path3.c_str());
//
//    string path5="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/crack_extension";
//    ofstream bfile;
//    bfile.open(path5.c_str());
//
//    string path6="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/reaction_force";
//    ofstream RFfile;
//    RFfile.open(path6.c_str());
//
//
//
//    Energy <2> energy(1.0, 1.0);
//    Force <2> force(1.0, 1.0);
//    Bond <2> bond;
//    int atom_id_old=10000000;
//
//    int min_load=0;
//    int max_load=260;
//
//    int min_unload=261;
//    int max_unload=521;
//
//    int min_reload=522;
//    int max_reload=782;


//    while (load_step < max_load_step)
//    {
//
//
//
//    	int crack_length=0.;
//
//    	ceres::GradientProblem problem (new MolecularStatic(atoms, 1.0, 1.0, num_optimized_atoms));
//    	ceres::GradientProblemSolver::Options options;
//    	options.minimizer_progress_to_stdout=true;
//    	options.max_num_iterations=2000;
//    	options.gradient_tolerance=1e-10;
//    	options.line_search_direction_type=ceres::LBFGS;
//    	ceres::GradientProblemSolver::Summary summary;
//    	ceres::Solve(options, problem, parameters, &summary);
//
//    	string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/dump.crack";
//    	string path1="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/spatial/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/dump.crackspatial" ;
//    	string path4="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/spatial/screening/energy_release_model/CR1.5_bond_elimination_max_eng_release/dump.bonds" ;
//
//    	ofstream file;
//    	ofstream ofile;
//
//    	load_step+=1;
//    	string out_string;
//    	stringstream ss;
//    	ss << load_step;
//    	out_string=ss.str();
//
//    	string file_name=path+out_string;
//    	string file_name1=path1+out_string;
//
//    	file.open(file_name.c_str());
//    	ofile.open(file_name1.c_str());
//
//    	if (load_step==1)
//    	{
//
//    		for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//    		{
//
//    			Point <2> spatialP= (*atom)->GetSpatialPosition();
//    			(*atom)->SetMaterialPosition(spatialP);
//
//    		}
//
//    	}
//
//
//    	double external_work=0;
//
//    	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//
//    	{
//
//    		if((*atom)->GetAtomRegion()==2)
//    		{
//
//    			Point <2> atomEmpiricalForce=force.ResultantForce(*atom);
//    			double force_y=-1.0*atomEmpiricalForce.GetYCoord();
//    			force_y=floor(force_y*100.0)/100.0;
//    			external_work+=(force_y)*(0.005/2.0);
//
//    		}
//
//    		else if((*atom)->GetAtomRegion()==1)
//    		{
//
//    			Point <2> atomEmpiricalForce=force.ResultantForce(*atom);
//    			double force_y=-1.0*atomEmpiricalForce.GetYCoord();
//    			force_y=floor(force_y*100.0)/100.0;
//    			external_work+=(force_y)*(-0.005/2.0);
//
//    		}
//
//    	}
//
//    	wfile << load_step << " " << floor(external_work*100.0)/100.0 << endl;
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    		typedef typename vector < Atom <2>* >::iterator At;
//
//    		if ((load_step >= min_load && load_step <= max_load) ||
//    			(load_step >= min_reload && load_step <= max_reload))
//    		{
//
//				for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//
//				{
//					int atom_id=(*atom)->GetID();
//
//					Atom <2> *atm;
//					atm=*atom;
//
//					vector <Atom <2>*> update_neighbor_list;
//					vector <Atom <2>*> update_bond_neighbor_list;
//
//					Point <2> atomMaterialPosition= atm->GetMaterialPosition();
//
//					double atomMaterialPositionX=atomMaterialPosition.GetXCoord();
//					double atomMaterialPositionY=atomMaterialPosition.GetYCoord();
//
//					double crackRegion_MinX=14.57;
//					double crackRegion_MinY=22.86;
//
//					double crackRegion_MaxX=56.59;
//					double crackRegion_MaxY=23.85;
//
//					if ( atomMaterialPositionX >= crackRegion_MinX  && atomMaterialPositionY >= crackRegion_MinY
//							&& atomMaterialPositionX <= crackRegion_MaxX  && atomMaterialPositionY <= crackRegion_MaxY  )
//					{
//						vector < Atom <2>* > neighbors=atm->Neighbor();
//						vector < Atom <2>* > bond_neighbors=atm->BondNeighbor();
//
//						vector <Atom <2>*> update_neighbor_list;
//						vector <Atom <2>*> update_bond_neighbor_list;
//
//						for (At neighbor=neighbors.begin(); neighbor!=neighbors.end(); ++neighbor)
//						{
//							int neighbor_id=(*neighbor)->GetID();
//
//							Point <2> neighborMaterialPosition= (*neighbor) -> GetMaterialPosition();
//
//							double neighborMaterialPositionX=neighborMaterialPosition.GetXCoord();
//							double neighborMaterialPositionY=neighborMaterialPosition.GetYCoord();
//
//							Atom <2> *neigh;
//							neigh=*neighbor;
//
//							vector <Atom <2>*> update_neighbor_bond_neighbor_list;
//
//							vector < Atom <2>* > neighbor_neighbors=neigh->BondNeighbor();
//							int neigh_size=neighbor_neighbors.size();
//
//
//							double material_distance=bond.MaterialBondDistance(*atm,*neigh);
//							double initial_distance=bond.InitialBondDistance(*atm,*neigh);
//							double interatomic_stretch_0=bond.InitialBondStretch(*atm,*neigh);
//
//							double sigma0=1./(material_distance*initial_distance);
//							double epsilon0=1./(material_distance*initial_distance);
//							double critical_stretch=1.84;
//							double critical_spatial_distance=1.5;
//							double spatial_bond_distance=bond.SpatialBondDistance(*atm, *neigh);
//
//							double term=sigma0*(1./critical_stretch)*interatomic_stretch_0;
//
////							double criticalConfigForce=4*epsilon0*(13*pow(term,12)-7*pow(term,6));
//							double distance_spatial=bond.SpatialBondStretch(*atm, *neigh);
//
//							double configForce=0.0;
//
//							configForce=abs(force.interatomicConfigForceValue(*atm,*neigh));
//
//	    					double criticalConfigForce=(3.76661921*sigma0);
//
//	//    					|| configForce > abs(criticalConfigForce) || distance_spatial < critical_stretch
////							spatial_bond_distance < critical_spatial_distance
//
//							if ( configForce < abs(criticalConfigForce)  )
//
//							{
//
//								update_neighbor_list.push_back(*neighbor);
//
//							}
//
//	//    					 && configForce <= abs(criticalConfigForce)  && distance_spatial >= critical_stretch
////							spatial_bond_distance >= critical_spatial_distance
//
//							else if ( (configForce >= abs(criticalConfigForce))
//									  && neighborMaterialPositionX >= crackRegion_MinX  && neighborMaterialPositionY >= crackRegion_MinY
//									  && neighborMaterialPositionX<= crackRegion_MaxX  && neighborMaterialPositionY <= crackRegion_MaxY
//									  && abs(neighborMaterialPositionY-atomMaterialPositionY) > 0.6 )
//							{
//
//								cout << "criticalConfigForce: " << criticalConfigForce << " configForce: "<< abs(configForce) << " neighbor id : "<< neigh->GetID()<< " atom id: "<< atm->GetID() << endl;
//
//								neighbor_neighbors.erase(find(neighbor_neighbors.begin(), neighbor_neighbors.end(),(*atom)));
//								Point <2> interatomic_config_force=force.MaterialBondForce(*atm,*neigh);
//	//    						Point <2> material_distance=bond.MaterialBondVec(*atm,*neigh);
//
//								double interatomic_config_force_x=interatomic_config_force.GetXCoord();
//	//    						double material_bond_distance_x=material_distance.GetXCoord();
//
//
//								double material_distance=bond.MaterialBondDistance(*atm,*neigh);
//								double initial_distance=bond.InitialBondDistance(*atm,*neigh);
//
//								double sigma0=1./(material_distance*initial_distance);
//								double epsilon0=1./(material_distance*initial_distance);
//
//								Point <2> initial_bond=bond.InitialBondVec(*atm,*neigh);
//								double initial_bond_x=initial_bond.GetXCoord();
//
//								Point <2> material_bond=bond.MaterialBondVec(*atm, *neigh);
//								double material_bond_x=material_bond.GetXCoord();
//
//								Point <2> spatial_bond=bond.SpatialBondVec(*atm,*neigh);
//								double spatial_bond_x=spatial_bond.GetXCoord();
//
//								double energy_release=0.;
//
//								// \lambda_\mathcal{R}
//								double lambda_release=2.0*(1./interatomic_stretch_0);
//								// \lambda_{\mathcal{R}_0}
//								double lambda_initial_release=1.10865*(1./interatomic_stretch_0);
//
//								double first_term=pow((1./lambda_initial_release),11)-pow((1./lambda_release),11);
//								double first_coeff=(13./11.)*pow(sigma0,12);
//
//								double second_term=pow((1./lambda_release),5)-pow((1./lambda_initial_release),5);
//								double second_coeff=(7./5.)*pow(sigma0,6);
//
//								energy_release=4*epsilon0*(first_term*first_coeff+second_term*second_coeff);
//
//								// \lambda_\mathcal{B}
//								double lambda_breaking=1.75*(1./interatomic_stretch_0);
//								// \lambda_{\mathcal{B}_0}
//								double lambda_initial_breaking=1.223*(1./interatomic_stretch_0);
//
//								double bond_breaking_energy=0.0;
//
//								double first_div_breaking=(1./lambda_breaking);
//								double second_div_breaking=(1./lambda_initial_breaking);
//
//								double first_coeff_breaking=pow(sigma0,12);
//								double first_div_breaking_12=pow(first_div_breaking,12);
//								double second_div_breaking_12=pow(second_div_breaking,12);
//								double first_term_breaking=first_div_breaking_12-second_div_breaking_12;
//
//								double second_coeff_breaking=pow(sigma0,6);
//								double first_div_breaking_6=pow(first_div_breaking,6);
//								double second_div_breaking_6=pow(second_div_breaking,6);
//
//								double second_term_breaking=second_div_breaking_6-first_div_breaking_6;
//
//								bond_breaking_energy=4*epsilon0*(first_coeff_breaking*first_term_breaking-second_coeff_breaking*second_term_breaking);
//
//	//    						double bond_breaking_energy=0.;
//	//
//	//    						bond_breaking_energy=4*epsilon0*(pow((sigma0/(1.75*(1./interatomic_stretch_0))),12)-pow((sigma0/(1.75*(1./interatomic_stretch_0))),6))-
//	//    								             4*epsilon0*(pow((sigma0/(1.223*(1./interatomic_stretch_0))),12)-pow((sigma0/(1.223*(1./interatomic_stretch_0))),6));
//	//
//	//    						double energy_release=0.;
//	//
//	//    						energy_release=(4*epsilon0*(1./1.9))*
//	//    								       (pow(sigma0*(1./1.9),12)-pow(sigma0*(1./1.9),6))-
//	//    								       (4*epsilon0* (1./1.10865))*
//	//										   (pow(sigma0*(1./1.10865),12)-pow(sigma0*(1./1.10865),6));
//
//								bfile << load_step << " " << neighbor_id << " " << atom_id << " " << energy_release << " " << bond_breaking_energy
//									 << " " << abs(configForce) << " "<< initial_bond_x << " "<< material_bond_x << " " <<
//									 spatial_bond_x << " " << material_bond_x*initial_bond_x <<endl;
//
//							}
//
//	//    					 && configForce <= abs(criticalConfigForce) && distance_spatial >= critical_stretch
////							spatial_bond_distance >= critical_spatial_distance
//							else if ( configForce >= abs(criticalConfigForce)
//									  && (neighborMaterialPositionX <= crackRegion_MinX  || neighborMaterialPositionY <= crackRegion_MinY
//									  || neighborMaterialPositionX>= crackRegion_MaxX  || neighborMaterialPositionY >= crackRegion_MaxY
//									  || abs(neighborMaterialPositionY-atomMaterialPositionY) < 0.6))
//							{
//
//								update_neighbor_list.push_back(*neighbor);
//
//							}
//
//							(*neighbor)->SetBondNeighbors(neighbor_neighbors);
//
//							if (neighbor_neighbors.size()!= neigh_size)
//							{
//
//								cout << "neighbor_neighbors size: " << neighbor_neighbors.size() << " neighbor size: " <<  neigh_size << endl;
//
//							}
//
//						}
//
//
//						for (At bond_neighbor=bond_neighbors.begin(); bond_neighbor!=bond_neighbors.end(); ++bond_neighbor)
//						{
//
//							int bond_neighbor_id=(*bond_neighbor)->GetID();
//
//							Point <2> bond_neighborMaterialPosition= (*bond_neighbor) -> GetMaterialPosition();
//
//							double bond_neighborMaterialPositionX=bond_neighborMaterialPosition.GetXCoord();
//							double bond_neighborMaterialPositionY=bond_neighborMaterialPosition.GetYCoord();
//
//
//							Atom <2> *bond_neigh;
//							bond_neigh=*bond_neighbor;
//
//
//							vector <Atom <2>*> update_bond_neighbor_neighbor_list;
//
//							vector < Atom <2>* > bond_neighbor_neighbors=bond_neigh->Neighbor();
//
//							int bond_neigh_size=bond_neighbor_neighbors.size();
//
//							double material_distance=bond.MaterialBondDistance(*atm,*bond_neigh);
//							double initial_distance=bond.InitialBondDistance(*atm,*bond_neigh);
//							double interatomic_stretch_0=bond.InitialBondStretch(*atm,*bond_neigh);
//
//							double sigma0=1./(material_distance*initial_distance);
//							double epsilon0=1./(material_distance*initial_distance);
//							double critical_stretch=1.84;
//
//							double critical_spatial_distance=1.5;
//							double spatial_bond_distance=bond.SpatialBondDistance(*atm, *bond_neigh);
//
//							double term=sigma0*(1./critical_stretch)*interatomic_stretch_0;
//
////							double criticalConfigForce=4*epsilon0*(13*pow(term,12)-7*pow(term,6));
//	    					double criticalConfigForce=(3.76661921*sigma0);
//
//
//							double distance_spatial=bond.SpatialBondStretch(*atm, *bond_neigh);
//
//
//							double configForce=0.0;
//
//							configForce=abs(force.interatomicConfigForceValue(*atm,*bond_neigh));
//	//    					|| configForce > abs(criticalConfigForce) || distance_spatial < critical_stretch
////							spatial_bond_distance < critical_spatial_distance
//
//							if (configForce < abs(criticalConfigForce))
//
//							{
//
//								update_bond_neighbor_list.push_back(*bond_neighbor);
//
//							}
//
//	//    					&& configForce <= abs(criticalConfigForce) && distance_spatial >= critical_stretch
////							spatial_bond_distance >= critical_spatial_distance
//
//							else if ( configForce >= abs(criticalConfigForce)
//									  && bond_neighborMaterialPositionX >= crackRegion_MinX  && bond_neighborMaterialPositionY >= crackRegion_MinY
//									  && bond_neighborMaterialPositionX <= crackRegion_MaxX  && bond_neighborMaterialPositionY <= crackRegion_MaxY
//									  &&  abs(bond_neighborMaterialPositionY-atomMaterialPositionY) > 0.6 )
//							{
//
//								cout << "criticalConfigForce: " << criticalConfigForce << " configForce: "<< abs(configForce) << " bond_neighbor id : "<< bond_neigh->GetID()<< " atom id: "<< atm->GetID() << endl;
//
//								bond_neighbor_neighbors.erase(find(bond_neighbor_neighbors.begin(), bond_neighbor_neighbors.end(),(*atom)));
//								Point <2> interatomic_config_force=force.MaterialBondForce(*atm,*bond_neigh);
//	//    						Point <2> material_distance=bond.MaterialBondVec(*atm,*bond_neigh);
//
//								double interatomic_config_force_x=interatomic_config_force.GetXCoord();
//	//    						double material_bond_distance_x=material_distance.GetXCoord();
//								double material_distance=bond.MaterialBondDistance(*atm,*bond_neigh);
//								double initial_distance=bond.InitialBondDistance(*atm,*bond_neigh);
//
//								double sigma0=1./(material_distance*initial_distance);
//								double epsilon0=1./(material_distance*initial_distance);
//
//								Point <2> initial_bond=bond.InitialBondVec(*atm,*bond_neigh);
//								double initial_bond_x=initial_bond.GetXCoord();
//
//								Point <2> material_bond=bond.MaterialBondVec(*atm, *bond_neigh);
//								double material_bond_x=material_bond.GetXCoord();
//
//								Point <2> spatial_bond=bond.SpatialBondVec(*atm,*bond_neigh);
//								double spatial_bond_x=spatial_bond.GetXCoord();
//
//								double energy_release=0.;
//
//								// \lambda_\mathcal{R}
//								double lambda_release=2.0*(1./interatomic_stretch_0);
//								// \lambda_{\mathcal{R}_0}
//								double lambda_initial_release=1.10865*(1./interatomic_stretch_0);
//
//								double first_term=pow((1./lambda_initial_release),11)-pow((1./lambda_release),11);
//								double first_coeff=(13./11.)*pow(sigma0,12);
//
//								double second_term=pow((1./lambda_release),5)-pow((1./lambda_initial_release),5);
//								double second_coeff=(7./5.)*pow(sigma0,6);
//
//								energy_release=4*epsilon0*(first_term*first_coeff+second_term*second_coeff);
//
//								// \lambda_\mathcal{B}
//								double lambda_breaking=1.75*(1./interatomic_stretch_0);
//								// \lambda_{\mathcal{B}_0}
//								double lambda_initial_breaking=1.223*(1./interatomic_stretch_0);
//
//								double bond_breaking_energy=0.0;
//
//								double first_div_breaking=(1./lambda_breaking);
//								double second_div_breaking=(1./lambda_initial_breaking);
//
//								double first_coeff_breaking=pow(sigma0,12);
//								double first_div_breaking_12=pow(first_div_breaking,12);
//								double second_div_breaking_12=pow(second_div_breaking,12);
//								double first_term_breaking=first_div_breaking_12-second_div_breaking_12;
//
//								double second_coeff_breaking=pow(sigma0,6);
//								double first_div_breaking_6=pow(first_div_breaking,6);
//								double second_div_breaking_6=pow(second_div_breaking,6);
//
//								double second_term_breaking=second_div_breaking_6-first_div_breaking_6;
//
//								bond_breaking_energy=4*epsilon0*(first_coeff_breaking*first_term_breaking-second_coeff_breaking*second_term_breaking);
//
//
//	//    						double bond_breaking_energy=0.;
//	//    						bond_breaking_energy=4*epsilon0*(pow((sigma0/(1.75*(1./interatomic_stretch_0))),12)-pow((sigma0/(1.75*(1./interatomic_stretch_0))),6))-
//	//    								             4*epsilon0*(pow((sigma0/(1.223*(1./interatomic_stretch_0))),12)-pow((sigma0/(1.223*(1./interatomic_stretch_0))),6));
//	//
//	//    						double energy_release=0.;
//	//
//	//    						energy_release=(4*epsilon0*(1./1.9))*
//	//    								       (pow(sigma0*(1./1.9),12)-pow(sigma0*(1./1.9),6))-
//	//    								       (4*epsilon0* (1./1.10865))*
//	//										   (pow(sigma0*(1./1.10865),12)-pow(sigma0*(1./1.10865),6));
//
//								bfile << load_step << " " << bond_neighbor_id << " " << atom_id << " "
//									  << energy_release << " " << bond_breaking_energy << " "<< abs(configForce)
//									  << " "<< initial_bond_x << " "<< material_bond_x << " " << spatial_bond_x
//									  << " " << material_bond_x*initial_bond_x <<endl;
//
//							}
//
//	//    					&& configForce <= abs(criticalConfigForce) && distance_spatial >= critical_stretch
////							spatial_bond_distance >= critical_spatial_distance
//							else if ( configForce >= abs(criticalConfigForce)
//									  && (bond_neighborMaterialPositionX <= crackRegion_MinX  || bond_neighborMaterialPositionY <= crackRegion_MinY
//									  || bond_neighborMaterialPositionX>= crackRegion_MaxX  || bond_neighborMaterialPositionY >= crackRegion_MaxY
//									  || abs(bond_neighborMaterialPositionY-atomMaterialPositionY) < 0.6) )
//							{
//
//								update_bond_neighbor_list.push_back(*bond_neighbor);
//
//							}
//
//							(*bond_neighbor)->SetNeighbors(bond_neighbor_neighbors);
//
//							if (bond_neighbor_neighbors.size()!= bond_neigh_size)
//							{
//
//								cout << "bond_neighbor_neighbors size: " << bond_neighbor_neighbors.size() << " bond_neighbor size: " <<  bond_neigh_size << endl;
//
//							}
//
//						}
//
//						if (update_neighbor_list.size()!= neighbors.size()
//							|| update_bond_neighbor_list.size()!=bond_neighbors.size())
//						{
//
//							(*atom)->SetNeighbors(update_neighbor_list);
//							(*atom)->SetBondNeighbors(update_bond_neighbor_list);
//
//							cout << "update_neighbor_list size: " << update_neighbor_list.size() << " neighbors size: " << neighbors.size() << " update_bond_neighbor_list size: " <<
//									update_bond_neighbor_list.size() << " bond_neighbors size: " << bond_neighbors.size() << endl;
//
//						}
//
//					}
//
//				}
//    		}
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//    	Point <2> *reactionForce=new Point<2>;
//    	reactionForce -> SetXCoord(0.0);
//        reactionForce -> SetYCoord(0.0);
//
//		double stress_11=0.;
//		double stress_12=0.;
//		double stress_21=0.;
//		double stress_22=0.;
//		double vonMisesStress=0.;
//
//        AtomicStress(atoms);
//
//    	if ( file.is_open() || ofile.is_open())
//    	{
//
//    		file <<"ITEM: TIMESTEP" <<endl;
//    		file << load_step <<endl;
//    		file <<"ITEM: NUMBER OF ATOMS" << endl;
//    		file << atoms.size() << endl;
////			file << "2306" << endl;
//    		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
//    		file << "0" << " " << "1.6" << endl;
//    		file << "0" << " " << "1.6" << endl;
//    		file << "0" << " " << "1.6" << endl;
//    		file <<"ITEM: ATOMS id type x y z kx ky kz fx fy fz" <<endl;
//
//    		ofile <<"ITEM: TIMESTEP" <<endl;
//    		ofile << load_step <<endl;
//    		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
//    		ofile << atoms.size() << endl;
////			ofile << "2306" << endl;
//    		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
//    		ofile << "0" << " " << "1.6" << endl;
//    		ofile << "0" << " " << "1.6" << endl;
//    		ofile << "0" << " " << "1.6" << endl;
//    		ofile <<"ITEM: ATOMS id type x y z kx ky kz fx fy fz" <<endl;
//
//    		for (At at=atoms.begin(); at!=atoms.end(); ++at)
//    		{
//
//
//    			if((*at)->GetAtomRegion()!=3 )
//
//    			{
//
//    				Point <2> config_force =force.ConfigurationalForce(*at);
//    				Point <2> interaction_force=force.ResultantForce(*at);
//
//    				Point <2> matposition = (*at)-> GetMaterialPosition();
//    				Point <2> spatialposition = (*at)-> GetSpatialPosition();
//
//    				double matposition_x=matposition.GetXCoord();
//    				double matposition_y=matposition.GetYCoord();
//
//    				int atom_id=(*at)->GetID();
//
//    				vector <double> atomic_stress=(*at)->getAtomicStress();
//
//    				stress_11=atomic_stress[0];
//    				stress_12=atomic_stress[1];
//    				stress_21=atomic_stress[2];
//    				stress_22=atomic_stress[3];
//    				vonMisesStress=atomic_stress[4];
//
//
//
//    				double materialp_x=matposition.GetXCoord();
//    				double materialp_y=matposition.GetYCoord();
//
//    				double spatialp_x=spatialposition.GetXCoord();
//    				double spatialp_y=spatialposition.GetYCoord();
//
//    				double config_force_x=0.;
//    				double config_force_y=0.;
//
//    				double resultant_force_x=0;
//    				double resultant_force_y=0;
//
//    				resultant_force_x=interaction_force.GetXCoord();
//    				resultant_force_y=interaction_force.GetYCoord();
//
//    				config_force_x=config_force.GetXCoord();
//    				config_force_y=config_force.GetYCoord();
//
//    				double AtomPotEng=0.;
//    				AtomPotEng=energy.atomTotalPotentialEnergy(*at);
//
//    				if ((*at)->GetAtomRegion()==4 )
//    				{
//
//    					file << atom_id<< " " << "4" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0"<<endl;
//
//    					ofile << atom_id<< " " << "4" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y<<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==5 )
//    				{
//
//    					file << atom_id<< " " << "5" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0"<<endl;
//
//    					ofile << atom_id<< " " << "5" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==6 )
//    				{
//
//    					file << atom_id<< " " << "6" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    					ofile << atom_id<< " " << "6" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==1 && (*at)->getAtomZone()==1 )
//    				{
//
//    					file << atom_id<< " " << "1" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    					ofile << atom_id<< " " << "1" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==1 && (*at)->getAtomZone()==0 )
//    				{
//
//    					file << atom_id<< " " << "1" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    					ofile << atom_id<< " " << "1" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//
//    				else if ((*at)->GetAtomRegion()==2 && (*at)->getAtomZone()==1 )
//    				{
//
//        				Point <2> atomEmpiricalForce=force.ResultantForce(*at);
//
//        				double Ex=atomEmpiricalForce.GetXCoord();
//        				double Ey=atomEmpiricalForce.GetYCoord();
//
//        				double Rx=reactionForce->GetXCoord();
//        				double Ry=reactionForce->GetYCoord();
//
//        				reactionForce->SetXCoord(Ex+Rx);
//        				reactionForce->SetYCoord(Ey+Ry);
//
//
//    					file << atom_id<< " " << "2" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    					ofile << atom_id<< " " << "2" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==2 && (*at)->getAtomZone()==0 )
//    				{
//
//        				Point <2> atomEmpiricalForce=force.ResultantForce(*at);
//
//        				double Ex=atomEmpiricalForce.GetXCoord();
//        				double Ey=atomEmpiricalForce.GetYCoord();
//
//        				double Rx=reactionForce->GetXCoord();
//        				double Ry=reactionForce->GetYCoord();
//
//        				reactionForce->SetXCoord(Ex+Rx);
//        				reactionForce->SetYCoord(Ey+Ry);
//
//
//    					file << atom_id<< " " << "2" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" << " " <<
//								AtomPotEng<<endl;
//
//    					ofile << atom_id<< " " << "2" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//
//    				else if ((*at)->GetAtomRegion()==0 && (*at)->getAtomZone()==1)
//    				{
//
//                			file << atom_id<< " " << "0" <<" " << materialp_x<<" "
//                					<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//									" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//									<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//                			ofile << atom_id<< " " << "0" <<" " << spatialp_x<<" "
//                					<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//									" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//									<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//    				else if ((*at)->GetAtomRegion()==0 && (*at)->getAtomZone()==0 )
//    				{
//
//    					file << atom_id<< " " << "0" <<" " << materialp_x<<" "
//    							<< setprecision(5)<< materialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    					ofile << atom_id<< " " << "0" <<" " << spatialp_x<<" "
//    							<< setprecision(5)<< spatialp_y <<" " << "0.0"<<
//								" "<< config_force_x << " " << config_force_y <<" " << "0.0"
//								<< " "<< resultant_force_x<< " "<< resultant_force_y<< " " << "0.0" <<endl;
//
//    				}
//
//
//    			}
//
//    		}
//
//    	}
//
//    	file.close();
//    	ofile.close();
//
//    	if (RFfile.is_open() ){
//
//    	double Forcex=reactionForce-> GetXCoord();
//    	double Forcey=reactionForce-> GetYCoord();
//
//    	RFfile << load <<" " << reactionForce->PointNorm() << endl;
//
//    	}
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    	writeDataFile (atoms, load_step, path4);
//
//
//    	double totalEnergy=0.0;
//    	totalEnergy=energy.TotPotentialEnergy(atoms);
//
//    	efile << load_step << " " << totalEnergy << endl;
//
//    	if ( (load_step >= min_load && load_step <= max_load) )
//    	{
//
//    		load+=0.005;
//
//			for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//			{
//				if( (*atom)->GetAtomRegion()==2 )
//				{
//					Point <2> Poldt = (*atom) -> GetMaterialPosition();
//					double Yold = Poldt.GetYCoord();
//					double Xold = Poldt.GetXCoord();
//
//					Point <2> Pnewt;
//					double Ynew= Yold + (load/2);
//					Pnewt.SetXCoord(Xold);
//					Pnewt.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewt);
//
//				}
//
//				else if( (*atom)->GetAtomRegion()==1 )
//				{
//
//					Point <2> Poldb = (*atom) -> GetMaterialPosition();
//					double Yold = Poldb.GetYCoord();
//					double Xold = Poldb.GetXCoord();
//
//					Point <2> Pnewb;
//					double Ynew= Yold - (load/2) ;
//					Pnewb.SetXCoord(Xold);
//					Pnewb.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewb);
//
//				}
//
//				else if ( (*atom) -> GetAtomRegion()==3)
//				{
//					Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//					(*atom) -> SetSpatialPosition(Poldf);
//
//				}
//
//			}
//
//		}
//
//
//    	if ( (load_step >= min_unload && load_step <= max_unload))
//    	{
//
//    		unload+=0.005;
//
//			for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//			{
//				if( (*atom)->GetAtomRegion()==2 )
//				{
//					Point <2> Poldt = (*atom) -> GetMaterialPosition();
//					double Yold = Poldt.GetYCoord();
//					double Xold = Poldt.GetXCoord();
//
//					Point <2> Pnewt;
//					double Ynew= Yold +(load/2) - (unload/2);
//					Pnewt.SetXCoord(Xold);
//					Pnewt.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewt);
//
//				}
//
//				else if( (*atom)->GetAtomRegion()==1 )
//				{
//
//					Point <2> Poldb = (*atom) -> GetMaterialPosition();
//					double Yold = Poldb.GetYCoord();
//					double Xold = Poldb.GetXCoord();
//
//					Point <2> Pnewb;
//					double Ynew= Yold -(load/2) + (unload/2) ;
//					Pnewb.SetXCoord(Xold);
//					Pnewb.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewb);
//
//				}
//
//				else if ( (*atom) -> GetAtomRegion()==3)
//				{
//					Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//					(*atom) -> SetSpatialPosition(Poldf);
//
//				}
//
//			}
//
//		}
//
//
//    	if ((load_step >= min_reload && load_step <= max_reload))
//    	{
//
//    		reload+=0.005;
//
//			for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//			{
//				if( (*atom)->GetAtomRegion()==2 )
//				{
//					Point <2> Poldt = (*atom) -> GetMaterialPosition();
//					double Yold = Poldt.GetYCoord();
//					double Xold = Poldt.GetXCoord();
//
//					Point <2> Pnewt;
//					double Ynew= Yold + (reload/2);
//					Pnewt.SetXCoord(Xold);
//					Pnewt.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewt);
//
//				}
//
//				else if( (*atom)->GetAtomRegion()==1 )
//				{
//
//					Point <2> Poldb = (*atom) -> GetMaterialPosition();
//					double Yold = Poldb.GetYCoord();
//					double Xold = Poldb.GetXCoord();
//
//					Point <2> Pnewb;
//					double Ynew= Yold - (reload/2) ;
//					Pnewb.SetXCoord(Xold);
//					Pnewb.SetYCoord(Ynew);
//
//					(*atom) -> SetSpatialPosition(Pnewb);
//
//				}
//
//				else if ( (*atom) -> GetAtomRegion()==3)
//				{
//					Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//					(*atom) -> SetSpatialPosition(Poldf);
//
//				}
//
//			}
//
//		}
//
//    }
//
//	efile.close();
//	wfile.close();
//	bfile.close();
//	RFfile.close();

	}






