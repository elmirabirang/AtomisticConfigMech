/*
 * test_silicon.cpp
 *
 *  Created on: May 30, 2020
 *      Author: S.Elmira Birang.O
 *      Central Institute for Scientific Computing
 *
 *
 * Stillinger-Weber parameters for various elements and mixtures
 *  these entries are in LAMMPS "metal" units:
 *  epsilon = eV; sigma = Angstroms
 *  other quantities are unitless

 *  format of a single entry (one or more lines):
 *  element 1, element 2, element 3,
 *  epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol

 *   Here are the original parameters in metal units, for Silicon from:

 *   Stillinger and Weber,  Phys. Rev. B, v. 31, p. 5262, (1985)


 *   Si Si Si 2.1683  2.0951  1.80  21  1.20  -0.333333333333
              7.049556277  0.6022245584  4.0  0.0 0.0
*/

#include "atom.h"
#include "point.h"
#include "energy.h"
#include "force.h"
#include "boundary.h"
#include "bond.h"
#include "matrix.h"

#include "unrelaxed_config.cpp"
#include "energy.cpp"
#include "force.cpp"
#include "boundary.cpp"
#include "assign_region.cpp"
#include "bond.cpp"
#include "write_data_file.cpp"
#include "matrix.cpp"

#include <math.h>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <list>
#include <algorithm>

#include <eigen3/Eigen/Core>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "gflags/gflags.h"


class MolecularStatic : public ceres::FirstOrderFunction
{

private:

	const vector < Atom <3>* > Atoms;

    double sigma_AlphaBeta, sigma_AlphaGamma,
	       epsilon_AlphaBetaGamma,Gamma, lambda_AlphaBetaGamma
	       ,cosine_teta0, a_AlphaBeta, a_AlphaGamma,
	       A_AlphaBeta, B_AlphaBeta, p_AlphaBeta, q_AlphaBeta;

	int dof;

public:


    MolecularStatic(const vector < Atom <3>* > atoms, double sigma_AlphaBeta0,double sigma_AlphaGamma0,double epsilon_AlphaBetaGamma0,
    		                                          double Gamma0, double lambda_AlphaBetaGamma0
 	                                                , double cosine_teta0_0, double a_AlphaBeta0, double a_AlphaGamma0,
													  double A_AlphaBeta0, double B_AlphaBeta0, double p_AlphaBeta0, double q_AlphaBeta0, int dof0):

    	                                      Atoms(atoms), sigma_AlphaBeta(sigma_AlphaBeta0), sigma_AlphaGamma(sigma_AlphaGamma0),
				                              epsilon_AlphaBetaGamma(epsilon_AlphaBetaGamma0),Gamma(Gamma0), lambda_AlphaBetaGamma(lambda_AlphaBetaGamma0)
			                                  , cosine_teta0(cosine_teta0_0), a_AlphaBeta(a_AlphaBeta0), a_AlphaGamma(a_AlphaGamma0),
			                                  A_AlphaBeta(A_AlphaBeta0), B_AlphaBeta(B_AlphaBeta0), p_AlphaBeta(p_AlphaBeta0), q_AlphaBeta(q_AlphaBeta0), dof(dof0){}

    typedef typename vector < Atom <3>* >::const_iterator At;


	virtual bool Evaluate (const double* parameters, double* cost, double* gradient) const
	{

	    Energy <3> energy;
	    Force <3> force;

	    int iter=0;

	    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
	    {

	    	int atom_region=(*atom)->GetAtomRegion();
	    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

	    	{

        		Point <3> spatial_position_3d(0.,0.,0.);

        		spatial_position_3d.SetXCoord(parameters[iter]);
        		spatial_position_3d.SetYCoord(parameters[iter+1]);
        		spatial_position_3d.SetZCoord(parameters[iter+2]);

	        	(*atom)->SetSpatialPosition(spatial_position_3d);

	        	iter+=3;

	    	}

	    }

		for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
		{


			 force.ResultantSWThreeBodyForce(*atom,sigma_AlphaBeta, sigma_AlphaGamma, epsilon_AlphaBetaGamma,
						                           Gamma, lambda_AlphaBetaGamma
					                               ,cosine_teta0, a_AlphaBeta, a_AlphaGamma);


		}

	    typedef typename vector < Atom <3>* >::const_reverse_iterator A;

        int iter_force=0;
        if (gradient!=NULL)
        {

        	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
        	{

        		int atom_region=(*atom)->GetAtomRegion();
        		if( atom_region!=1 && atom_region!=2 && atom_region!=3 )
        		{

        			Point <3> resultant_force2body(0.,0.,0.);

            		resultant_force2body=force.ResultantSWTwoBodyForce(*atom,sigma_AlphaBeta, epsilon_AlphaBetaGamma,
                                                                                 A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
                                                                                 q_AlphaBeta, a_AlphaBeta);


                    Point <3> atomThreeBodyForce(0.,0.,0.);


                    atomThreeBodyForce=(*atom)-> GetForce();

            		gradient[iter_force]=-resultant_force2body.GetXCoord()-atomThreeBodyForce.GetXCoord();
            		gradient[iter_force+1]=-resultant_force2body.GetYCoord()-atomThreeBodyForce.GetYCoord();
            		gradient[iter_force+2]=-resultant_force2body.GetZCoord()-atomThreeBodyForce.GetZCoord();

            		iter_force+=3;



        		}


        	}

        }

		for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
		{

			Point <3> zero_force(0.,0.,0.);

			(*atom)-> SetForce(zero_force);


		}

        double total_energy=0.;

    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{

    		total_energy=total_energy+energy.totalStillingerWeberEnergy(*atom, sigma_AlphaBeta,
                                                            sigma_AlphaGamma, epsilon_AlphaBetaGamma,
					                                        Gamma, lambda_AlphaBetaGamma,
					                                        cosine_teta0, a_AlphaBeta, a_AlphaGamma,
					                                        A_AlphaBeta, B_AlphaBeta, p_AlphaBeta,
					                                        q_AlphaBeta);

    	}

    	cost[0]=total_energy;

        return true;

	}

	virtual int NumParameters() const
	{
		return{dof};

	}

};


int main()
{
    using namespace std;
    vector < Atom <3>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <3> (1, 20, 10, 5.4310, 3, "Si");


    typedef typename vector < Atom <3>* >::iterator At;
    Energy <3> energy;
    Force  <3> force;

    Boundary <3> boundary_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;
    double bottom_BC_z_min=0.0;

    double bottom_BC_x_max=110.;
    double bottom_BC_y_max=110.;
    double bottom_BC_z_max=0.1;

    int bottom_BC_id=1;

    boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min, bottom_BC_z_min,
    		                          bottom_BC_x_max, bottom_BC_y_max, bottom_BC_z_max,
									  bottom_BC_id);

    Boundary <3> boundary_top;

    double top_BC_x_min=0.0;
    double top_BC_y_min=0.0;
    double top_BC_z_min=54.;

    double top_BC_x_max=110.0;
    double top_BC_y_max=110.;
    double top_BC_z_max=55.0;

    int top_BC_id=2;

    boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min, top_BC_z_min,
    		                       top_BC_x_max, top_BC_y_max, top_BC_z_max, top_BC_id);

    Boundary <3> crack_top;

    double crack_top_x_min=0.;
    double crack_top_y_min=0.;
    double crack_top_z_min=28.4;

    double crack_top_x_max=11.;
    double crack_top_y_max=25.8;
    double crack_top_z_max=28.7;


    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min,crack_top_z_min,
        		                crack_top_x_max, crack_top_y_max,crack_top_z_max, crack_top_id);

    Boundary <3> crack_bottom;

    double crack_bottom_x_min=0.;
    double crack_bottom_y_min=0.0;
    double crack_bottom_z_min=27.;

    double crack_bottom_x_max=11.;
    double crack_bottom_y_max=27.2;
    double crack_bottom_z_max=27.2;

    int crack_bottom_id=4;

    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min, crack_bottom_z_min,
                                   crack_bottom_x_max, crack_bottom_y_max,crack_bottom_z_max, crack_bottom_id);


    vector < Boundary <3> > boundaries;
    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);

    boundaries.push_back(crack_bottom);
    boundaries.push_back(crack_top);

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,2.6);

    string path4="/calculate/elmira/atomistic_config_mech/results/dump.silicon" ;

    writeDataFile (atoms, 1, path4);

	vector < Atom <3>* > optimized_atoms;
	int i=0;
    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ( (*atom) -> GetAtomRegion()==0 || (*atom) -> GetAtomRegion()==4 || (*atom) -> GetAtomRegion()==5 )
    	{
            optimized_atoms.push_back(*atom);
    	}

    }

    int num_optimized_atoms=optimized_atoms.size()*3;

    cout << "num_optimized_atoms: " << num_optimized_atoms << endl;

	int iterator=0;
	int atom_num=0;
	double parameters [num_optimized_atoms];
	while (iterator<num_optimized_atoms)
	{

		Point <3> optimized_atom=optimized_atoms[atom_num]->GetMaterialPosition();

		double x=optimized_atom.GetXCoord();
		double y=optimized_atom.GetYCoord();
		double z=optimized_atom.GetZCoord();

		parameters[iterator]=x;
		parameters[iterator+1]=y;
		parameters[iterator+2]=z;

		iterator+=3;
		atom_num+=1;

	}

    double tot_energy=0;
    int number_load_steps=1;
    int load_step=0;
    double load=0.0;

	ceres::GradientProblem problem (new MolecularStatic(atoms,2.0951, 2.0951,
														2.1683, 1.2, 21.,-0.333333,
														1.8, 1.8,
														7.049556277, 0.6022245584,
														4.0, 0.0, num_optimized_atoms));

	ceres::GradientProblemSolver::Options options;
	options.minimizer_progress_to_stdout=true;
//		options.max_lbfgs_rank=20;
	options.max_num_iterations=2000;
	options.gradient_tolerance=1e-6;
	options.line_search_direction_type=ceres::LBFGS;
//	options.line_search_type=ceres::WOLFE;
	ceres::GradientProblemSolver::Summary summary;
	ceres::Solve(options, problem, parameters, &summary);
	std::cout << summary.FullReport()<< "\n";

    while (load_step< number_load_steps)
    {

		string path2="/calculate/elmira/Implementations/atomistic_config_mech/results/dump.Si";
		ofstream file;
		file.open(path2.c_str());


    	ofstream MFile;
    	ofstream SFile;

    	string MPath="/calculate/elmira/Implementations/atomistic_config_mech/results/three_dim/silicon/material/dump.material";
    	string SPath="/calculate/elmira/Implementations/atomistic_config_mech/results/three_dim/silicon/spatial/dump.spatial";

    	load_step+=1;

    	string out_string;
    	stringstream ss;
    	ss << load_step;
    	out_string=ss.str();

    	string Mfile_name=MPath+out_string;
    	string Sfile_name=SPath+out_string;

    	MFile.open(Mfile_name.c_str());
    	SFile.open(Sfile_name.c_str());

    	if (load_step==1)
    	{
    		for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    		{

    			Point <3> spatialP= (*atom)->GetSpatialPosition();
    			(*atom)->SetMaterialPosition(spatialP);

    		}

    	}


		if ( MFile.is_open() && SFile.is_open())
		{

			MFile <<"ITEM: TIMESTEP" <<endl;
			MFile << load_step <<endl;
			MFile <<"ITEM: NUMBER OF ATOMS" << endl;
			MFile << atoms.size() << endl;
			MFile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
			MFile << "0" << " " << "1.6" << endl;
			MFile << "0" << " " << "1.6" << endl;
			MFile << "0" << " " << "1.6" << endl;
			MFile <<"ITEM: ATOMS id type x y z fx fy fz Kx Ky Kz" <<endl;

			SFile <<"ITEM: TIMESTEP" <<endl;
			SFile << load_step <<endl;
			SFile <<"ITEM: NUMBER OF ATOMS" << endl;
			SFile << atoms.size() << endl;
			SFile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
			SFile << "0" << " " << "1.6" << endl;
			SFile << "0" << " " << "1.6" << endl;
			SFile << "0" << " " << "1.6" << endl;
			SFile <<"ITEM: ATOMS id type x y z fx fy fz Kx Ky Kz" <<endl;


			for(At atom=atoms.begin(); atom!=atoms.end(); ++atom)
			{

				Point <3> spatial_position=(*atom)->GetSpatialPosition();
				Point <3> material_position=(*atom)->GetMaterialPosition();

				double spatialp_x=spatial_position.GetXCoord();
				double spatialp_y=spatial_position.GetYCoord();
				double spatialp_z=spatial_position.GetZCoord();

				double materialp_x=material_position.GetXCoord();
				double materialp_y=material_position.GetYCoord();
				double materialp_z=material_position.GetZCoord();

				Point <3> atomic_force=force.ResultantSWTwoBodyForce((*atom), 2.0951, 2.1683,
																	 7.049556277, 0.6022245584,
																	 4.0, 0.0, 1.8);

				Point <3> atomic_config_force=force.ResultantConfigSWTwoBodyForce((*atom), 2.0951, 2.1683,
																	        7.049556277, 0.6022245584,
																	        4.0, 0.0, 1.8);

				Point <3> atomic_force3body=force.ResultantSWThreeBodyForce((*atom),2.0951, 2.0951,
																			 2.1683, 1.2, 21.,-0.333333,
																			 1.8, 1.8);

//				cout << "atomic_force3body: " << atomic_force3body.PointNorm() << "\n";

				Point <3> atomic_config_force3body=force.ResultantConfigSWThreeBodyForce((*atom),2.0951, 2.0951,
																			              2.1683, 1.2, 21.,-0.333333,
																			              1.8, 1.8);

				double energy_atom_a=0.;

				energy_atom_a=energy.totalStillingerWeberEnergy((*atom),
																2.0951, 2.0951,
																2.1683, 1.2, 21.,-0.33,
																1.8, 1.8,
																7.049556277, 0.6022245584,
																4.0, 0.0);

				tot_energy+=energy_atom_a;

				double atomic_force_x=atomic_force.GetXCoord();
				double atomic_force_y=atomic_force.GetYCoord();
				double atomic_force_z=atomic_force.GetZCoord();

				double atomic_config_force_x=atomic_config_force.GetXCoord();
				double atomic_config_force_y=atomic_config_force.GetYCoord();
				double atomic_config_force_z=atomic_config_force.GetZCoord();

				double atomic_force3body_x=atomic_force3body.GetXCoord();
				double atomic_force3body_y=atomic_force3body.GetYCoord();
				double atomic_force3body_z=atomic_force3body.GetZCoord();

				double atomic_config_force3body_x=atomic_config_force3body.GetXCoord();
				double atomic_config_force3body_y=atomic_config_force3body.GetYCoord();
				double atomic_config_force3body_z=atomic_config_force3body.GetZCoord();

				int atom_id=(*atom)->GetID();


				if ((*atom)->GetAtomRegion()==0)
				{

				SFile << atom_id << " " << "0" <<" " << spatialp_x<<" "
					 << setprecision(5)<< spatialp_y <<" " << spatialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;

				MFile << atom_id << " " << "0" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << materialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;
				}

				else if ((*atom)->GetAtomRegion()==1)
				{

				SFile << atom_id << " " << "1" <<" " << spatialp_x<<" "
					 << setprecision(5)<< spatialp_y <<" " << spatialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;

				MFile << atom_id << " " << "1" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << materialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;
				}

				else if ((*atom)->GetAtomRegion()==2)
				{

				SFile << atom_id << " " << "2" <<" " << spatialp_x<<" "
					 << setprecision(5)<< spatialp_y <<" " << spatialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;

				MFile << atom_id << " " << "2" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << materialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;
				}

				else if ((*atom)->GetAtomRegion()==4)
				{

				SFile << atom_id << " " << "4" <<" " << spatialp_x<<" "
					 << setprecision(5)<< spatialp_y <<" " << spatialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;

				MFile << atom_id << " " << "4" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << materialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;
				}

				else
				{

				SFile << atom_id << " " << "5" <<" " << spatialp_x<<" "
					 << setprecision(5)<< spatialp_y <<" " << spatialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;

				MFile << atom_id << " " << "5" <<" " << materialp_x<<" "
					 << setprecision(5)<< materialp_y <<" " << materialp_z<<
					 " "<< atomic_force_x+atomic_force3body_x<< " "
					 << atomic_force_y+atomic_force3body_y<< " " << atomic_force_z+atomic_force3body_z
					 << " " << atomic_config_force_x+atomic_config_force3body_x
					 << " " << atomic_config_force_y+atomic_config_force3body_y
					 << " " << atomic_config_force_z+atomic_config_force3body_z<<endl;
				}
			}

		}

		MFile.close();
		SFile.close();

		load+=0.001;

    	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    	{
    		if( (*atom)->GetAtomRegion()==2 )
    		{
    			Point <3> Poldt = (*atom) -> GetMaterialPosition();
    			double Yold = Poldt.GetYCoord();
    			double Xold = Poldt.GetXCoord();
    			double Zold = Poldt.GetZCoord();

    			Point <3> Pnewt;
    			double Znew= Zold + (load/2);
    			Pnewt.SetXCoord(Xold);
    			Pnewt.SetYCoord(Yold);
    			Pnewt.SetZCoord(Znew);

    			(*atom) -> SetSpatialPosition(Pnewt);

    		}

    		else if( (*atom)->GetAtomRegion()==1 )
    		{

    			Point <3> Poldb = (*atom) -> GetMaterialPosition();
    			double Yold = Poldb.GetYCoord();
    			double Xold = Poldb.GetXCoord();
    			double Zold = Poldb.GetZCoord();

    			Point <3> Pnewb;
    			double Znew= Zold - (load/2) ;
    			Pnewb.SetXCoord(Xold);
    			Pnewb.SetYCoord(Yold);
    			Pnewb.SetZCoord(Znew);


    			(*atom) -> SetSpatialPosition(Pnewb);

    		}

    	}

		ceres::GradientProblem problem (new MolecularStatic(atoms,2.0951, 2.0951,
															2.1683, 1.2, 21.,-0.333333333,
															1.8, 1.8,
															7.049556277, 0.6022245584,
															4.0, 0.0, num_optimized_atoms));

		ceres::GradientProblemSolver::Options options;
		options.minimizer_progress_to_stdout=true;
		options.max_num_iterations=2000;
		options.gradient_tolerance=1e-6;
		options.line_search_direction_type=ceres::LBFGS;
//		options.line_search_type=ceres::WOLFE;
//		options.max_num_line_search_step_size_iterations=40;
		ceres::GradientProblemSolver::Summary summary;
		ceres::Solve(options, problem, parameters, &summary);
		std::cout << summary.FullReport()<< "\n";


    }


    cout << "total energy: " << tot_energy   << endl;

    return 0;
}


