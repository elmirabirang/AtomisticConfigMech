/*
 * test_eam_cu.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: S.Elmira Birang.O
 *      Chair of Applied Mechanics
 *      Central Institute of Scientific Computing
 *      University of Erlangen-Nuremberg
 */



#include "atom.h"
#include "point.h"
#include "energy.h"
#include "force.h"
#include "boundary.h"
#include "bond.h"

#include "unrelaxed_config.cpp"
#include "energy.cpp"
#include "force.cpp"
#include "boundary.cpp"
#include "assign_region.cpp"
#include "bond.cpp"
#include "write_data_file.cpp"

#include "broken_bonds.cpp"
#include "config_force_criterion.cpp"

#include <vector>
#include <string>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

#include <eigen3/Eigen/Core>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "gflags/gflags.h"

#include "modification_interaction_list_3d.cpp"

using namespace std;

class MolecularStatic : public ceres::FirstOrderFunction
{

private:

	const vector < Atom <3>* > Atoms;
	int dof;

public:


    MolecularStatic(const vector < Atom <3>* > atoms, int dof0):
    	           Atoms(atoms), dof(dof0){}

    typedef typename vector < Atom <3>* >::const_iterator At;


	virtual bool Evaluate (const double* parameters, double* cost, double* gradient) const
	{

	    Energy <3> energy;
	    Force <3> force;

		double a= 3.80362;

		double beta1= 0.17394;
		double beta2= 5.3566e2;

		double r_03= -2.19885;
		double r_04= -2.61984e2;

		double F0= -2.28235;
		double F2= 1.35535;

		double q1= -1.27775 ;
		double q2=-0.86074;
		double q3= 1.78804;
		double q4= 2.97571;

		double Q1= 0.4 ;
		double Q2= 0.3 ;

        double alpha1=2.97758;
    	double alpha2=1.54927;

        double r_01=0.83591;
    	double r_02=4.46867;

    	double E1=2.01458e2;
    	double E2=6.59288e-3;

    	double delta=0.86225e-2;
    	double h= 0.50037;
    	double r_cut=5.50678;

    	double r_s1=2.24;
    	double r_s2=1.8;
    	double r_s3=1.2;

        double s1=4;
    	double s2=40;
    	double s3=1.150e3;

	    int iter=0;

	    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
	    {

	    	int atom_region=(*atom)->GetAtomRegion();

	    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

	    	{

        		Point <3> spatial_position_3d;

        		spatial_position_3d.SetXCoord(parameters[iter]);
        		spatial_position_3d.SetYCoord(parameters[iter+1]);
        		spatial_position_3d.SetZCoord(parameters[iter+2]);

	        	(*atom)->SetSpatialPosition(spatial_position_3d);

	        	iter+=3;

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

        			Point <3> atom_morse_force=force.derivativeMorseFunction(*atom,alpha1, alpha2,
        		                                                             r_01, r_02, E1, E2, delta,
        		                                                             r_cut, h, r_s1, r_s2, r_s3,
        		                                                             s1, s2, s3);

        			Point <3> atom_embedding_force=force.embeddingForce(*atom, a, beta1,
        			                                                     beta2, r_03, r_04, F0, F2,
        			                                                     q1, q2, q3, q4, Q1, Q2,
        					                                             h, r_cut);

            		gradient[iter_force]=-1*atom_morse_force.GetXCoord()-1*atom_embedding_force.GetXCoord();
            		gradient[iter_force+1]=-1*atom_morse_force.GetYCoord()-1*atom_embedding_force.GetYCoord();
            		gradient[iter_force+2]=-1*atom_morse_force.GetZCoord()-1*atom_embedding_force.GetZCoord();


            		iter_force+=3;

        		}

        	}

        }


    	cost[0]=(1*energy.morseFunction(Atoms, alpha1, alpha2,
                                      r_01, r_02, E1, E2, delta,
                                      r_cut, h, r_s1, r_s2, r_s3,
                                      s1, s2, s3)+energy.embeddingEnergy ( Atoms, a, beta1,
                         	                                               beta2, r_03, r_04, F0, F2,
                         			                                       q1, q2, q3, q4, Q1, Q2,
                         					                                h, r_cut));



    	cout << "eng:" << cost[0] << endl;

        return true;

	}

	virtual int NumParameters() const
	{
		return{dof};

	}

};

int main ()
{

	vector < Atom <3>* > unrelax_atoms;
	unrelax_atoms=UnrelaxedConfigGenerator <3> (20, 9, 39, 3.615, 3);

	Boundary <3> boundary_bottom;
	Boundary <3> boundary_top;
	Boundary <3> crack_top;
	Boundary <3> crack_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;
    double bottom_BC_z_min=0.0;

    double bottom_BC_x_max=70;
    double bottom_BC_y_max=70;
    double bottom_BC_z_max=0.1;

	int bottom_BC_id=1;

	boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min,bottom_BC_z_min, bottom_BC_x_max,
			                          bottom_BC_y_max,bottom_BC_z_max, bottom_BC_id);

    double top_BC_x_min=0.0;
    double top_BC_y_min=0.0;
    double top_BC_z_min=68;

    double top_BC_x_max=70;
    double top_BC_y_max=70;
    double top_BC_z_max=69;

	int top_BC_id=2;

	boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min,top_BC_z_min, top_BC_x_max,
			                          top_BC_y_max,top_BC_z_max, top_BC_id);

    double crack_top_x_min=0;
    double crack_top_y_min=0;
    double crack_top_z_min=9.5*3.615;

    double crack_top_x_max=8*3.615;
    double crack_top_y_max=4.6*3.615;
    double crack_top_z_max=11.5*3.615;
    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min, crack_top_z_min,
        		                crack_top_x_max, crack_top_y_max,crack_top_z_max,crack_top_id);

    double crack_bottom_x_min=0;
    double crack_bottom_y_min=0;
    double crack_bottom_z_min=7.5*3.615;

    double crack_bottom_x_max=8*3.615;
    double crack_bottom_y_max=4.6*3.615;
    double crack_bottom_z_max=9.5*3.615;
    int crack_bottom_id=4;


    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min,crack_bottom_z_min,
                                   crack_bottom_x_max, crack_bottom_y_max,crack_bottom_z_max ,crack_bottom_id);



    vector < Boundary <3> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);
    boundaries.push_back(crack_bottom);
    boundaries.push_back(crack_top);

    typedef typename vector < Atom <3>* >::iterator At;

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {

        AssignRegion(*atom,boundaries);

    }

    vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,3);


    double alpha1=2.97758;
	double alpha2=1.54927;

    double r_01=0.83591;
	double r_02=4.46867;

	double E1=2.01458e2;
	double E2=6.59288e-3;

	double delta=0.86225e-2;
	double h= 0.50037;
	double r_cut=5.50679;

	double r_s1=2.24;
	double r_s2=1.8;
	double r_s3=1.2;

    double s1=4;
	double s2=40;
	double s3=1150;


	double a= 3.80362;

	double beta1= 0.17394;
	double beta2= 5.3566e2;

	double r_03= -2.19885;
	double r_04= -261.984;

	double F0= -2.28235;
	double F2= 1.35535;

	double q1= -1.27775 ;
	double q2=-0.86074;
	double q3= 1.78804;
	double q4= 2.97571;

	double Q1= 0.4 ;
	double Q2= 0.3 ;

	Force <3> forcei;
	double force_norm=0.0;

	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)

	    {

	           Point <3> atom_embedding_force=forcei.embeddingForce(*atom, a, beta1,
                                                        beta2, r_03, r_04, F0, F2,
                                                        q1, q2, q3, q4, Q1, Q2,
                                                       h, r_cut);

	           force_norm=force_norm+atom_embedding_force.PointNorm();
	    }


	cout << "atom embedding force: " << force_norm <<endl;



	vector < Atom <3>* > optimized_atoms;
	int i=0;
    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ( (*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
    	{
            optimized_atoms.push_back(*atom);
    	}

    }

    int num_optimized_atoms=optimized_atoms.size()*3;

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

	cout << "number of optimized atoms: " << num_optimized_atoms << endl;

    string path2="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/screening/irreversible_copper/dump.cracktip";
    ofstream ffile;
	ffile.open(path2.c_str());

    string path3="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/screening/irreversible_copper/dump.totalenergy";
    ofstream Efile;
	Efile.open(path3.c_str());

    string path4="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/screening/irreversible_copper/dump.reactionforce";
    ofstream EFfile;
	EFfile.open(path4.c_str());

    double max_load=24.0;

    double load=0;
    double unload=0;
    double load_secondary=0;
    bool check_broken_bond=false;

    int load_step=0;

    while (load < max_load)
    {

        load+=0.05;
        load_step+=1;

        Force <3> force;

    	ceres::GradientProblem problem (new MolecularStatic(atoms, num_optimized_atoms));
    	ceres::GradientProblemSolver::Options options;
    	options.minimizer_progress_to_stdout=true;
    	options.max_num_iterations=2000;
    	options.gradient_tolerance=1e-6;
    	options.line_search_direction_type=ceres::LBFGS;

//    	options.parameter_tolerance=1e-10;
//    	options.max_lbfgs_rank=20;

    	ceres::GradientProblemSolver::Summary summary;
    	ceres::Solve(options, problem, parameters, &summary);
    	std::cout << summary.FullReport()<< "\n";

//    	string path= "/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/spatial/screening/large_displacement_step/configurational_force/copper/";

//    	writeDataFile (atoms, load_step, path);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        string out_string;
        stringstream ss;
        ss << load_step;
        out_string=ss.str();


        string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/screening/irreversible_copper/dump.crack";
        string path1="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/spatial/screening/irreversible_copper/dump.crackspatial" ;
        string path5="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/screening/irreversible_copper/dump.temporalmaterial";


        string file_name=path+out_string;
        string file_name1=path1+out_string;
        string file_name2=path5+out_string;

    	ofstream file;
    	ofstream ofile;
    	ofstream Mfile;

    	file.open(file_name.c_str());
    	ofile.open(file_name1.c_str());
    	Mfile.open(file_name2.c_str());

    	double id=0;

    	Point <3> *reactionForce=new Point<3>;
    	reactionForce -> SetXCoord(0.0);
    	reactionForce -> SetYCoord(0.0);
    	reactionForce -> SetZCoord(0.0);

    	if ( file.is_open() && ffile.is_open() && Efile.is_open() )
    	{

    		file <<"ITEM: TIMESTEP" <<endl;
    		file << load_step <<endl;
    		file <<"ITEM: NUMBER OF ATOMS" << endl;
    		file << atoms.size() << endl;
    		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
    		file << "0" << " " << "1.6" << endl;
    		file << "0" << " " << "1.6" << endl;
    		file << "0" << " " << "1.6" << endl;
    		file <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

    		ofile <<"ITEM: TIMESTEP" <<endl;
    		ofile <<load_step <<endl;
    		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
    		ofile << atoms.size() << endl;
    		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
    		ofile << "0" << " " << "1.6" << endl;
    		ofile << "0" << " " << "1.6" << endl;
    		ofile << "0" << " " << "1.6" << endl;
    		ofile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

    		Mfile <<"ITEM: TIMESTEP" <<endl;
    		Mfile <<load_step <<endl;
    		Mfile <<"ITEM: NUMBER OF ATOMS" << endl;
    		Mfile << atoms.size() << endl;
    		Mfile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
    		Mfile << "0" << " " << "1.6" << endl;
    		Mfile << "0" << " " << "1.6" << endl;
    		Mfile << "0" << " " << "1.6" << endl;
    		Mfile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

    		Energy <3> energy;
    		Force <3> force;
    		double totalEnergy=energy.TotPotentialEnergy(atoms);

    		Efile << load << " " << totalEnergy << endl;



    		for (At at=atoms.begin(); at!=atoms.end(); ++at)
    		{

    			if ((*at)->GetAtomRegion()!=3)

    			{

    				Point <3> config_force_pair = force.configForceEamPair(*at, alpha1, alpha2,
                                                                           r_01, r_02, E1, E2, delta,
                                                                           r_cut, h, r_s1, r_s2, r_s3,
                                                                           s1, s2, s3);

    				Point <3> config_force_embedding = force.configForceEamEmbedding(*at, a, beta1,
                                                                                     beta2, r_03, r_04, F0, F2,
                                                                                     q1, q2, q3, q4, Q1, Q2,
                                                                                     h, r_cut);

    				Point <3>  config_force=(config_force_pair+config_force_embedding)*-1;

        			Point <3> matposition = (*at)-> GetMaterialPosition();
        			Point <3> spatialposition = (*at)-> GetSpatialPosition();

        			double materialp_x=matposition.GetXCoord();
        			double materialp_y=matposition.GetYCoord();
        			double materialp_z=matposition.GetZCoord();

        			double spatialp_x=spatialposition.GetXCoord();
        			double spatialp_y=spatialposition.GetYCoord();
        			double spatialp_z=spatialposition.GetZCoord();

        			double config_force_x=config_force.GetXCoord();
        			double config_force_y=config_force.GetYCoord();
        			double config_force_z=config_force.GetZCoord();

        			double temporalp_x=materialp_x-0.05*config_force_x;
        			double temporalp_y=materialp_y-0.05*config_force_y;
        			double temporalp_z=materialp_z-0.05*config_force_z;

        			int atom_id=(*at)->GetID();


        			if ((*at)->GetAtomRegion()==4)
        			{

            			file << id<< " " << "4" <<" " << materialp_x<<" "
            		    << setprecision(5)<< materialp_y <<" " <<  materialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			ofile << id<< " " << "4" <<" " << spatialp_x<<" "
            		    << setprecision(5)<< spatialp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			Mfile << id<< " " << "4" <<" " << temporalp_x<<" "
            		    << setprecision(5)<< temporalp_y <<" " << temporalp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

        			}



        			else if ((*at)->GetAtomRegion()==5)
        			{

            			file << id<< " " << "5" <<" " << materialp_x<<" "
            		    << setprecision(5)<< materialp_y <<" " << materialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			ofile << id<< " " << "5" <<" " << spatialp_x<<" "
            		    << setprecision(5)<< spatialp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			Mfile << id<< " " << "5" <<" " << temporalp_x<<" "
            		    << setprecision(5)<< temporalp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

        			}

        			else if ((*at)->GetAtomRegion()==1)
        			{

            			file << id<< " " << "1" <<" " << materialp_x<<" "
            		    << setprecision(5)<< materialp_y <<" " << materialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			ofile << id<< " " << "1" <<" " << spatialp_x<<" "
            		    << setprecision(5)<< spatialp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y<<" " << config_force_z <<endl;

            			Mfile << id<< " " << "1" <<" " << temporalp_x<<" "
            		    << setprecision(5)<< temporalp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

        			}


        			else if ((*at)->GetAtomRegion()==2)
        			{



            			Point <3> atom_morse_force=force.derivativeMorseFunction(*at,alpha1, alpha2,
            		                                                             r_01, r_02, E1, E2, delta,
            		                                                             r_cut, h, r_s1, r_s2, r_s3,
            		                                                             s1, s2, s3);

            			Point <3> atom_embedding_force=force.embeddingForce(*at, a, beta1,
            			                                                     beta2, r_03, r_04, F0, F2,
            			                                                     q1, q2, q3, q4, Q1, Q2,
            					                                             h, r_cut);


        				double Ex=(atom_morse_force+atom_embedding_force).GetXCoord();
        				double Ey=(atom_morse_force+atom_embedding_force).GetYCoord();
        				double Ez=(atom_morse_force+atom_embedding_force).GetZCoord();

        				double Rx=reactionForce->GetXCoord();
        				double Ry=reactionForce->GetYCoord();
        				double Rz=reactionForce->GetZCoord();

        				reactionForce->SetXCoord(Ex+Rx);
        				reactionForce->SetYCoord(Ey+Ry);
        				reactionForce->SetZCoord(Ez+Rz);


            			file << id<< " " << "2" <<" " << materialp_x<<" "
            		    << setprecision(5)<< materialp_y <<" " << materialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			ofile << id<< " " << "2" <<" " << spatialp_x<<" "
            		    << setprecision(5)<< spatialp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			Mfile << id<< " " << "2" <<" " << temporalp_x<<" "
            		    << setprecision(5)<< temporalp_y <<" " << temporalp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

        			}

    				else
    				{

            			file << id<< " " << "0" <<" " << materialp_x<<" "
            		    << setprecision(5)<< materialp_y <<" " << materialp_z<<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			ofile << id<< " " << "0" <<" " << spatialp_x<<" "
            		    << setprecision(5)<< spatialp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

            			Mfile << id<< " " << "0" <<" " << temporalp_x<<" "
            		    << setprecision(5)<< temporalp_y <<" " << spatialp_z <<
    				    " "<< config_force_x << " " << config_force_y <<" " << config_force_z <<endl;

    				}


    			}
    			id+=1;

    		}

    	}

//    	double crack_line=33.5;
//
//    	ModifyNeighborList <3> (atoms, alpha1, alpha2,
//    	                   r_01, r_02, E1, E2, delta,
//    	                   r_cut, h, r_s1, r_s2, r_s3,
//    	                   s1, s2, s3, a, beta1,
//    					   beta2, r_03, r_04, F0, F2,
//    					   q1, q2, q3, q4, Q1, Q2, crack_line);

		string write_data_path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests"
				               "/results/dim3/spatial/screening/irreversible_copper/trajectory/trajectory.data";

		writeDataFile <3> (atoms, load_step, write_data_path);

		if (EFfile.is_open() ){

			double Rforce_x=reactionForce->GetXCoord();
			double Rforce_y=reactionForce->GetYCoord();
			double Rforce_z=reactionForce->GetZCoord();

		    EFfile << load<<" " << reactionForce->PointNorm() << endl;

		}

	    else
		{
		    cerr << "Couldn't open dump.reactionforce for writing." << std::endl;
		}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //loading

    	if (load <= 8.0)
    	{

    		if (load_step >120)

    			{

    				check_broken_bond=ConfigForceCriterion <3> (atoms);

    			}

            for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
            {
            	if( (*atom)->GetAtomRegion()==2 )
            	{
                	Point <3> Poldt = (*atom) -> GetMaterialPosition();
                	double Yold = Poldt.GetYCoord();
                	double Xold = Poldt.GetXCoord();
                	double Zold=Poldt.GetZCoord();

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
                    double Zold=Poldb.GetZCoord();

                    Point <3> Pnewb;
                    double Znew= Zold - (load/2) ;

                    Pnewb.SetXCoord(Xold);
                    Pnewb.SetYCoord(Yold);
                    Pnewb.SetZCoord(Znew);

                    (*atom) -> SetSpatialPosition(Pnewb);

            	}

            	else if ( (*atom) -> GetAtomRegion()==3)
            	{

            		Point <3> Poldf= (*atom) -> GetSpatialPosition();

            		(*atom) -> SetSpatialPosition(Poldf);

            	}


            }


    	}

    	//unloading

    	else if (load > 8.0 && load <= 16.0)
    	{

    		unload=0.05;

            for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
            {
            	if( (*atom)->GetAtomRegion()==2 )
            	{
                	Point <3> Poldt = (*atom) -> GetSpatialPosition();
                	double Yold = Poldt.GetYCoord();
                	double Xold = Poldt.GetXCoord();
                	double Zold=Poldt.GetZCoord();

                	Point <3> Pnewt;

                	double Znew= Zold - (unload/2);
                	Pnewt.SetXCoord(Xold);
                	Pnewt.SetYCoord(Yold);
                	Pnewt.SetZCoord(Znew);

                    (*atom) -> SetSpatialPosition(Pnewt);

            	}

            	else if( (*atom)->GetAtomRegion()==1 )
            	{

                    Point <3> Poldb = (*atom) -> GetSpatialPosition();

                    double Yold = Poldb.GetYCoord();
                    double Xold = Poldb.GetXCoord();
                    double Zold=Poldb.GetZCoord();

                    Point <3> Pnewb;
                    double Znew= Zold + (unload/2) ;

                    Pnewb.SetXCoord(Xold);
                    Pnewb.SetYCoord(Yold);
                    Pnewb.SetZCoord(Znew);

                    (*atom) -> SetSpatialPosition(Pnewb);

            	}

            	else if ( (*atom) -> GetAtomRegion()==3)
            	{

            		Point <3> Poldf= (*atom) -> GetSpatialPosition();

            		(*atom) -> SetSpatialPosition(Poldf);

            	}


            }
    	}

    	//reloading

    	else if (load > 16.0 && load <= 24.0)
    	{

    		load_secondary=0.05;

            for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
            {
            	if( (*atom)->GetAtomRegion()==2 )
            	{
                	Point <3> Poldt = (*atom) -> GetSpatialPosition();
                	double Yold = Poldt.GetYCoord();
                	double Xold = Poldt.GetXCoord();
                	double Zold=Poldt.GetZCoord();

                	Point <3> Pnewt;

                	double Znew= Zold + (load_secondary/2);
                	Pnewt.SetXCoord(Xold);
                	Pnewt.SetYCoord(Yold);
                	Pnewt.SetZCoord(Znew);

                    (*atom) -> SetSpatialPosition(Pnewt);

            	}

            	else if( (*atom)->GetAtomRegion()==1 )
            	{

                    Point <3> Poldb = (*atom) -> GetSpatialPosition();

                    double Yold = Poldb.GetYCoord();
                    double Xold = Poldb.GetXCoord();
                    double Zold=Poldb.GetZCoord();

                    Point <3> Pnewb;
                    double Znew= Zold - (load_secondary/2) ;

                    Pnewb.SetXCoord(Xold);
                    Pnewb.SetYCoord(Yold);
                    Pnewb.SetZCoord(Znew);

                    (*atom) -> SetSpatialPosition(Pnewb);

            	}

            	else if ( (*atom) -> GetAtomRegion()==3)
            	{

            		Point <3> Poldf= (*atom) -> GetSpatialPosition();

            		(*atom) -> SetSpatialPosition(Poldf);

            	}


            }

       }

    }

    ffile.close();
    Efile.close();
    EFfile.close();

}




