/*
 * optimizer_ceres_3d.cpp
 *
 *  Created on: May 20, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
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


class MolecularStatic : public ceres::FirstOrderFunction
{
private:

	const vector < Atom <3>* > Atoms;
	double sigma;
	double epsilon;
	int dofs;


public:


    MolecularStatic(const vector < Atom <3>* > atoms, double sigma0, double epsilon0, int dofs0):
    	           Atoms(atoms), sigma(sigma0), epsilon(epsilon0), dofs(dofs0){}

    typedef typename vector < Atom <3>* >::const_iterator At;


	virtual bool Evaluate (const double* parameters, double* cost, double* gradient) const
	{

	    Energy <3> energy(1.0, 1.0);
	    Force <3> force(1.0, 1.0);

	    int iter=0;

	    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
	    {

	    	int atom_region=(*atom)->GetAtomRegion();
	    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

	    	{
	    		if (atom_region==6)
	    		{
		    		Point <3> material_position=(*atom)->GetMaterialPosition();
		    		double fixed_x= material_position.GetXCoord() ;

		    		Point <3> spatial_position_3d_region6;

		    		spatial_position_3d_region6.SetXCoord(fixed_x);
		    		spatial_position_3d_region6.SetYCoord(parameters[iter+1]);
		    		spatial_position_3d_region6.SetZCoord(parameters[iter+2]);

		    		(*atom)->SetSpatialPosition(spatial_position_3d_region6);

		    		iter+=3;

	    		}
	    		else
	    		{

	        		Point <3> spatial_position_3d;

	        		spatial_position_3d.SetXCoord(parameters[iter]);
	        		spatial_position_3d.SetYCoord(parameters[iter+1]);
	        		spatial_position_3d.SetZCoord(parameters[iter+2]);

		        	(*atom)->SetSpatialPosition(spatial_position_3d);

		        	iter+=3;

	    		}



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

            		Point <3> resultant_force=force.ResultantForce(*atom);

            		gradient[iter_force]=-1*resultant_force.GetXCoord();
            		gradient[iter_force+1]=-1*resultant_force.GetYCoord();
            		gradient[iter_force+2]=-1*resultant_force.GetZCoord();

            		iter_force+=3;

        		}


        	}

        }

    	cost[0]=energy.TotPotentialEnergy(Atoms);

        return true;

	}

	virtual int NumParameters() const
	{
		return{dofs};

	}

};

int main ()
{

    vector < Atom <3>* > unrelax_atoms;
//    unrelax_atoms=UnrelaxedConfigGenerator <3> (20, 9, 29, 1.6, 3);
//    unrelax_atoms=UnrelaxedConfigGenerator <3> (22, 15, 35, 1.6, 3);
    unrelax_atoms=UnrelaxedConfigGenerator <3> (25, 11, 33, 1.1, 3);
    Boundary <3> boundary_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;
    double bottom_BC_z_min=0.0;

    double bottom_BC_x_max=30;
    double bottom_BC_y_max=10;
    double bottom_BC_z_max=0.1;
    int bottom_BC_id=1;

    boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min, bottom_BC_z_min, bottom_BC_x_max,
    		                          bottom_BC_y_max, bottom_BC_z_max, bottom_BC_id);

    Boundary <3> boundary_top;

    double top_BC_x_min=0.0;
    double top_BC_y_min=0.0;
    double top_BC_z_min=29.5;

    double top_BC_x_max=30;
    double top_BC_y_max=10;
    double top_BC_z_max=30;

    int top_BC_id=2;

    boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min, top_BC_z_min, top_BC_x_max,
    		                       top_BC_y_max, top_BC_z_max, top_BC_id);

    Boundary <3> crack_top;

    double crack_top_x_min=0;
    double crack_top_y_min=0;
    double crack_top_z_min=15.1;

    double crack_top_x_max=10.3;
    double crack_top_y_max=10;
    double crack_top_z_max=19.5;


    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min, crack_top_z_min,
        		                    crack_top_x_max, crack_top_y_max, crack_top_z_max, crack_top_id);

    Boundary <3> crack_bottom;

    double crack_bottom_x_min=0;
    double crack_bottom_y_min=0;
    double crack_bottom_z_min=10.5;

    double crack_bottom_x_max=10.3;
    double crack_bottom_y_max=10;
    double crack_bottom_z_max=14.5;

    int crack_bottom_id=4;


    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min, crack_bottom_z_min,
                                       crack_bottom_x_max, crack_bottom_y_max, crack_bottom_z_max, crack_bottom_id);

    Boundary <3> crack;

    double crack_x_min=0;
    double crack_y_min=0;
    double crack_z_min=14.6;

    double crack_x_max=10.3;
    double crack_y_max=11.3;
    double crack_z_max=15.0;

    int crack_id=3;

    crack.SetBoundaryRegion(crack_x_min, crack_y_min, crack_z_min,crack_x_max, crack_y_max, crack_z_max, crack_id);

    vector < Boundary <3> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);
    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);
    boundaries.push_back(crack);


    typedef typename vector < Atom <3>* >::iterator At;
    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,2.5);

    int id=0;

 	ofstream ufile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/removed atoms/dump.unrelax");

 	if (ufile.is_open())
 	{

		ufile <<"ITEM: TIMESTEP" <<endl;
		ufile << "0" <<endl;
		ufile <<"ITEM: NUMBER OF ATOMS" << endl;
		ufile << atoms.size() << endl;
		ufile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		ufile << "0" << " " << "50" << endl;
		ufile << "0" << " " << "50" << endl;
		ufile << "0" << " " << "50" << endl;
		ufile <<"ITEM: ATOMS id type x y z" <<endl;

 		for (At at=atoms.begin(); at!=atoms.end(); ++at)
 		{

 			if((*at)->GetAtomRegion()!=3)

 			{
     			Point <3> matposition = (*at)-> GetMaterialPosition();

     			double materialp_x=matposition.GetXCoord();
     			double materialp_y=matposition.GetYCoord();
     			double materialp_z=matposition.GetZCoord();

     			if ((*at)->GetAtomRegion()==4)
     			{

         			ufile << id<< " " << "4" <<" " << materialp_x << " "
         		    << setprecision(5)<< materialp_y << " "<< materialp_z << endl;

     			}



     			else if ((*at)->GetAtomRegion()==5)
     			{

         			ufile << id<< " " << "5" <<" " << materialp_x <<" "
         		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

     			}

     			else if ((*at)->GetAtomRegion()==1)
     			{

         			ufile << id<< " " << "1" <<" " << materialp_x << " "
         		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

     			}

     			else if ((*at)->GetAtomRegion()==2)
     			{

         			ufile << id << " " << "2" <<" " << materialp_x <<" "
         		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

     			}

     			else
     			{

         			ufile << id << " " << "0" <<" " << materialp_x <<" "
         		    << setprecision(5)<< materialp_y << " "<< materialp_z << endl;

     			}

     			id+=1;

 			}

 		}

 	}

     else
 	{
 	    cerr << "Couldn't open dump.unrelax for writing." << std::endl;
 	}


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

    double max_load=2.0;
    double load=0;
    int load_step=0;

    while (load < max_load)
    {

	ceres::GradientProblem problem (new MolecularStatic(atoms, 1.0, 1.0,num_optimized_atoms));
	ceres::GradientProblemSolver::Options options;
	options.minimizer_progress_to_stdout=true;
	options.max_num_iterations=2000;
	options.gradient_tolerance=1e-10;
	options.parameter_tolerance=1e-15;
	options.line_search_direction_type=ceres::LBFGS;
	ceres::GradientProblemSolver::Summary summary;
	ceres::Solve(options, problem, parameters, &summary);


	Force <3> force(1.0, 1.0);



    if (load==0)
    {

    	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    	{

    		Point <3> spatialP= (*atom)->GetSpatialPosition();
    		(*atom)->SetMaterialPosition(spatialP);

    	}

    }

    load+=0.05;
    load_step+=1;
    string out_string;
    stringstream ss;
    ss << load_step;
    out_string=ss.str();

    string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/removed atoms/dump.crack";
    string path1="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/spatial/removed atoms/dump.crackspatial" ;

    string file_name=path+out_string;
    string file_name1=path1+out_string;

	ofstream file;
	ofstream ofile;
	file.open(file_name.c_str());
	ofile.open(file_name1.c_str());

	id=0;

	if (file.is_open())
	{

		file <<"ITEM: TIMESTEP" <<endl;
		file << load_step <<endl;
		file <<"ITEM: NUMBER OF ATOMS" << endl;
		file << atoms.size() << endl;
		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		file << "0" << " " << "50" << endl;
		file << "0" << " " << "50" << endl;
		file << "0" << " " << "50" << endl;
		file <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

		ofile <<"ITEM: TIMESTEP" <<endl;
		ofile <<load_step <<endl;
		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
		ofile << atoms.size() << endl;
		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		ofile << "0" << " " << "50" << endl;
		ofile << "0" << " " << "50" << endl;
		ofile << "0" << " " << "50" << endl;
		ofile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

		for (At at=atoms.begin(); at!=atoms.end(); ++at)
		{

			if((*at)->GetAtomRegion()!=3)

			{
				Point <3> config_force =force.ConfigurationalForce(*at);
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

    			if ((*at)->GetAtomRegion()==4)
    			{

        			file << id<< " " << "4" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

        			ofile << id<< " " << "4" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y << " "<< spatialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;


    			}



    			else if ((*at)->GetAtomRegion()==5)
    			{

        			file << id<< " " << "5" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

        			ofile << id<< " " << "5" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y << " "<< spatialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

    			}

    			else if ((*at)->GetAtomRegion()==1)
    			{

        			file << id<< " " << "1" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

        			ofile << id<< " " << "1" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y << " "<< spatialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

    			}


    			else if ((*at)->GetAtomRegion()==2)
    			{

        			file << id<< " " << "2" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

        			ofile << id<< " " << "2" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y << " "<< spatialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

    			}


    			else
    			{

        			file << id<< " " << "0" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;

        			ofile << id<< " " << "0" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y << " "<< spatialp_z <<
				    " "<< config_force_x << " " << config_force_y << " " << config_force_z <<endl;



    			}

    			id+=1;

			}

		}

	}

    else
	{
	    cerr << "Couldn't open dump.crack for writing." << std::endl;
	}

    Point <3> top_force;
    top_force.SetXCoord(0.0);
    top_force.SetYCoord(0.0);
    top_force.SetZCoord(0.0);

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

            Point <3> resultant_force=force.ResultantForce(*atom);
            top_force=top_force+resultant_force;

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

	std::cout << summary.FullReport()<< "\n";

    }

return 0;

}



