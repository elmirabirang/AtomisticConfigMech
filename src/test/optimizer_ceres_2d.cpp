/*
 * optimizer_ceres_2d.cpp
 *
 *  Created on: May 22, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
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
#include "write_data_file.cpp"

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

// save the configurational forces of atom 704 (693)


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

int main ()
{

    vector < Atom <2>* > unrelax_atoms;
//    unrelax_atoms=UnrelaxedConfigGenerator <3> (20, 9, 29, 1.6, 3);
    unrelax_atoms=UnrelaxedConfigGenerator <2> (50, 29, 0, 1.1, 2);
    Boundary <2> boundary_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;

    double bottom_BC_x_max=55.0;
    double bottom_BC_y_max=0.1;
    int bottom_BC_id=1;

    boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min, bottom_BC_x_max,
    		                          bottom_BC_y_max, bottom_BC_id);

    Boundary <2> boundary_top;

    double top_BC_x_min=0.0;
    double top_BC_y_min=26.6;

    double top_BC_x_max=55.0;
    double top_BC_y_max=26.7;
    int top_BC_id=2;

    boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min,top_BC_x_max,
    		                       top_BC_y_max, top_BC_id);

    Boundary <2> crack_top;

    double crack_top_x_min=0;
    double crack_top_y_min=14.0;

    double crack_top_x_max=12.0;
    double crack_top_y_max=16.2;
    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min,
        		               crack_top_x_max, crack_top_y_max,crack_top_id);

    Boundary <2> crack_bottom;

    double crack_bottom_x_min=0;
    double crack_bottom_y_min=10.0;

    double crack_bottom_x_max=12.0;
    double crack_bottom_y_max=13.3;
    int crack_bottom_id=4;


    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min,
                                   crack_bottom_x_max, crack_bottom_y_max, crack_bottom_id);

    Boundary <2> crack;

    double crack_x_min=0;
    double crack_y_min=7*1.1*sqrt(3)-0.01;

    double crack_x_max=12.0;
    double crack_y_max=7*1.1*sqrt(3)+0.01;
    int crack_id=3;

    crack.SetBoundaryRegion(crack_x_min, crack_y_min, crack_x_max, crack_y_max,crack_id);


    vector < Boundary <2> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);
    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);
    boundaries.push_back(crack);

    typedef typename vector < Atom <2>* >::iterator At;
    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,1.5);

    int id=0;

 	ofstream ufile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/dump.unrelax");

 	if (ufile.is_open())
 	{

 		ufile <<"ITEM: TIMESTEP" <<endl;
 		ufile <<"0" <<endl;
 		ufile <<"ITEM: NUMBER OF ATOMS" << endl;
// 		ufile << atoms.size() << endl;
 		ufile << "1425" << endl;
 		ufile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;

		ufile << "0" << " " << "1.6" << endl;
		ufile << "0" << " " << "1.6" << endl;
		ufile << "0" << " " << "1.6" << endl;

 		ufile <<"ITEM: ATOMS id type x y" <<endl;

 		for (At at=atoms.begin(); at!=atoms.end(); ++at)
 		{

 			if((*at)->GetAtomRegion()!=3)

 			{
     			Point <2> matposition = (*at)-> GetMaterialPosition();

     			double materialp_x=matposition.GetXCoord();
     			double materialp_y=matposition.GetYCoord();

     			if ((*at)->GetAtomRegion()==4)
     			{

         			ufile << id<< " " << "4" <<" " << materialp_x << " "
         		    << setprecision(5)<< materialp_y << " " << "0.0"<< endl;

     			}



     			else if ((*at)->GetAtomRegion()==5)
     			{

         			ufile << id<< " " << "5" <<" " << materialp_x <<" "
         		    << setprecision(5) << materialp_y << " " << "0.0" << endl;

     			}

     			else if ((*at)->GetAtomRegion()==1)
     			{

         			ufile << id<< " " << "1" <<" " << materialp_x << " "
         		    << setprecision(5) << materialp_y << " " << "0.0"<< endl;

     			}

     			else if ((*at)->GetAtomRegion()==2)
     			{

         			ufile << id << " " << "2" <<" " << materialp_x <<" "
         		    << setprecision(5) << materialp_y << " " << "0.0"<<  endl;

     			}


     			else
     			{

         			ufile << id << " " << "0" <<" " << materialp_x <<" "
         		    << setprecision(5)<< materialp_y << " " << "0.0"<<  endl;

     			}

     			id+=1;

 			}

 		}

 	}

     else
 	{
 	    cerr << "Couldn't open dump.unrelax for writing." << std::endl;
 	}


	vector < Atom <2>* > optimized_atoms;
	int i=0;
    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ( (*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
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

    cout << "number of optimized atoms: " << num_optimized_atoms << endl;

    double max_load=7.5;
    double load=0;
    double unload=0;
    double load_secondary=0;
    int load_step=0;
    int max_load_step=35;

    string path2="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/blunting/large_displacement_step/configurational_force/dump.cracktip";
    ofstream ffile;
	ffile.open(path2.c_str());

    string path3="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/blunting/large_displacement_step/configurational_force/dump.totalenergy";
    ofstream Efile;
	Efile.open(path3.c_str());

    string path4="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/blunting/large_displacement_step/configurational_force/dump.reactionforce";
    ofstream EFfile;
	EFfile.open(path4.c_str());

    while (load_step < max_load_step)
    {

	ceres::GradientProblem problem (new MolecularStatic(atoms, 1.0, 1.0, num_optimized_atoms));
	ceres::GradientProblemSolver::Options options;
	options.minimizer_progress_to_stdout=true;
	options.max_num_iterations=1000;
	options.gradient_tolerance=1e-6;
	options.line_search_direction_type=ceres::LBFGS;
	ceres::GradientProblemSolver::Summary summary;
	ceres::Solve(options, problem, parameters, &summary);

	Force <2> force(1.0, 1.0);



    if (load==0)
    {

    	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    	{

    		Point <2> spatialP= (*atom)->GetSpatialPosition();
    		(*atom)->SetMaterialPosition(spatialP);

    	}

    }


    load_step+=1;
    string out_string;
    stringstream ss;
    ss << load_step;
    out_string=ss.str();


    string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/blunting/large_displacement_step/configurational_force/dump.crack";
    string path1="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/spatial/blunting/large_displacement_step/configurational_force/dump.crackspatial" ;
    string path5="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/blunting/large_displacement_step/configurational_force/dump.temporalmaterial";


    string file_name=path+out_string;
    string file_name1=path1+out_string;
    string file_name2=path5+out_string;

	ofstream file;
	ofstream ofile;
	ofstream Mfile;


	file.open(file_name.c_str());
	ofile.open(file_name1.c_str());
	Mfile.open(file_name2.c_str());

	id=0;

	Point <2> *reactionForce=new Point<2>;
	reactionForce -> SetXCoord(0.0);
	reactionForce -> SetYCoord(0.0);


	if ( file.is_open() && ffile.is_open() && Efile.is_open() )
	{

		file <<"ITEM: TIMESTEP" <<endl;
		file << load_step <<endl;
		file <<"ITEM: NUMBER OF ATOMS" << endl;
//		file << atoms.size() << endl;
		file << "1425" << endl;
		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		file << "0" << " " << "1.6" << endl;
		file << "0" << " " << "1.6" << endl;
		file << "0" << " " << "1.6" << endl;
		file <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

		ofile <<"ITEM: TIMESTEP" <<endl;
		ofile <<load_step <<endl;
		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
//		ofile << atoms.size() << endl;
		ofile << "1425" << endl;
		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		ofile << "0" << " " << "1.6" << endl;
		ofile << "0" << " " << "1.6" << endl;
		ofile << "0" << " " << "1.6" << endl;
		ofile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

		Mfile <<"ITEM: TIMESTEP" <<endl;
		Mfile <<load_step <<endl;
		Mfile <<"ITEM: NUMBER OF ATOMS" << endl;
//		ofile << atoms.size() << endl;
		Mfile << "1425" << endl;
		Mfile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
		Mfile << "0" << " " << "1.6" << endl;
		Mfile << "0" << " " << "1.6" << endl;
		Mfile << "0" << " " << "1.6" << endl;
		Mfile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;

		Energy <2> energy(1.0, 1.0);
		Force <2> force(1.0, 1.0);
		double totalEnergy=energy.TotPotentialEnergy(atoms);
		Efile << load << " " << totalEnergy << endl;


		for (At at=atoms.begin(); at!=atoms.end(); ++at)
		{

			if((*at)->GetAtomRegion()!=3)

			{
				Point <2> config_force =force.ConfigurationalForce(*at);
    			Point <2> matposition = (*at)-> GetMaterialPosition();
    			Point <2> spatialposition = (*at)-> GetSpatialPosition();

    			double materialp_x=matposition.GetXCoord();
    			double materialp_y=matposition.GetYCoord();


    			double spatialp_x=spatialposition.GetXCoord();
    			double spatialp_y=spatialposition.GetYCoord();


    			double config_force_x=config_force.GetXCoord();
    			double config_force_y=config_force.GetYCoord();

    			double temporalp_x=materialp_x-0.05*config_force_x;
    			double temporalp_y=materialp_y-0.05*config_force_y;

    			int atom_id=(*at)->GetID();

    			if (atom_id==704)
    			{

    				double crack_tip_potential_energy=energy.atomTotalPotentialEnergy(*at);

    				double force_norm=sqrt(pow(config_force_x,2) + pow (config_force_y,2));

    				ffile << load << " "  << force_norm << " " << crack_tip_potential_energy << endl;

    			}


    			if ((*at)->GetAtomRegion()==4)
    			{

        			file << id<< " " << "4" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			ofile << id<< " " << "4" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			Mfile << id<< " " << "4" <<" " << temporalp_x<<" "
        		    << setprecision(5)<< temporalp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

    			}



    			else if ((*at)->GetAtomRegion()==5)
    			{

        			file << id<< " " << "5" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			ofile << id<< " " << "5" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			Mfile << id<< " " << "5" <<" " << temporalp_x<<" "
        		    << setprecision(5)<< temporalp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

    			}

    			else if ((*at)->GetAtomRegion()==1)
    			{

        			file << id<< " " << "1" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y<<" " << "0.0" <<endl;

        			ofile << id<< " " << "1" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y<<" " << "0.0" <<endl;

        			Mfile << id<< " " << "1" <<" " << temporalp_x<<" "
        		    << setprecision(5)<< temporalp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

    			}


    			else if ((*at)->GetAtomRegion()==2)
    			{


    				Point <2> atomEmpiricalForce=force.ResultantForce(*at);

    				double Ex=atomEmpiricalForce.GetXCoord();
    				double Ey=atomEmpiricalForce.GetYCoord();

    				double Rx=reactionForce->GetXCoord();
    				double Ry=reactionForce->GetYCoord();

    				reactionForce->SetXCoord(Ex+Rx);
    				reactionForce->SetYCoord(Ey+Ry);


        			file << id<< " " << "2" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0" <<endl;

        			ofile << id<< " " << "2" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			Mfile << id<< " " << "2" <<" " << temporalp_x<<" "
        		    << setprecision(5)<< temporalp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

    			}


    			else
    			{

        			file << id<< " " << "0" <<" " << materialp_x<<" "
        		    << setprecision(5)<< materialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			ofile << id<< " " << "0" <<" " << spatialp_x<<" "
        		    << setprecision(5)<< spatialp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;

        			Mfile << id<< " " << "0" <<" " << temporalp_x<<" "
        		    << setprecision(5)<< temporalp_y <<" " << "0.0"<<
				    " "<< config_force_x << " " << config_force_y <<" " << "0.0"<<endl;



    			}

    			id+=1;

			}

		}

	}



    else
	{
	    cerr << "Couldn't open dump.crack for writing." << std::endl;
	}

	//writeDataFile <2> (atoms, load_step);

	if (EFfile.is_open() ){

	double Forcex=reactionForce-> GetXCoord();
	double Forcey=reactionForce-> GetYCoord();

	EFfile << load <<" " << reactionForce->PointNorm() << endl;

	}


    else
	{
	    cerr << "Couldn't open dump.reactionforce for writing." << std::endl;
	}



	if (load_step < 50){

		load+=0.05;


	    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	    {
	    	if( (*atom)->GetAtomRegion()==2 )
	    	{
	        	Point <2> Poldt = (*atom) -> GetMaterialPosition();
	        	double Yold = Poldt.GetYCoord();
	        	double Xold = Poldt.GetXCoord();

	        	Point <2> Pnewt;
	        	double Ynew= Yold + ( load/2);
	        	Pnewt.SetXCoord(Xold);
	        	Pnewt.SetYCoord(Ynew);

	            (*atom) -> SetSpatialPosition(Pnewt);
	            Point <2> resultant_force=force.ResultantForce(*atom);

	    	}

	    	else if( (*atom)->GetAtomRegion()==1 )
	    	{

	            Point <2> Poldb = (*atom) -> GetMaterialPosition();
	            double Yold = Poldb.GetYCoord();
	            double Xold = Poldb.GetXCoord();

	            Point <2> Pnewb;
	            double Ynew= Yold - (load/2) ;
	            Pnewb.SetXCoord(Xold);
	            Pnewb.SetYCoord(Ynew);

	            (*atom) -> SetSpatialPosition(Pnewb);

	    	}

	    	else if ( (*atom) -> GetAtomRegion()==3)
	    	{
	    		Point <2> Poldf= (*atom) -> GetSpatialPosition();

	    		(*atom) -> SetSpatialPosition(Poldf);

	    	}

	    }


	}

//	else if (load_step >49 && load_step <100){
//
//		unload=0.05;
//
//	    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//	    {
//	    	if( (*atom)->GetAtomRegion()==2 )
//	    	{
//	        	Point <2> Poldt = (*atom) -> GetSpatialPosition();
//	        	double Yold = Poldt.GetYCoord();
//	        	double Xold = Poldt.GetXCoord();
//
//	        	Point <2> Pnewt;
//	        	double Ynew= Yold - ( unload/2);
//	        	Pnewt.SetXCoord(Xold);
//	        	Pnewt.SetYCoord(Ynew);
//
//	            (*atom) -> SetSpatialPosition(Pnewt);
//	            Point <2> resultant_force=force.ResultantForce(*atom);
//
//	    	}
//
//	    	else if( (*atom)->GetAtomRegion()==1 )
//	    	{
//
//	            Point <2> Poldb = (*atom) -> GetSpatialPosition();
//	            double Yold = Poldb.GetYCoord();
//	            double Xold = Poldb.GetXCoord();
//
//	            Point <2> Pnewb;
//	            double Ynew= Yold + (unload/2) ;
//	            Pnewb.SetXCoord(Xold);
//	            Pnewb.SetYCoord(Ynew);
//
//	            (*atom) -> SetSpatialPosition(Pnewb);
//
//	    	}
//
//	    	else if ( (*atom) -> GetAtomRegion()==3)
//	    	{
//	    		Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//	    		(*atom) -> SetSpatialPosition(Poldf);
//
//	    	}
//
//	    }
//
//	}

//	else if (load_step >99 && load_step <150){
//
//		load_secondary=0.05;
//
//	    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//	    {
//	    	if( (*atom)->GetAtomRegion()==2 )
//	    	{
//	        	Point <2> Poldt = (*atom) -> GetSpatialPosition();
//	        	double Yold = Poldt.GetYCoord();
//	        	double Xold = Poldt.GetXCoord();
//
//	        	Point <2> Pnewt;
//	        	double Ynew= Yold + (load_secondary/2);
//	        	Pnewt.SetXCoord(Xold);
//	        	Pnewt.SetYCoord(Ynew);
//
//	            (*atom) -> SetSpatialPosition(Pnewt);
//	            Point <2> resultant_force=force.ResultantForce(*atom);
//
//	    	}
//
//	    	else if( (*atom)->GetAtomRegion()==1 )
//	    	{
//
//	            Point <2> Poldb = (*atom) -> GetSpatialPosition();
//	            double Yold = Poldb.GetYCoord();
//	            double Xold = Poldb.GetXCoord();
//
//	            Point <2> Pnewb;
//	            double Ynew= Yold - (load_secondary/2) ;
//	            Pnewb.SetXCoord(Xold);
//	            Pnewb.SetYCoord(Ynew);
//
//	            (*atom) -> SetSpatialPosition(Pnewb);
//
//	    	}
//
//	    	else if ( (*atom) -> GetAtomRegion()==3)
//	    	{
//	    		Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//	    		(*atom) -> SetSpatialPosition(Poldf);
//
//	    	}
//
//	    }
//
//	}


	std::cout << summary.FullReport()<< "\n";

    }

ffile.close();
Efile.close();
EFfile.close();


return 0;

}






