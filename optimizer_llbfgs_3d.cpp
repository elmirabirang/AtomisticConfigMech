/*
 * optimizer_llbfgs_3d.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nuremberg
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


#include <lbfgs.h>



using namespace std;

struct InputData
{
    vector < Atom <3>* > Atoms;
    double sigma;
    double epsilon;

};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{

	InputData data_tmp = *((InputData *)instance);


    double sigma=data_tmp.sigma;
    double epsilon=data_tmp.epsilon;


    vector < Atom <3>* > Atoms;
    Atoms=data_tmp.Atoms;

    Energy <3> energy(sigma,epsilon);
    Force <3> force(sigma,epsilon);

    int dimension=3;

    typedef typename vector < Atom <3>* >::iterator At;

    switch (dimension){

    case 2:{

    int iter=0;
    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
    {
    	//only for atoms that are labeled

    	int atom_region=(*atom)->GetAtomRegion();


    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

    	{

        	Point <3> spatial_position;
        	spatial_position.SetXCoord(x[iter]);
        	spatial_position.SetYCoord(x[iter+1]);

        	(*atom)->SetSpatialPosition(spatial_position);

        	iter+=2;

    	}

    }
    break;
    }

    case 3:{

        int iter=0;
    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{
    		int atom_region=(*atom)->GetAtomRegion();
    		if( atom_region!=1 && atom_region!=2 && atom_region!=3 )
    		{
        		Point <3> spatial_position_3d;

        		spatial_position_3d.SetXCoord(x[iter]);
        		spatial_position_3d.SetYCoord(x[iter+1]);
        		spatial_position_3d.SetZCoord(x[iter+2]);

        		(*atom)->SetSpatialPosition(spatial_position_3d);

        		iter+=3;

    		}


    	}

    	break;
    }

    }



    switch (dimension)
    {


    case 2:{

    int iter=0;
    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
    {
    	//only for atoms that are labeled

    	int atom_region=(*atom)->GetAtomRegion();

    	if(atom_region!=1 && atom_region!=2 && atom_region!=3)

    	{

        	Point <3> resultant_force=force.ResultantForce(*atom);

        	lbfgsfloatval_t gx=-1*resultant_force.GetXCoord();
        	lbfgsfloatval_t gy=-1*resultant_force.GetYCoord();

        	g[iter]=gx;
        	g[iter+1]=gy;

        	iter+=2;

    	}

    }
    break;
    }


    case 3:{

        int iter_force=0;
    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{

    		int atom_region=(*atom)->GetAtomRegion();
    		if( atom_region!=1 && atom_region!=2 && atom_region!=3 )
    		{

        		Point <3> resultant_force=force.ResultantForce(*atom);

        		lbfgsfloatval_t gx=-1*resultant_force.GetXCoord();
        		lbfgsfloatval_t gy=-1*resultant_force.GetYCoord();
        		lbfgsfloatval_t gz=-1*resultant_force.GetZCoord();

        		g[iter_force]=gx;

        		g[iter_force+1]=gy;

        		g[iter_force+2]=gz;

        		iter_force+=3;

    		}


    	}
    	break;
    }

    }


    double eng=energy.TotPotentialEnergy(Atoms);
    lbfgsfloatval_t energy2d=eng;


    cout << "energy: " << eng <<endl;

    return(energy2d);

}

static int progress (
		void *instance,
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
		)
{
	ofstream file("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/energy_force.tex");

	file << "energy: " << fx << endl;
	file << "force: " << gnorm << endl;

    printf("Iteration %d:\n", k);
    printf("  fx = %f= ", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;

}


int main(){

//update the material position inside load step loop

    vector < Atom <3>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <3> (32, 9, 27, 1.1, 3);

    cout << "number of atoms: "<< unrelax_atoms.size() <<endl;

    typedef typename vector < Atom <3>* >::iterator At;
    //send material_position to real_1d_array x

    Boundary <3> boundary_bottom;

    double bottom_BC_x_min=0.0;
    double bottom_BC_y_min=0.0;
    double bottom_BC_z_min=0.0;

    double bottom_BC_x_max=44;
    double bottom_BC_y_max=9;
    double bottom_BC_z_max=0.1;
    int bottom_BC_id=1;

    boundary_bottom.SetBoundaryRegion(bottom_BC_x_min, bottom_BC_y_min, bottom_BC_z_min, bottom_BC_x_max,
    		                          bottom_BC_y_max, bottom_BC_z_max, bottom_BC_id);

    Boundary <3> boundary_top;

    double top_BC_x_min=0.0;
    double top_BC_y_min=0.0;
    double top_BC_z_min=24;

    double top_BC_x_max=44;
    double top_BC_y_max=9;
    double top_BC_z_max=25;
    int top_BC_id=2;

    boundary_top.SetBoundaryRegion(top_BC_x_min, top_BC_y_min, top_BC_z_min, top_BC_x_max,
    		                       top_BC_y_max, top_BC_z_max, top_BC_id);

    Boundary <3> crack_top;

    double crack_top_x_min=0;
    double crack_top_y_min=0;
    double crack_top_z_min=12.2;

    double crack_top_x_max=12.3;
    double crack_top_y_max=9;
    double crack_top_z_max=16.5;

    int crack_top_id=5;

    crack_top.SetBoundaryRegion(crack_top_x_min, crack_top_y_min, crack_top_z_min,
    		                    crack_top_x_max, crack_top_y_max, crack_top_z_max, crack_top_id);

    Boundary <3> crack_bottom;

    double crack_bottom_x_min=0;
    double crack_bottom_y_min=0;
    double crack_bottom_z_min=8.5;

    double crack_bottom_x_max=12.3;
    double crack_bottom_y_max=9;
    double crack_bottom_z_max=12.16;

    int crack_bottom_id=4;


    crack_bottom.SetBoundaryRegion(crack_bottom_x_min, crack_bottom_y_min, crack_bottom_z_min,
                                   crack_bottom_x_max, crack_bottom_y_max, crack_bottom_z_max, crack_bottom_id);

//    Boundary <3> crack;

//    double xa=0;
//    double ya=0;
//    double za=5.6;
//
//    double xb=2.5;
//    double yb=8.5;
//    double zb=5.62;
//
//    crack.SetBoundaryRegion(xa,ya,za,xb,yb,zb,3);


    vector < Boundary <3> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);

//    boundaries.push_back(crack);


    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,2.5);

    vector < Atom <3>* > fixed_atoms;
    vector < Atom <3>* > bottom_atoms;
    vector < Atom <3>* > top_atoms;
    vector < Atom <3>* > optimized_atoms;

    int i, ret=0;
    int N=(7655-568)*3;

    lbfgsfloatval_t *x=lbfgs_malloc(N);
    lbfgsfloatval_t fx;
    lbfgs_parameter_t param;

    int count=0;
    int count_bottom=0;
    int count_top=0;

    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ( (*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
    	{
            Point <3> mposition;

            mposition=(*atom)->GetMaterialPosition();

            x[i]=mposition.GetXCoord();
            x[i+1]=mposition.GetYCoord();
            x[i+2]=mposition.GetZCoord();

            optimized_atoms.push_back(*atom);
            i+=3;
    	}

    	else if((*atom) -> GetAtomRegion()==1)
    	{
    		bottom_atoms.push_back(*atom);
    		count_bottom+=1;
    	}

    	else if ((*atom) -> GetAtomRegion()==2)
    	{

    		top_atoms.push_back(*atom);
    		count_top+=1;
    	}

    	else if ((*atom)-> GetAtomRegion()==3)
    	{

    		fixed_atoms.push_back(*atom);
    		count+=1;

    	}

    }

    cout << "count: " << count << endl;
    cout << "count bottom: " << count_bottom << endl;
    cout << "count top: " << count_top << endl;

    int id=0;

	ofstream ofile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/dump.unrelax");

	if (ofile.is_open())
	{

		ofile <<"ITEM: TIMESTEP" <<endl;
		ofile <<"0" <<endl;
		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
		ofile << atoms.size() << endl;
		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;

		ofile <<"ITEM: ATOMS id type x y z" <<endl;

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

        			ofile << id<< " " << "4" <<" " << materialp_x << " "
        		    << setprecision(5)<< materialp_y << " "<< materialp_z << endl;

    			}



    			else if ((*at)->GetAtomRegion()==5)
    			{

        			ofile << id<< " " << "5" <<" " << materialp_x <<" "
        		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

    			}

    			else if ((*at)->GetAtomRegion()==1)
    			{

        			ofile << id<< " " << "1" <<" " << materialp_x << " "
        		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

    			}

    			else if ((*at)->GetAtomRegion()==2)
    			{

        			ofile << id << " " << "2" <<" " << materialp_x <<" "
        		    << setprecision(5) << materialp_y << " "<< materialp_z << endl;

    			}


    			else
    			{

        			ofile << id << " " << "0" <<" " << materialp_x <<" "
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


	ofstream lfile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/crack.data");
	id=1;

	if (lfile.is_open())
	{

		for (At at=atoms.begin(); at!=atoms.end(); ++at)
		{

			if((*at)->GetAtomRegion()!=3)

			{

    			Point <3> matposition = (*at)-> GetMaterialPosition();

    			double materialp_x=matposition.GetXCoord();
    			double materialp_y=matposition.GetYCoord();
    			double materialp_z=matposition.GetZCoord();

    			lfile << id << " "<< id<< " " << "1"<< " " << "1.0" <<" "<< materialp_x << " "
        		      << setprecision(5)<< materialp_y << " "<< materialp_z << endl;

    			id+=1;


			}

		}

	}

    else
	{
	    cerr << "Couldn't open crack.data for writing." << std::endl;
	}


    double max_load=2.5;
    double load=0;
    int load_step=0;

    while (load < max_load)
    {

        InputData data;
        data.Atoms = atoms;
        data.sigma=1.0;
        data.epsilon=1.0;

        Force <3> force(data.sigma,data.epsilon);

        lbfgs_parameter_init(&param);
        param.epsilon=1e-6;

        ret=lbfgs(N,x,&fx,evaluate, progress, (void*)&data, &param);

        string file_name=std::to_string(load);
        int id=0;
        load_step+=1;

        if (load==0)
        {

        	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
        	{

        		Point <3> spatialP= (*atom)->GetSpatialPosition();
        		(*atom)->SetMaterialPosition(spatialP);

        	}

        }


    	ofstream file("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/material/dump.crack"+file_name);
    	ofstream ofile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/spatial/dump.crackspatial"+file_name);

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




        load+=0.05;
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

//                (*atom) -> SetMaterialPosition(Pnewt);
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

//                (*atom) -> SetMaterialPosition(Pnewb);
                (*atom) -> SetSpatialPosition(Pnewb);

        	}

        	else if ( (*atom) -> GetAtomRegion()==3)
        	{
        		Point <3> Poldf= (*atom) -> GetSpatialPosition();

        		(*atom) -> SetSpatialPosition(Poldf);

        	}


        }

        double top_forces_norm=top_force.PointNorm();
        file << top_forces_norm << endl;

    }

    return 0;


    }














