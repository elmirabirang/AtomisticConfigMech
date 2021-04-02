/**
 * test_performance_LJ.cpp
 *
 *  Created on: Oct 26, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
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

#include <eigen3/Eigen/Core>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "gflags/gflags.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;
#include <algorithm>

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

//    auto start0 = high_resolution_clock::now();

	virtual bool Evaluate (const double* parameters, double* cost, double* gradient) const
	{



	    Energy <2> energy(1.0, 1.0);
	    Force <2> force(1.0, 1.0);

	    auto start1 = high_resolution_clock::now();

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

	    auto stop1 = high_resolution_clock::now();
	    auto duration1 =duration_cast<microseconds>(stop1 - start1);

	    cout << "Time taken by function1: " << duration1.count() << " microseconds" << std::endl;

	    auto start2 = high_resolution_clock::now();

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

	    auto stop2 = high_resolution_clock::now();
	    auto duration2 =duration_cast<microseconds>(stop2 - start2);

	    cout << "Time taken by function2: " << duration2.count() << " microseconds" << std::endl;

	    auto start3 = high_resolution_clock::now();

    	cost[0]=energy.TotPotentialEnergy(Atoms);

	    auto stop3 = high_resolution_clock::now();
	    auto duration3 =duration_cast<microseconds>(stop3 - start3);

	    cout << "Time taken by function3: " << duration3.count() << " microseconds" << std::endl;

        return true;

	}

	virtual int NumParameters() const
	{
		return{dof};

	}

};

int main()
{

    vector < Atom <2>* > unrelax_atoms;

//    auto start = high_resolution_clock::now();

    unrelax_atoms=UnrelaxedConfigGenerator <2> (53, 49, 0, 1.1, 2, "model");

//    auto stop = high_resolution_clock::now();
//    auto duration =duration_cast<microseconds>(stop - start);
//
//    cout << "Time taken by function: "
//    		<< duration.count() << " microseconds" << std::endl;

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

    vector < Boundary <2> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);

    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);

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

	auto start0 = high_resolution_clock::now();

	ceres::GradientProblem problem (new MolecularStatic(atoms, 1.0, 1.0, num_optimized_atoms));
	ceres::GradientProblemSolver::Options options;
	options.minimizer_progress_to_stdout=true;
	options.max_num_iterations=2000;
	options.gradient_tolerance=1e-10;
	options.line_search_direction_type=ceres::LBFGS;
	ceres::GradientProblemSolver::Summary summary;
	ceres::Solve(options, problem, parameters, &summary);

    auto stop0 = high_resolution_clock::now();
    auto duration0 =duration_cast<microseconds>(stop0 - start0);

    cout << "Time taken by function0: " << duration0.count() << " microseconds" << std::endl;

}
