/*
 * calculate_critical_EAM_force.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nuremberg
 *
 *
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

using namespace std;

int main ()
{

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

		Force <2> force;
		Energy <2> energy;
		double atom_force_norm=0.0;
		typedef typename vector < Atom <2>* >::const_iterator At;

	    string path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/eam_force.txt";
	    ofstream file;
		file.open(path.c_str());

		double max_interatomic_distance=2.4;
		double interatomic_distance=2.25;
		vector < Atom <2>* > unrelax_atoms;

		vector <double> interatomic_distance_vector;

		while (interatomic_distance < max_interatomic_distance)

		{
			interatomic_distance_vector.push_back(interatomic_distance);
			interatomic_distance+=0.005;

        }


            for (int i=0; i<interatomic_distance_vector.size() ;++i)
            {
             	double distance=interatomic_distance_vector[i];

            	cout << distance << endl;

    			unrelax_atoms=UnrelaxedConfigGenerator <2> (2, 1, 0, distance, 2);
    			vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,10);


	    	for(At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	    	{

    			    Point <2> atom_morse_force(0.,0.);
    			    Point <2> atom_embedding_force(0.,0.);

	    			atom_morse_force=force.derivativeMorseFunction(*atom,alpha1, alpha2,
	    		                                                             r_01, r_02, E1, E2, delta,
	    		                                                             r_cut, h, r_s1, r_s2, r_s3,
	    		                                                             s1, s2, s3);



	    			atom_embedding_force=force.embeddingForce(*atom, a, beta1,
	    			                                                     beta2, r_03, r_04, F0, F2,
	    			                                                     q1, q2, q3, q4, Q1, Q2,
	    					                                             h, r_cut);

	    			Point <2>  atom_force(0.,0.);

            		atom_force.SetXCoord(atom_morse_force.GetXCoord()+atom_embedding_force.GetXCoord());
            		atom_force.SetYCoord(atom_morse_force.GetYCoord()+atom_embedding_force.GetYCoord());
            		atom_force.SetZCoord(atom_morse_force.GetZCoord()+atom_embedding_force.GetZCoord());

            		double force_length=0;
            		force_length=atom_force.PointNorm();

            		Point <2> materialposition=(*atom)->GetMaterialPosition();


    				Point <2> config_force_pair = force.configForceEamPair(*atom, alpha1, alpha2,
                                                                           r_01, r_02, E1, E2, delta,
                                                                           r_cut, h, r_s1, r_s2, r_s3,
                                                                           s1, s2, s3);

    				Point <2> config_force_embedding = force.configForceEamEmbedding(*atom, a, beta1,
                                                                                     beta2, r_03, r_04, F0, F2,
                                                                                     q1, q2, q3, q4, Q1, Q2,
                                                                                     h, r_cut);

    				Point <2>  config_force=(config_force_pair+config_force_embedding)*-1;

        			double config_force_x=config_force.GetXCoord();
        			double config_force_y=config_force.GetYCoord();
        			double config_force_z=config_force.GetZCoord();

        			if ((*atom)->GetID()==1)
        			{
        			     file << distance << " " << -1*atom_force.GetXCoord() << " Config Force: " << config_force_x <<endl;

        			}


	    	}

    		double Eng=(1*energy.morseFunction(atoms, alpha1, alpha2,
    		                                      r_01, r_02, E1, E2, delta,
    		                                      r_cut, h, r_s1, r_s2, r_s3,
    		                                      s1, s2, s3)+energy.embeddingEnergy ( atoms, a, beta1,
    		                         	                                               beta2, r_03, r_04, F0, F2,
    		                         			                                       q1, q2, q3, q4, Q1, Q2,
    		                         					                                h, r_cut));

//    		file << distance << " " << Eng <<endl;

            }


}


