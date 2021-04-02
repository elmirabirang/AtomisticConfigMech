/*
 * optimizer.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
 */


#include "alglib/optimization.h"
#include "alglib/optimization.cpp"
#include "atom.h"
#include "point.h"
#include "energy.h"
#include "force.h"
#include "boundary.h"

#include "unrelaxed_config.cpp"
#include "energy.cpp"
#include "force.cpp"
#include "boundary.cpp"
#include "assign_region.cpp"

#include <vector>
#include <string>

#include "alglib/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "alglib/alglibinternal.cpp"
#include "alglib/alglibmisc.cpp"
#include "alglib/ap.cpp"
#include "alglib/dataanalysis.cpp"
#include "alglib/diffequations.cpp"
#include "alglib/fasttransforms.cpp"
#include "alglib/integration.cpp"
#include "alglib/interpolation.cpp"
#include "alglib/linalg.cpp"
#include "alglib/solvers.cpp"
#include "alglib/specialfunctions.cpp"
#include "alglib/statistics.cpp"


#include "alglib/alglibinternal.h"
#include "alglib/alglibmisc.h"
#include "alglib/ap.h"
#include "alglib/dataanalysis.h"
#include "alglib/diffequations.h"
#include "alglib/fasttransforms.h"
#include "alglib/integration.h"
#include "alglib/interpolation.h"
#include "alglib/linalg.h"
#include "alglib/solvers.h"
#include "alglib/specialfunctions.h"
#include "alglib/statistics.h"

#include <nlopt.hpp>
#include <string>

using namespace std;
using namespace alglib;

struct InputData
{
    vector < Atom <2>* > Atoms;
    double sigma;
    double epsilon;

};

double optimize_atomistic(const std::vector<double> &material_position, std::vector<double> &grad, void *my_func_data)
{

//	for (int i=0; i< x.size(); ++i)
//	{
//		cout << x[i] <<"\n";
//	}


    InputData data_tmp = *((InputData *)my_func_data);


    double sigma=data_tmp.sigma;
    double epsilon=data_tmp.epsilon;

    vector < Atom <2>* > Atoms;
    Atoms=data_tmp.Atoms;
    typedef typename vector < Atom <2>* >::iterator At;

    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    {
    	Point <2> position=(*atom) -> GetSpatialPosition();
    	double px=position.GetXCoord();
    	double py=position.GetYCoord();

    	cout << "x: " << px << " y: "<< py << "\n";

    	Point <2> positionm=(*atom)-> GetMaterialPosition();
    	double pmx=positionm.GetXCoord();
    	double pmy=positionm.GetYCoord();

    	cout << "xm: " << pmx << " ym: "<<  pmy << "\n";

    }

    Energy <2> energy(sigma,epsilon);
    Force <2> force(sigma,epsilon);

    int dimension=2;



    if (dimension==2){

    int iter=0;
    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
    {
    	//only for atoms that are labeled

    	int atom_id=(*atom)->GetID();

//    	&& atom_id!=603 && atom_id!=604 && atom_id!=605
//    	        	&& atom_id!=606 && atom_id!=607 && atom_id!=608 && atom_id!=609 && atom_id!=610 && atom_id!=611
//    	    		&& atom_id!=612 && atom_id!=613 && atom_id!=614 && atom_id!=615 && atom_id!=616 && atom_id!=617


    	if( atom_id > 2 && atom_id < 5 )

    	{

        	Point <2> spatial_position;
        	spatial_position.SetXCoord(material_position[iter]);
        	spatial_position.SetYCoord(material_position[iter+1]);

        	(*atom)->SetSpatialPosition(spatial_position);

        	iter+=2;

    	}

    }

    }

    else
    {
        int iter=0;
    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{
    		int atom_id=(*atom)->GetID();

    		if( atom_id > 2 && atom_id < 10 )
    		{
        		Point <2> spatial_position;
        		spatial_position.SetXCoord(material_position[iter]);
        		spatial_position.SetYCoord(material_position[iter+1]);
        		spatial_position.SetZCoord(material_position[iter+2]);

        		(*atom)->SetSpatialPosition(spatial_position);

        		iter+=3;

    		}


    	}


    }


    if (dimension==2){

    int iter=0;
    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
    {
    	//only for atoms that are labeled

    	int atom_id=(*atom)->GetID();
//    	&& atom_id!=603 && atom_id!=604 && atom_id!=605
//    	    		&& atom_id!=606 && atom_id!=607 && atom_id!=608 && atom_id!=609 && atom_id!=610 && atom_id!=611
//    				&& atom_id!=612 && atom_id!=613 && atom_id!=614 && atom_id!=615 && atom_id!=616 && atom_id!=617

    	if( atom_id > 2 && atom_id < 5 )

    	{

        	Point <2> resultant_force=force.ResultantForce(*atom);

        	grad[iter]=-1*resultant_force.GetXCoord();
        	grad[iter+1]=-1*resultant_force.GetYCoord();

        	iter+=2;

    	}

    }

    }

    else
    {
        int iter=0;
    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{
    		int atom_id=(*atom)->GetID();

    		if( atom_id > 2 && atom_id < 10 )
    		{

        		Point <2> resultant_force=force.ResultantForce(*atom);

        		grad[iter]=resultant_force.GetXCoord();
        		grad[iter+1]=resultant_force.GetYCoord();
        		grad[iter+2]=resultant_force.GetZCoord();

        		iter+=3;

    		}


    	}


    }


    double eng=energy.TotPotentialEnergy(Atoms);
//    cout << "energy: " << eng <<endl;
//
//    double sum=0;
//
//    for (int i=0; i<grad.size() ; ++i)
//    {
//    	sum+=(grad[i]*grad[i]);
//    }
//
//    cout << "norm of force is: " << sqrt(sum) << endl;

    return(eng);

}


int main(){

//update the material position inside load step loop

    vector < Atom <2>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <2>(3, 3, 0, 1.1, 2);

    cout << "number of atoms: "<< unrelax_atoms.size() <<endl;


    typedef typename vector < Atom <2>* >::iterator At;
    //send material_position to real_1d_array x

    Boundary <2> boundary_bottom;
    boundary_bottom.SetBoundaryRegion(0.0,0.0,2.21,0.1, 1);

    Boundary <2> boundary_top;
    boundary_top.SetBoundaryRegion(0.0, 1.9, 2.21, 2.0, 2);

//    Boundary <2> crack;
//    crack.SetBoundaryRegion(0.0, 17.1, 15.5, 17.2, 3);

//    Boundary <2> crack_top;
//    crack_top.SetBoundaryRegion(0.0, 17.2, 15.5, 26.7, 4);
//
//    Boundary <2> crack_bottom;
//    crack_bottom.SetBoundaryRegion(0.0, 9.5, 15.5, 17.1, 5);

    vector < Boundary <2> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);

//    boundaries.push_back(crack);

//    boundaries.push_back(crack_top);
//    boundaries.push_back(crack_bottom);

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,10);



    vector < Atom <2>* > fixed_atoms;
    vector < Atom <2>* > bottom_atoms;
    vector < Atom <2>* > top_atoms;
    vector < Atom <2>* > optimized_atoms;

    vector <double> material_position;


    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ((*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
    	{
            Point <2> mposition;

            mposition=(*atom)->GetMaterialPosition();

            material_position.push_back(mposition.GetXCoord());
            material_position.push_back(mposition.GetYCoord());
            optimized_atoms.push_back(*atom);
    	}

    	else if((*atom) -> GetAtomRegion()==1)
    	{
    		bottom_atoms.push_back(*atom);
    	}

    	else if ((*atom) -> GetAtomRegion()==2)
    	{

    		top_atoms.push_back(*atom);
    	}

    	else if ((*atom)-> GetAtomRegion()==3)
    	{

    		fixed_atoms.push_back(*atom);

    	}

    }

    double max_load=0.2;
    double load=0;

    string vec_to_string;


    while (load < max_load)
    {

        InputData data;
        data.Atoms = atoms;
        data.sigma=1.0;
        data.epsilon=1.0;

        nlopt::opt opt(nlopt::LD_LBFGS, material_position.size());
        opt.set_min_objective(optimize_atomistic,(void*)&data);

        opt.set_xtol_rel(1e-10);
//        opt.set_maxeval(35);

        double minf;
        opt.optimize(material_position, minf);

//        vector <double> material_position;




        string file_name=to_string(load);

//        const char *path="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/relaxed_config.xtl";

    	ofstream ofile(file_name+"relaxed_config_n.xtl");

    	if (ofile.is_open())
    	{

    		ofile <<"TITLE XCrySDen XSF file" <<endl;
    		ofile <<"CELL" << endl;

    		ofile <<" "<<"100.0"<<" " <<"100.0" <<" "<< "100.0" <<
    				" " << "90.0" << " " << "90.0" << " " << "90.0"<<endl;

    		ofile << "SYMMETRY NUMBER 1" <<endl;
    		ofile << "SYMMETRY LABEL  P1" << endl;
    		ofile << "ATOMS" <<endl;
    		ofile <<"NAME" <<"    " << "X" <<"     "<<"Y" << "    " << "Z" <<endl;
    		int iteration=0;

    		for (At at=atoms.begin(); at!=atoms.end(); ++at)
    		{
    			Point <2> atomP=(*at)->GetSpatialPosition();
    			double atomP_x=atomP.GetXCoord();
    			double atomP_y=atomP.GetYCoord();

    			cout << "x: " << atomP_x << "y: " <<  atomP_y << endl;

    			if( (*at)-> GetAtomRegion()==0 || (*at)->GetAtomRegion()==4 || (*at)->GetAtomRegion()==5 )
    			{
        			Point <2> spatial_position =(*at)->GetSpatialPosition();
        			double spatialp_x=spatial_position.GetXCoord();
        			double spatialp_y=spatial_position.GetYCoord();

        			ofile << "Mg"<< "      " <<setprecision(5)<< material_position[iteration]<<"    "
        					<< setprecision(5)<< material_position[iteration+1] << "    "<< "0.0"<<endl;
        			iteration+=2;

    			}
    			else{

        			Point <2> spatial_position =(*at)->GetSpatialPosition();
        			double spatialp_x=spatial_position.GetXCoord();
        			double spatialp_y=spatial_position.GetYCoord();

        			ofile << "Mg"<< "      " <<setprecision(5)<< spatialp_x<<"    "
        					<< setprecision(5)<< spatialp_y << "    "<< "0.0"<<endl;

    			}

    		}
    		ofile << "EOF" <<endl;

    	}
    	else
    	{
    	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
    	}

        load+=0.1;


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

        	}

        }

        for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
        {
        	if( (*atom)->GetAtomRegion()==1 )
        	{

                Point <2> Poldb = (*atom) -> GetMaterialPosition();
                double Yold = Poldb.GetYCoord();
                double Xold = Poldb.GetXCoord();

                Point <2> Pnewb;
                double Ynew= Yold - ( load/2) ;
                Pnewb.SetXCoord(Xold);
                Pnewb.SetYCoord(Ynew);

                (*atom) -> SetSpatialPosition(Pnewb);

        	}

        }


    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	return 0;

}






