///*
// * optimizer.cpp
// *
// *  Created on: 20th Feb, 2018
// *      Author: S.Elmira Birang.O
// *      Mechanical Engineering PhD Candidate
// *      Chair of Applied Mechanics
// *      University of Erlangen-Nürnberg
// */
//
//
//#include "atom.h"
//#include "point.h"
//#include "energy.h"
//#include "force.h"
//#include "boundary.h"
//#include "bond.h"
//
//#include "unrelaxed_config.cpp"
//#include "energy.cpp"
//#include "force.cpp"
//#include "boundary.cpp"
//#include "assign_region.cpp"
//#include "bond.cpp"
//
//#include <vector>
//#include <string>
//
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//
//#include <string>
//
//#include <lbfgs.h>
//
//
//
//using namespace std;
//
//struct InputData
//{
//    vector < Atom <2>* > Atoms;
//    double sigma;
//    double epsilon;
//
//};
//
//static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
//{
//
//	InputData data_tmp = *((InputData *)instance);
//
//
//    double sigma=data_tmp.sigma;
//    double epsilon=data_tmp.epsilon;
//
//    vector < Atom <2>* > Atoms;
//    Atoms=data_tmp.Atoms;
//
//    Energy <2> energy(sigma,epsilon);
//    Force <2> force(sigma,epsilon);
//
//    int dimension=2;
//
//    typedef typename vector < Atom <2>* >::iterator At;
//
//    switch (dimension){
//
//    case 2:{
//
//    int iter=0;
//    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
//    {
//    	//only for atoms that are labeled
//
//    	int atom_region=(*atom)->GetAtomRegion();
//
//
//    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)
//
//    	{
//
//        	Point <2> spatial_position;
//        	spatial_position.SetXCoord(x[iter]);
//        	spatial_position.SetYCoord(x[iter+1]);
//
//        	(*atom)->SetSpatialPosition(spatial_position);
//
//        	iter+=2;
//
//    	}
//
//    }
//    break;
//    }
//
//    case 3:{
//
//        int iter=0;
//    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
//    	{
//    		int atom_id=(*atom)->GetID();
//        	Atom <2> *Atm;
//        	Atm=*atom;
//
//    		if( atom_id > 2 && atom_id < 10 )
//    		{
//        		Point <2> spatial_position;
//        		spatial_position.SetXCoord(x[iter]);
//        		spatial_position.SetYCoord(x[iter+1]);
//        		spatial_position.SetZCoord(x[iter+2]);
//
//        		(*atom)->SetSpatialPosition(spatial_position);
//
//        		iter+=3;
//
//    		}
//
//
//    	}
//
//    	break;
//    }
//
//    }
//
//
//
//    switch (dimension)
//    {
//
//
//    case 2:{
//
//    int iter=0;
//    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
//    {
//    	//only for atoms that are labeled
//
//    	int atom_region=(*atom)->GetAtomRegion();
//
//    	if(atom_region!=1 && atom_region!=2 && atom_region!=3)
//
//    	{
//
////        	Point <2> resultant_force=force.ResultantForce(*atom);
//
//        	Point <2> resultant_force=force.ResultantConfigForce(*atom);
//
//        	lbfgsfloatval_t gx=-1*resultant_force.GetXCoord();
//        	lbfgsfloatval_t gy=-1*resultant_force.GetYCoord();
//
//        	g[iter]=gx;
//        	g[iter+1]=gy;
//
//        	iter+=2;
//
//    	}
//
//    }
//    break;
//    }
//
//
//    case 3:{
//
//        int iter=0;
//    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
//    	{
//    		int atom_id=(*atom)->GetID();
//        	Atom <2> *Atm;
//        	Atm=*atom;
//
//    		if( atom_id > 2 && atom_id < 10 )
//    		{
//
//        		Point <2> resultant_force=force.ResultantForce(*atom);
//
//        		lbfgsfloatval_t gx=resultant_force.GetXCoord();
//        		lbfgsfloatval_t gy=resultant_force.GetYCoord();
//        		lbfgsfloatval_t gz=resultant_force.GetZCoord();
//
//        		g[iter]=gx;
//        		g[iter+1]=gy;
//        		g[iter+2]=gz;
//
//        		iter+=3;
//
//    		}
//
//
//    	}
//    	break;
//    }
//
//    }
//
//
////    double eng=energy.TotPotentialEnergy(Atoms);
//
//    double eng=energy.ConfigTotPotentialEnergy(Atoms);
//
//    lbfgsfloatval_t energy2d=eng;
//
//
//    cout << "energy: " << eng <<endl;
//
//    return(energy2d);
//
//}
//
//static int progress (
//		void *instance,
//		const lbfgsfloatval_t *x,
//		const lbfgsfloatval_t *g,
//		const lbfgsfloatval_t fx,
//		const lbfgsfloatval_t xnorm,
//		const lbfgsfloatval_t gnorm,
//		const lbfgsfloatval_t step,
//		int n,
//		int k,
//		int ls
//		)
//{
//    printf("Iteration %d:\n", k);
//    printf("  fx = %f= ", fx);
//    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//    printf("\n");
//    return 0;
//
//}
//
//
//int main(){
//
////update the material position inside load step loop
//
//    vector < Atom <2>* > unrelax_atoms;
//    unrelax_atoms=UnrelaxedConfigGenerator <2>(50, 29, 0, 1.1, 2);
//
//    cout << "number of atoms: "<< unrelax_atoms.size() <<endl;
//
//
//    typedef typename vector < Atom <2>* >::iterator At;
//    //send material_position to real_1d_array x
//
//    Boundary <2> boundary_bottom;
//    boundary_bottom.SetBoundaryRegion(0.0, 0.0, 54.0, 0.1, 1);
//
//    Boundary <2> boundary_top;
//    boundary_top.SetBoundaryRegion(0.0, 26.6, 54.0, 26.8, 2);
//
//    Boundary <2> crack_top;
//    crack_top.SetBoundaryRegion(0.0, 12.0, 12, 16.2, 4);
//
//    Boundary <2> crack_bottom;
//    crack_bottom.SetBoundaryRegion(0.0, 6.0, 12, 11.9, 5);
//
////    Boundary <2> crack;
//
////    double xa=0;
////    double ya=11/2*1.9035;
////    double xb=6*1.1;
////    double yb=15/2*1.9035;
////    double xc=0;
////    double yc=19/2*1.9035;
//
////    double xa=0;
////    double ya=11.4;
////    double xb=10*1.1+0.01;
////    double yb=12.4;
////
////    crack.SetBoundaryRegion(xa,ya,xb,yb,3);
//
////    crack.SetTriangBoundaryRegion(xa, ya, xb, yb, xc, yc,11);
//
//    vector < Boundary <2> > boundaries;
//
//    boundaries.push_back(boundary_bottom);
//    boundaries.push_back(boundary_top);
//    boundaries.push_back(crack_top);
//    boundaries.push_back(crack_bottom);
//
////    boundaries.push_back(crack);
//
//
//
//
//    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
//    {
//        AssignRegion(*atom,boundaries);
//
//    }
//
//    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,2.5);
//
//    vector < Atom <2>* > fixed_atoms;
//    vector < Atom <2>* > bottom_atoms;
//    vector < Atom <2>* > top_atoms;
//    vector < Atom <2>* > optimized_atoms;
//
//    int i, ret=0;
//    int N=(1436 - 90)*2;
//
//    lbfgsfloatval_t *x=lbfgs_malloc(N);
//    lbfgsfloatval_t fx;
//    lbfgs_parameter_t param;
//
//    double count=0;
//
//    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//    {
//
//    	if ((*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
//    	{
//            Point <2> mposition;
//
//            mposition=(*atom)->GetMaterialPosition();
//
//            x[i]=mposition.GetXCoord();
//            x[i+1]=mposition.GetYCoord();
//
//            optimized_atoms.push_back(*atom);
//            i+=2;
//            cout << "here " << endl;
//
//    	}
//
//    	else if((*atom) -> GetAtomRegion()==1)
//    	{
//    		bottom_atoms.push_back(*atom);
//    		count+=1;
//    	}
//
//    	else if ((*atom) -> GetAtomRegion()==2)
//    	{
//
//    		top_atoms.push_back(*atom);
//    		count+=1;
//    	}
//
//    	else if ((*atom)-> GetAtomRegion()==3)
//    	{
//
//    		fixed_atoms.push_back(*atom);
//
//    	}
//
//    }
//
//    cout << "count: " << count << endl;
//
////	ofstream ofile("config.xtl");
////
////	if (ofile.is_open())
////	{
////
////		ofile <<"TITLE XCrySDen XSF file" <<endl;
////		ofile <<"CELL" << endl;
////
////		ofile <<" "<<"100.0"<<" " <<"100.0" <<" "<< "100.0" <<
////				" " << "90.0" << " " << "90.0" << " " << "90.0"<<endl;
////
////		ofile << "SYMMETRY NUMBER 1" <<endl;
////		ofile << "SYMMETRY LABEL  P1" << endl;
////		ofile << "ATOMS" <<endl;
////		ofile <<"NAME" <<"    " << "X" <<"     "<<"Y" << "    " << "Z" <<endl;
////
////		for (At at=atoms.begin(); at!=atoms.end(); ++at)
////		{
////			int atom_region=(*at)-> GetAtomRegion();
////
////
////				Point <2> spatial_position =(*at)->GetSpatialPosition();
////				double spatialp_x=spatial_position.GetXCoord();
////				double spatialp_y=spatial_position.GetYCoord();
////				if ( atom_region!=3)
////				{
////
////				ofile << "Mg"<< "      " <<setprecision(5)<< spatialp_x<<"    "
////						<< setprecision(5)<< spatialp_y << "    "<< "0.0"<<endl;
////
////			}
////
////			else
////			{
////
////				ofile << "Na"<< "      " <<setprecision(5)<< spatialp_x<<"    "
////						<< setprecision(5)<< spatialp_y << "    "<< "0.0"<<endl;
////
////			}
////
////
////
////		}
////		ofile << "EOF" <<endl;
////
////	}
////	else
////	{
////	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
////	}
//
//
//    double max_load=0.1;
//    double load=0;
//
//    string vec_to_string;
//
//
//    while (load < max_load)
//    {
//
//        InputData data;
//        data.Atoms = atoms;
//        data.sigma=1;
//        data.epsilon=1;
//        int id=0;
//        Force <2> force(data.sigma,data.epsilon);
//
//        lbfgs_parameter_init(&param);
//
//        ret=lbfgs(N,x,&fx,evaluate, progress, (void*)&data, &param);
//
//        string file_name=to_string(load);
//
////    	ofstream ofile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/dump.crack"+file_name);
////
////    	if (ofile.is_open())
////    	{
////
////    		ofile <<"ITEM: TIMESTEP" <<endl;
////    		ofile <<"0" <<endl;
////    		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
////    		ofile << atoms.size() << endl;
////    		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
////
////    		ofile << "SYMMETRY NUMBER 1" <<endl;
////    		ofile << "SYMMETRY LABEL  P1" << endl;
////    		ofile << "ATOMS" <<endl;
////    		ofile <<"ITEM: ATOMS id type xs ys zs fx fy fz" <<endl;
////
////    		for (At at=atoms.begin(); at!=atoms.end(); ++at)
////    		{
////
////    			if((*at)->GetAtomRegion()!=3)
////
////    			{
////    				Point <2> config_force =force.ConfigurationalForce(*at);
////
////        			Point <2> spatial_position =(*at)->GetSpatialPosition();
////
////        			double spatialp_x=spatial_position.GetXCoord();
////        			double spatialp_y=spatial_position.GetYCoord();
////
////        			double config_force_x=config_force.GetXCoord();
////        			double config_force_y=config_force.GetYCoord();
////
////        			ofile << id<< " " << "4" << spatialp_x<<" "
////        					<< setprecision(5)<< spatialp_y << " "<< "0.0"
////							" "<< config_force_x << " " << config_force_y << " " << "0.0" <<endl;
////        			id+=1;
////
////    			}
////
////    		}
////    		ofile << "EOF" <<endl;
////
////    	}
////    	else
////    	{
////    	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
////    	}
//
//
//        load+=0.1;
//
//        for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
//        {
//        	if( (*atom)->GetAtomRegion()==2 )
//        	{
//            	Point <2> Poldt = (*atom) -> GetMaterialPosition();
//            	double Yold = Poldt.GetYCoord();
//            	double Xold = Poldt.GetXCoord();
//
//            	Point <2> Pnewt;
//            	double Ynew= Yold + ( load/2);
//            	Pnewt.SetXCoord(Xold);
//            	Pnewt.SetYCoord(Ynew);
//
////                (*atom) -> SetMaterialPosition(Pnewt);
//                (*atom) -> SetSpatialPosition(Pnewt);
//
//        	}
//
//        	else if( (*atom)->GetAtomRegion()==1 )
//        	{
//
//                Point <2> Poldb = (*atom) -> GetMaterialPosition();
//                double Yold = Poldb.GetYCoord();
//                double Xold = Poldb.GetXCoord();
//
//                Point <2> Pnewb;
//                double Ynew= Yold - (load/2) ;
//                Pnewb.SetXCoord(Xold);
//                Pnewb.SetYCoord(Ynew);
//
////                (*atom) -> SetMaterialPosition(Pnewb);
//                (*atom) -> SetSpatialPosition(Pnewb);
//
//        	}
//
//        	else if ( (*atom) -> GetAtomRegion()==3)
//        	{
//        		Point <2> Poldf= (*atom) -> GetSpatialPosition();
//
//        		(*atom) -> SetSpatialPosition(Poldf);
//
//        	}
//
//
//        }
//
//    }
//
//    return 0;
//
//
//    }
//
//
//
//
//
//
//
//
//
/*
 * optimizer.cpp
 *
 *  Created on: 20th Feb, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nürnberg
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>

#include <lbfgs.h>



using namespace std;


struct InputData
{
    vector < Atom <2>* > Atoms;
    double sigma;
    double epsilon;

};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{

	InputData data_tmp = *((InputData *)instance);


    double sigma=data_tmp.sigma;
    double epsilon=data_tmp.epsilon;


    vector < Atom <2>* > Atoms;
    Atoms=data_tmp.Atoms;

    Energy <2> energy(sigma,epsilon);
    Force <2> force(sigma,epsilon);

    int dimension=2;

    typedef typename vector < Atom <2>* >::iterator At;

    switch (dimension){

    case 2:{

    int iter=0;
    for (At atom=Atoms.begin(); atom!=Atoms.end(); ++atom )
    {
    	//only for atoms that are labeled

    	int atom_region=(*atom)->GetAtomRegion();


    	if( atom_region!=1 && atom_region!=2 && atom_region!=3)

    	{

        	Point <2> spatial_position;
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
    		int atom_id=(*atom)->GetID();
        	Atom <2> *Atm;
        	Atm=*atom;

    		if( atom_id > 2 && atom_id < 10 )
    		{
        		Point <2> spatial_position;
        		spatial_position.SetXCoord(x[iter]);
        		spatial_position.SetYCoord(x[iter+1]);
        		spatial_position.SetZCoord(x[iter+2]);

        		(*atom)->SetSpatialPosition(spatial_position);

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

        	Point <2> resultant_force=force.ResultantForce(*atom);

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

        int iter=0;
    	for(At atom=Atoms.begin(); atom!=Atoms.end(); ++atom)
    	{
    		int atom_id=(*atom)->GetID();
        	Atom <2> *Atm;
        	Atm=*atom;

    		if( atom_id > 2 && atom_id < 10 )
    		{

        		Point <2> resultant_force=force.ResultantForce(*atom);

        		lbfgsfloatval_t gx=resultant_force.GetXCoord();
        		lbfgsfloatval_t gy=resultant_force.GetYCoord();
        		lbfgsfloatval_t gz=resultant_force.GetZCoord();

        		g[iter]=gx;
        		g[iter+1]=gy;
        		g[iter+2]=gz;

        		iter+=3;

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
    printf("Iteration %d:\n", k);
    printf("  fx = %f= ", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;

}


int main(){

//update the material position inside load step loop

    vector < Atom <2>* > unrelax_atoms;
    //unrelax_atoms=UnrelaxedConfigGenerator <2>(50, 29, 0, 1.1,2.5);
    unrelax_atoms=UnrelaxedConfigGenerator <2>(50, 29, 0, 1.1, 2);

    cout << "number of atoms: "<< unrelax_atoms.size() <<endl;


    typedef typename vector < Atom <2>* >::iterator At;
    //send material_position to real_1d_array x

    Boundary <2> boundary_bottom;
//    boundary_bottom.SetBoundaryRegion(0.0, 0.0, 54.0, 0.1, 1);
      boundary_bottom.SetBoundaryRegion(0.0, 0.0, 54.0, 0.1, 1);

    Boundary <2> boundary_top;
//    boundary_top.SetBoundaryRegion(0.0, 26.6, 54.0, 26.8, 2);
    boundary_top.SetBoundaryRegion(0.0, 26.6, 54.0, 26.8, 2);


    Boundary <2> crack_top;
//    crack_top.SetBoundaryRegion(0.0, 12.0, 12, 16.2, 4);
    crack_top.SetBoundaryRegion(0.0, 12.0, 12, 16.2, 4);

    Boundary <2> crack_bottom;
//    crack_bottom.SetBoundaryRegion(0.0, 6.0, 12, 11.9, 5);
    crack_bottom.SetBoundaryRegion(0.0, 6.0, 12, 11.9, 5);

    Boundary <2> crack;

//    double xa=0;
//    double ya=11/2*1.9035;
//    double xb=6*1.1;
//    double yb=15/2*1.9035;
//    double xc=0;
//    double yc=19/2*1.9035;

    double xa=0;
    double ya=6*1.9035+0.1;
    double xb=10*1.1+0.01;
    double yb=7.5*1.9035+0.1;

//    crack.SetBoundaryRegion(xa,ya,xb,yb,3);

//    crack.SetTriangBoundaryRegion(xa, ya, xb, yb, xc, yc,11);

    vector < Boundary <2> > boundaries;

    boundaries.push_back(boundary_bottom);
    boundaries.push_back(boundary_top);
    boundaries.push_back(crack);


    boundaries.push_back(crack_top);
    boundaries.push_back(crack_bottom);

    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
        AssignRegion(*atom,boundaries);

    }

    vector < Atom <2>* > atoms=FindNeighbors(unrelax_atoms,2.5);

    vector < Atom <2>* > fixed_atoms;
    vector < Atom <2>* > bottom_atoms;
    vector < Atom <2>* > top_atoms;
    vector < Atom <2>* > optimized_atoms;

    int i, ret=0;
    int N=(1436 - 90)*2;
    //int N=(1436 - 100)*2;

    lbfgsfloatval_t *x=lbfgs_malloc(N);
    lbfgsfloatval_t fx;
    lbfgs_parameter_t param;

    double count=0;

    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {

    	if ((*atom)->GetAtomRegion()==0 || (*atom)->GetAtomRegion()==4 || (*atom)->GetAtomRegion()==5)
    	{
            Point <2> mposition;

            mposition=(*atom)->GetMaterialPosition();

            x[i]=mposition.GetXCoord();
            x[i+1]=mposition.GetYCoord();

            optimized_atoms.push_back(*atom);
            i+=2;
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
    		count+=1;

    	}

    }

    cout << "count: " << count << endl;

//	ofstream ofile("config.xtl");
//
//	if (ofile.is_open())
//	{
//
//		ofile <<"TITLE XCrySDen XSF file" <<endl;
//		ofile <<"CELL" << endl;
//
//		ofile <<" "<<"100.0"<<" " <<"100.0" <<" "<< "100.0" <<
//				" " << "90.0" << " " << "90.0" << " " << "90.0"<<endl;
//
//		ofile << "SYMMETRY NUMBER 1" <<endl;
//		ofile << "SYMMETRY LABEL  P1" << endl;
//		ofile << "ATOMS" <<endl;
//		ofile <<"NAME" <<"    " << "X" <<"     "<<"Y" << "    " << "Z" <<endl;
//
//		for (At at=atoms.begin(); at!=atoms.end(); ++at)
//		{
//			int atom_region=(*at)-> GetAtomRegion();
//
//
//				Point <2> spatial_position =(*at)->GetSpatialPosition();
//				double spatialp_x=spatial_position.GetXCoord();
//				double spatialp_y=spatial_position.GetYCoord();
//				if ( atom_region!=3)
//				{
//
//				ofile << "Mg"<< "      " <<setprecision(5)<< spatialp_x<<"    "
//						<< setprecision(5)<< spatialp_y << "    "<< "0.0"<<endl;
//
//			}
//
//			else
//			{
//
//				ofile << "Na"<< "      " <<setprecision(5)<< spatialp_x<<"    "
//						<< setprecision(5)<< spatialp_y << "    "<< "0.0"<<endl;
//
//			}
//
//
//
//		}
//		ofile << "EOF" <<endl;
//
//	}
//	else
//	{
//	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
//	}


    double max_load=1.0;
    double load=0;
    int load_step=0;

    string vec_to_string;

    while (load < max_load)
    {

        InputData data;
        data.Atoms = atoms;
        data.sigma=1.0;
        data.epsilon=1.0;


        Force <2> force(data.sigma,data.epsilon);

        lbfgs_parameter_init(&param);
        param.epsilon=1e-7;

        ret=lbfgs(N,x,&fx,evaluate, progress, (void*)&data, &param);


        int id=0;
        load_step+=1;
        string file_name=std::to_string(load);

        if (load==0){

        	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
        	{

        		Point <2> spatialP= (*atom)->GetSpatialPosition();
        		(*atom)->SetMaterialPosition(spatialP);

        	}

        }


    	ofstream file("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/spatial/dump.crackspatial"+file_name);
    	ofstream ofile("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim2/material/dump.crack"+file_name);

        if (file.is_open())
            	{

            		ofile <<"ITEM: TIMESTEP" <<endl;
            		ofile << load_step <<endl;
            		ofile <<"ITEM: NUMBER OF ATOMS" << endl;
            		ofile << atoms.size() << endl;
            		ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
            		ofile <<"0" <<" " << "50" <<endl;
            		ofile <<"0" <<" " << "50" <<endl;
            		ofile <<"0" <<" " << "50" <<endl;
            		ofile <<"ITEM: ATOMS id type xs ys zs fx fy fz" <<endl;

            		file <<"ITEM: TIMESTEP" <<endl;
            		file << load_step <<endl;
            		file <<"ITEM: NUMBER OF ATOMS" << endl;
            		file << atoms.size() << endl;
            		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
            		file <<"0" <<" " << "50" <<endl;
            		file <<"0" <<" " << "50" <<endl;
            		file <<"0" <<" " << "50" <<endl;
            		file <<"ITEM: ATOMS id type xs ys zs fx fy fz" <<endl;

            		for (At at=atoms.begin(); at!=atoms.end(); ++at)
            		{

            			if((*at)->GetAtomRegion()!=3)

            			{
            				Point <2> config_force =force.ConfigurationalForce(*at);
                			Point <2> matposition = (*at)->GetMaterialPosition();
                			Point <2> spatial_position = (*at)-> GetSpatialPosition();

                			double materialp_x=matposition.GetXCoord();
                			double materialp_y=matposition.GetYCoord();

                			double spatialp_x=spatial_position.GetXCoord();
                			double spatialp_y=spatial_position.GetYCoord();

                			double config_force_x=config_force.GetXCoord();
                			double config_force_y=config_force.GetYCoord();

                			ofile << id<< " " << "4" <<" " << materialp_x<<" "
                				  << setprecision(5)<< materialp_y << " "<< "0.0"
        						  " "<< config_force_x << " " << config_force_y << " " << "0.0" <<endl;

                			file << id<< " " << "4" <<" " << spatialp_x<<" "
                				  << setprecision(5)<< spatialp_y << " "<< "0.0"<<
								  " "<< config_force_x << " " << config_force_y << " " << "0.0" <<endl;

                			id+=1;

            			}

            		}

            	}

          else
            	{
            	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
            	}

        load+=0.05;

        Point <2> top_force;
        top_force.SetXCoord(0.0);
        top_force.SetYCoord(0.0);

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

//                (*atom) -> SetMaterialPosition(Pnewt);
                (*atom) -> SetSpatialPosition(Pnewt);
                Point <2> resultant_force=force.ResultantForce(*atom);
                top_force=top_force+resultant_force;
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

//                (*atom) -> SetMaterialPosition(Pnewb);
                (*atom) -> SetSpatialPosition(Pnewb);

        	}

        	else if ( (*atom) -> GetAtomRegion()==3)
        	{
        		Point <2> Poldf= (*atom) -> GetSpatialPosition();

        		(*atom) -> SetSpatialPosition(Poldf);

        	}

        }

        double top_forces_norm=top_force.PointNorm();
        ofile << top_forces_norm << endl;

    }

    return 0;


    }










