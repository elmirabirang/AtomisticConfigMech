/*
 * test_SW.cpp
 *
 *  Created on: Aug 3, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#include "energy.h"
#include "energy.cpp"

#include "force.h"
#include "force.cpp"

#include "point.h"
//#include "point.cpp"

#include "atom.h"
//#include "atom.cpp"

#include "bond.h"
#include "bond.cpp"

#include "matrix.h"

#include "third_order_tensor.h"
#include "third_order_tensor.cpp"
#include "write_data_file.cpp"
#include "unrelaxed_config.cpp"

#include <math.h>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <list>
#include <algorithm>
#include <iomanip> 

int main()
{

//	Point <3> atom_alpha_matP(0.5,0.5,0.);
//	Point <3> atom_beta_matP(0.0,0.0,0.0);
//	Point <3> atom_gamma_matP(0.,1.,0.);
//
//	Point <3> atom_alpha_spaP(0.5,0.5,0.);
//	Point <3> atom_beta_spaP(0.0,0.0,0.0);
//	Point <3> atom_gamma_spaP(0.,1.,0.);
//
//	Atom <3> atom_alpha;
//	Atom <3> atom_beta;
//	Atom <3> atom_gamma;
//
//	atom_alpha.SetID(0);
//	atom_beta.SetID(1);
//	atom_gamma.SetID(2);
//
//	atom_alpha.SetMaterialPosition(atom_alpha_matP);
//	atom_beta.SetMaterialPosition(atom_beta_matP);
//	atom_gamma.SetMaterialPosition(atom_gamma_matP);
//
//	atom_alpha.SetSpatialPosition(atom_alpha_spaP);
//	atom_beta.SetSpatialPosition(atom_beta_spaP);
//	atom_gamma.SetSpatialPosition(atom_gamma_spaP);
//
//	Energy <3> energy;
//	
//	Force <3> force;
//
//	double deformational_eng=energy.deformationalSWtwoBody(atom_alpha,atom_beta,2.0951, 2.1683,
//								                           7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//
//	double normal_eng=energy.stillingerWeberTwoBody(atom_alpha,atom_beta,2.0951, 2.1683,
//								                    7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//
//	double deformational_3BodyEng=energy.deformationalSWthreeBody(atom_alpha,atom_beta,atom_gamma,
//			                                                      2.0951, 2.0951,2.1683, 1.2, 21.,
//																  -0.333,1.8, 1.8 );
//
//	double normal_3BodyEng=energy.stillingerWeberThreeBody((atom_beta),(atom_alpha), (atom_gamma),
//			                                                2.0951, 2.0951,2.1683, 1.2, 21.,
//									                        -0.333,1.8, 1.8 );
//	
//	double configurational_eng=energy.configurationalSWtwoBody(atom_alpha,atom_beta,2.0951, 2.1683,
//                                                               7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//	
//	double configurational_3BodyEng=energy.configurationalSWthreeBody(atom_alpha,atom_beta,atom_gamma,
//                                                                      2.0951, 2.0951,2.1683, 1.2, 21.,
//			                                                          -0.333,1.8, 1.8);
//
//	
//	Point <3> DefBondSWForce=force.deformationalBondSWForce(atom_alpha, atom_beta, 2.0951, 2.1683,
//									                        7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//	
//	Point <3> ConfigBondSWForce=force.configurationalBondSWForce(atom_alpha, atom_beta, 2.0951, 2.1683,
//									                             7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//
//
//	double DefBondSWForce_normal=DefBondSWForce.PointNorm();
//
//	Point <3> BondSWForce=force.StillingerWeberTwoBodyForce(atom_alpha, atom_beta, 2.0951, 2.1683,
//								                            7.049556277, 0.6022245584, 4.0, 0.0, 1.8);
//
//	double BondSWForce_normal=BondSWForce.PointNorm();
//
//	Point <3> forceH2H3=force.StillingerWeberThreeBodyForceH2H3(atom_alpha, atom_beta, atom_gamma, 2.0951, 2.0951,2.1683, 1.2, 21.,
//                                                              -0.333333333333,1.8, 1.8 );
//	Point <3> forceH1=force.StillingerWeberThreeBodyForceH1( atom_beta, atom_alpha,atom_gamma, 2.0951, 2.0951,2.1683, 1.2, 21.,
//                                                              -0.333333333333,1.8, 1.8 );
//	double norm_forceH1=forceH1.PointNorm();
//	double norm_forceH2H3=forceH2H3.PointNorm();
//
////	std::cout << "norm_forceH1: " << norm_forceH1 << "  norm_forceH2H3: " << norm_forceH2H3 << endl;
//	
//	std:: cout << "(" << forceH1.GetXCoord() << ", "<< forceH1.GetYCoord() << ", "<< forceH1.GetZCoord() << ") "<< "\n";
//
//	Point <3> forcecp=force.deformationalcpSW3BodyForce(atom_alpha, atom_beta, atom_gamma,2.0951, 2.0951,2.1683, 1.2, 21.,
//                                                        -0.333,1.8, 1.8 );
//	
//	std:: cout << "(" <<  forcecp.GetXCoord() << ", "<<  forcecp.GetYCoord() << ", "<<  forcecp.GetZCoord() << ") "<< "\n";
//
//	Point <3> forcepr=force.deformationalprSW3BodyForce(atom_alpha, atom_beta, atom_gamma,2.0951, 2.0951,2.1683, 1.2, 21.,
//                                                        -0.333,1.8, 1.8);
//	
//	std:: cout << "(" <<  forcepr.GetXCoord() << ", "<<  forcepr.GetYCoord() 
//			   << ", "<<  forcepr.GetZCoord() << ") "<< "\n";
//
//	double force_CPNorm=forcecp.PointNorm();
//	double forcepr_norm=forcepr.PointNorm();
//	
//	Point <3> configForceCP=force.configurationalcpSW3BodyForce(atom_alpha, atom_beta, atom_gamma,2.0951, 2.0951,2.1683, 1.2, 21.,-0.333,1.8, 1.8);
//	Point <3> configForcePR=force.configurationalprSW3BodyForce(atom_alpha, atom_beta, atom_gamma,2.0951, 2.0951,2.1683, 1.2, 21.,-0.333,1.8, 1.8);
//
////	std::cout << "force_CPNorm: " << force_CPNorm << endl;
////	std::cout << "forcepr_norm" << forcepr_norm << endl;
//
//
//	std::cout << "Def Eng 2Body: "<< deformational_eng << " Eng 2Body: "<< normal_eng 
//			  << " Def Eng 3Body: "<< deformational_3BodyEng << " Eng 3Body: "
//			  << normal_3BodyEng  << " Config Eng 2Body: " << configurational_eng << 
//			  " Config Eng 3Body: " << configurational_3BodyEng << endl;
//	
//	std::cout << "Def x: "      << DefBondSWForce.GetXCoord()
//			  << " Def y: "     << DefBondSWForce.GetYCoord()
//			  << " Def z: "     << DefBondSWForce.GetZCoord()
//			  
//			  << " Config x: "  << ConfigBondSWForce.GetXCoord()
//			  << " Config y: "  << ConfigBondSWForce.GetYCoord()
//			  << " Config z: "  << ConfigBondSWForce.GetZCoord()
//			  
//			  << " Config CP x: " << configForceCP.GetXCoord()
//			  << " Config CP y: " << configForceCP.GetYCoord()
//			  << " Config CP z: " << configForceCP.GetZCoord()
//			  
//			  << " Config PR x: " << configForcePR.GetXCoord()
//			  << " Config PR y: " << configForcePR.GetYCoord()
//			  << " Config PR z: " << configForcePR.GetZCoord()
//			  
//			  << "\n";
	
	typedef typename vector < Atom <3>* >::reverse_iterator A;
	typedef typename vector < Atom <3>* >::iterator At;
//	 std::reverse_iterator<std::string::iterator> r = s.rbegin();
	
	Energy <3> energy;	
	Force <3> force;
	
	ofstream SFile;

	string SPath="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/results/dim3/silicon_test/dump.spatial";

	stringstream ss;
	ss << 0;
	
	string out_string;
	out_string=ss.str();

	string Sfile_name=SPath+out_string;
	
	SFile.open(Sfile_name.c_str());
	
    using namespace std;
    vector < Atom <3>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <3> (2, 2, 2, 5.431, 3, "Si");
    
//    	Point <3> atom1_matP(0.,0.,0.);
//    	Point <3> atom2_matP(2.7155, 2.7155, 0.);
//    	Point <3> atom3_matP(1.35775, 1.35775, 1.6);
//    	Point <3> atom4_matP(0., 2.7155, 2.7155);
//    	Point <3> atom5_matP(2.7155, 0., 2.7155);
//    
//    	Point <3> atom1_spaP(0.,0.,0.);
//    	Point <3> atom2_spaP(2.7155, 2.7155, 0.);
//    	Point <3> atom3_spaP(1.35775, 1.35775, 1.6);
//    	Point <3> atom4_spaP(0., 2.7155, 2.7155);
//    	Point <3> atom5_spaP(2.7155, 0., 2.7155);
//    
//    	Atom <3> atom1;
//    	Atom <3> atom2;
//    	Atom <3> atom3;
//    	Atom <3> atom4;
//    	Atom <3> atom5;
//    
//    	atom1.SetID(0);
//    	atom2.SetID(1);
//    	atom3.SetID(2);
//    	atom4.SetID(3);
//    	atom5.SetID(4);
//    
//    	atom1.SetMaterialPosition(atom1_matP);
//    	atom2.SetMaterialPosition(atom2_matP);
//    	atom3.SetMaterialPosition(atom3_matP);
//    	atom4.SetMaterialPosition(atom4_matP);
//    	atom5.SetMaterialPosition(atom5_matP);
//    
//    	atom1.SetSpatialPosition(atom1_spaP);
//    	atom2.SetSpatialPosition(atom2_spaP);
//    	atom3.SetSpatialPosition(atom3_spaP);
//    	atom4.SetSpatialPosition(atom4_spaP);
//    	atom5.SetSpatialPosition(atom5_spaP);
//    	
//    	vector <Atom <3>*> unrelax_atoms;
//    	
//    	unrelax_atoms.push_back(&atom1);
//    	unrelax_atoms.push_back(&atom2);
//    	unrelax_atoms.push_back(&atom3);
//    	unrelax_atoms.push_back(&atom4);
//    	unrelax_atoms.push_back(&atom5);
    	
    for (At atom=unrelax_atoms.begin(); atom!=unrelax_atoms.end(); ++atom)
    {
    	
    	(*atom)->SetAtomRegion(0);
    	
    	if ((*atom)->GetID()==14)
    	{
    		Point <3> spatial_position(1.35775,1.35775,1.6);
    		(*atom)->SetSpatialPosition(spatial_position);
    	}

    }
    vector < Atom <3>* > atoms=FindNeighbors(unrelax_atoms,2.6);
    
    string path4="/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/dump.silicon" ;

    writeDataFile (atoms, 1, path4);
    
    double total_energy=0.;
    
	if ( SFile.is_open() )
	{
    
	SFile <<"ITEM: TIMESTEP" <<endl;
	SFile << 0 <<endl;
	SFile <<"ITEM: NUMBER OF ATOMS" << endl;
	SFile << atoms.size() << endl;
	SFile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;
	SFile << "0" << " " << "1.6" << endl;
	SFile << "0" << " " << "1.6" << endl;
	SFile << "0" << " " << "1.6" << endl;
	SFile <<"ITEM: ATOMS id type x y z fx fy fz" <<endl;
	
	for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
	{
		
		force.ResultantSWThreeBodyForce((*atom),2.0951, 2.0951,
								   2.1683, 1.2, 21.,-0.333333333333,
								   1.8, 1.8);
		
	}
	

    for (At atom=atoms.begin(); atom!=atoms.end(); ++atom)
    {
    	
    	double energy_atom_a=0.;
    	
		energy_atom_a=energy.totalStillingerWeberEnergy((*atom),
														2.0951, 2.0951,
														2.1683, 1.2, 21.,-0.33333333333,
														1.8, 1.8,
														7.049556277, 0.6022245584,
														4.0, 0.0);
		
		Point <3> spatial_position=(*atom)->GetSpatialPosition();
		
		double spatialp_x=spatial_position.GetXCoord();
		double spatialp_y=spatial_position.GetYCoord();
		double spatialp_z=spatial_position.GetZCoord();
		
		Point <3> atomic_force(0.,0.,0.);

		atomic_force=force.ResultantSWTwoBodyForce((*atom), 2.0951, 2.1683,
															 7.049556277, 0.6022245584,
															 4.0, 0.0, 1.8);
		double atomic_force_x=atomic_force.GetXCoord();
		double atomic_force_y=atomic_force.GetYCoord();
		double atomic_force_z=atomic_force.GetZCoord();
		
		Point <3> atomic_force3body(0.,0.,0.);

//		atomic_force3body=force.ResultantSWThreeBodyForce((*atom),2.0951, 2.0951,
//																  2.1683, 1.2, 21.,-0.333333333333,
//																1.8, 1.8);
		
		atomic_force3body=(*atom)-> GetForce();
		
		
		double atomic_force3body_x=atomic_force3body.GetXCoord();
		double atomic_force3body_y=atomic_force3body.GetYCoord();
		double atomic_force3body_z=atomic_force3body.GetZCoord();
		
		SFile << (*atom) -> GetID() << " " << "0" 
			 <<" "  <<  setprecision(5) << spatialp_x
			 << " " <<  setprecision(5) << spatialp_y 
			 << " " <<  setprecision(5) << spatialp_z
			 << " " <<  atomic_force_x+atomic_force3body_x
			 << " " <<  atomic_force_y+atomic_force3body_y
			 << " " <<  atomic_force_z+atomic_force3body_z
			 <<endl;
		
		total_energy+=energy_atom_a;
		
    }
	}
    
    SFile.close();
    
    std::cout << "total energy:" << total_energy << "\n";
    

			  
	
//	std::cout << "def Bond SW Force: " << DefBondSWForce_normal << "  Bond SW Force: " << BondSWForce_normal << endl;

    return 0;
}




