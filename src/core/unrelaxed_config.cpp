/*
 * unrelaxed_config.cpp
 * This file generates an unrelaxed configuration for  crystalline structures.
 * generate the points
 * set the points and their id to construct atoms.
 *
 *  Created on: Nov 29, 2018
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */


#include <vector>
#include <string>
#include "point.h"
#include "atom.h"
#include "bond.h"
#include "boundary.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

template <int dim>
vector < Atom <dim>* > UnrelaxedConfigGenerator(int number_x, int number_y, int number_z, double lattice_const, int dimension, string type)
{

//	number_x: number of atoms in x direction
//  number_y: number of atoms in y direction.
//	number_z: number of atoms in z direction.

	vector < Point <dim>* > unrelax_points;
	vector < Atom <dim>* > unrelax_atoms;


	if (dimension==2) {

	double x_lattice=lattice_const;
	double y_lattice=sqrt(3)*lattice_const;

	double x_max=x_lattice*number_x;
	double y_max=y_lattice*number_y;

	double y_e=0.00;
	double y_o=y_lattice/2;

	int y_iter=0;

	int ID=0;

	while(y_iter<number_y)
	{
		double x_e=0.00;
		double x_o=x_lattice/2;

		int x_iter=0;

		if (y_iter%2==0)
		{
			while (x_iter< number_x)
			{

	            Point <dim> *P=new Point <dim>;
	            P->SetXCoord(x_e);
	            P->SetYCoord(y_e);

	            Atom <dim> *atom=new Atom <dim>;
	            atom->SetInitialPosition(*P);
	            atom->SetMaterialPosition(*P);
	            atom->SetSpatialPosition(*P);
	            atom->SetID(ID);
	            atom->SetCellID(-1);

	            unrelax_points.push_back(P);
	            unrelax_atoms.push_back(atom);
	            ID+=1;

	            x_e+=x_lattice;
	            x_iter+=1;

			}
			y_e+=y_lattice;
			y_iter+=1;
		}

		else{
				while (x_iter< number_x && (x_o+x_lattice)<x_max)
				{

					Point <dim> *P=new Point <dim>;
					P->SetXCoord(x_o);
					P->SetYCoord(y_o);

					Atom <dim> *atom=new Atom <dim>;
					atom->SetInitialPosition(*P);
					atom->SetMaterialPosition(*P);
					atom->SetSpatialPosition(*P);
					atom->SetID(ID);
					atom->SetCellID(-1);

					unrelax_points.push_back(P);
					unrelax_atoms.push_back(atom);
					ID+=1;

					x_o+=x_lattice;
					x_iter+=1;


				}
				y_o+=y_lattice;
				y_iter+=1;

			}

		}


// uncomment it to generate 2d position file

	ofstream ofile("dump.unrelaxed_2dconfig");
	int id=0;

	if (ofile.is_open())
	{

	    ofile <<"ITEM: TIMESTEP" <<endl;
	    ofile <<"0" <<endl;
	    ofile <<"ITEM: NUMBER OF ATOMS" << endl;
	    ofile << unrelax_points.size() << endl;
	    ofile <<"ITEM: BOX BOUNDS mm mm pp"<<endl;

	    ofile <<"50" <<" " <<"50" << endl;
	    ofile <<"50" <<" " <<"50" << endl;
	    ofile <<"50" <<" " <<"50" << endl;

		ofile <<"ITEM: ATOMS id type x y z" <<endl;

		for (int i=0; i< unrelax_atoms.size(); ++i)

		    {

		       double materialp_x=unrelax_points[i]->GetXCoord();
		       double materialp_y=unrelax_points[i]->GetYCoord();


		       ofile << id<< " " << "0" <<" " << materialp_x<<" "
		            << setprecision(5)<< materialp_y << " "<< "0" <<endl;

		       id+=1;

		    }

	}
	else
	{
	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
	}

	}



	else if (dimension==3 && type!="Si"){


		int z_iter_3d=0;

		double x_lattice_3d=1.*lattice_const;
		double y_lattice_3d=1.*lattice_const;
		double z_lattice_3d=1.*lattice_const;

		double x_max_3d=x_lattice_3d*number_x;
		double y_max_3d=y_lattice_3d*number_y;
		double z_max_3d=z_lattice_3d*number_z ;

		double z_e_3d=0.0;
		double z_o_3d=z_lattice_3d/2;

		int ID_3d=0;

		while ( z_iter_3d < number_z)
		{

		if (z_iter_3d%2==0)
		{

		int y_iter_3d=0;

		double y_e_3d=0.00;
		double y_o_3d=y_lattice_3d/2;

		while(y_iter_3d <number_y)
		{
			double x_e_3d=0.00;
			double x_o_3d=x_lattice_3d/2;

			int x_iter_3d=0;

			if ( y_iter_3d%2==0 )
			{
				while (x_iter_3d < number_x)
				{

		            Point <dim> *P=new Point <dim>;
		            P->SetXCoord(x_e_3d);
		            P->SetYCoord(y_e_3d);
		            P-> SetZCoord(z_e_3d);

		            Atom <dim> *atom=new Atom <dim>;
		            atom->SetInitialPosition(*P);
		            atom->SetMaterialPosition(*P);
		            atom->SetSpatialPosition(*P);
		            atom->SetID(ID_3d);
		            atom->SetCellID(-1);

		            unrelax_points.push_back(P);
		            unrelax_atoms.push_back(atom);
		            ID_3d+=1;

		            x_e_3d+=x_lattice_3d;
		            x_iter_3d+=1;

				}
				y_e_3d+=y_lattice_3d;
				y_iter_3d+=1;
			}

			else{
					while (x_iter_3d< number_x && (x_o_3d+x_lattice_3d)<x_max_3d)
					{

						Point <dim> *P=new Point <dim>;
						P->SetXCoord(x_o_3d);
						P->SetYCoord(y_o_3d);
						P -> SetZCoord(z_e_3d);

						Atom <dim> *atom=new Atom <dim>;
						atom->SetInitialPosition(*P);
						atom->SetMaterialPosition(*P);
						atom->SetSpatialPosition(*P);
						atom->SetID(ID_3d);
						atom->SetCellID(-1);

						unrelax_points.push_back(P);
						unrelax_atoms.push_back(atom);
						ID_3d+=1;

						x_o_3d+=x_lattice_3d;
						x_iter_3d+=1;


					}
					y_o_3d+=y_lattice_3d;
					y_iter_3d+=1;

				}

			}

		z_e_3d += (z_lattice_3d);
		z_iter_3d +=1;

		}

		else{


			int y_iter_3d=0;

			double y_e_3d=0.00;
			double y_o_3d=y_lattice_3d/2;

			while(y_iter_3d <number_y)
			{
				double x_e_3d=x_lattice_3d/2;
				double x_o_3d=0.0;

				int x_iter_3d=0;

				if ( y_iter_3d%2==0 )
				{
					while (x_iter_3d < number_x && (x_e_3d+x_lattice_3d)<x_max_3d)
					{

			            Point <dim> *P=new Point <dim>;
			            P->SetXCoord(x_e_3d);
			            P->SetYCoord(y_e_3d);
			            P-> SetZCoord(z_o_3d);

			            Atom <dim> *atom=new Atom <dim>;
			            atom->SetInitialPosition(*P);
			            atom->SetMaterialPosition(*P);
			            atom->SetSpatialPosition(*P);
			            atom->SetID(ID_3d);
			            atom->SetCellID(-1);

			            unrelax_points.push_back(P);
			            unrelax_atoms.push_back(atom);
			            ID_3d+=1;

			            x_e_3d+=x_lattice_3d;
			            x_iter_3d+=1;

					}
					y_e_3d+=y_lattice_3d;
					y_iter_3d+=1;
				}

				else{
						while (x_iter_3d < number_x )
						{

							Point <dim> *P=new Point <dim>;
							P->SetXCoord(x_o_3d);
							P->SetYCoord(y_o_3d);
							P -> SetZCoord(z_o_3d);

							Atom <dim> *atom=new Atom <dim>;
							atom->SetInitialPosition(*P);
							atom->SetMaterialPosition(*P);
							atom->SetSpatialPosition(*P);
							atom->SetID(ID_3d);
							atom->SetCellID(-1);

							unrelax_points.push_back(P);
							unrelax_atoms.push_back(atom);
							ID_3d+=1;

							x_o_3d+=x_lattice_3d;
							x_iter_3d+=1;


						}
						y_o_3d+=y_lattice_3d;
						y_iter_3d+=1;

					}

				}

			z_o_3d += (z_lattice_3d);
			z_iter_3d +=1;

		}


		}
	}

		if(dimension==3 && type=="Si")
		{

			//construct unit cell
			vector < Point <dim>* > unit_cell_atoms;
			vector < Atom <dim>* > atoms_list;

			Point <dim> P0(0.,0.,0.);
			unit_cell_atoms.push_back(&P0);

			Point <dim> P1(lattice_const,0.,0.);
			unit_cell_atoms.push_back(&P1);

			Point <dim> P2(0.,lattice_const,0.);
			unit_cell_atoms.push_back(&P2);

			Point <dim> P3(0.,0.,lattice_const);
			unit_cell_atoms.push_back(&P3);

			Point <dim> P4(lattice_const,lattice_const,0.);
			unit_cell_atoms.push_back(&P4);

			Point <dim> P5(lattice_const,0,lattice_const);
			unit_cell_atoms.push_back(&P5);

            Point <dim> P6(0.,lattice_const,lattice_const);
			unit_cell_atoms.push_back(&P6);

			Point <dim> P7(lattice_const,lattice_const,lattice_const);
			unit_cell_atoms.push_back(&P7);

			//////////////////////////////

			Point <dim> P8(0.5*lattice_const,0.5*lattice_const,0.);
			unit_cell_atoms.push_back(&P8);

			Point <dim> P9(lattice_const,0.5*lattice_const,0.5*lattice_const);
			unit_cell_atoms.push_back(&P9);

			Point <dim> P10(0.5*lattice_const,0.,0.5*lattice_const);
			unit_cell_atoms.push_back(&P10);

			Point <dim> P11(0.5*lattice_const,lattice_const,0.5*lattice_const);
			unit_cell_atoms.push_back(&P11);

			Point <dim> P12(0,0.5*lattice_const,0.5*lattice_const);
			unit_cell_atoms.push_back(&P12);

			Point <dim> P13(0.5*lattice_const,0.5*lattice_const,lattice_const);
			unit_cell_atoms.push_back(&P13);

			//////////////////////////////

			Point <dim> P14(0.25*lattice_const,0.25*lattice_const,0.25*lattice_const);
			unit_cell_atoms.push_back(&P14);

			Point <dim> P15(0.75*lattice_const,0.75*lattice_const,0.25*lattice_const);
			unit_cell_atoms.push_back(&P15);

			Point <dim> P16(0.75*lattice_const,0.25*lattice_const,0.75*lattice_const);
			unit_cell_atoms.push_back(&P16);

			Point <dim> P17(0.25*lattice_const,0.75*lattice_const,0.75*lattice_const);
			unit_cell_atoms.push_back(&P17);


			int ID_3d=0;
			for (int i=0; i<unit_cell_atoms.size(); ++i)
			{

				Atom <dim> *atom=new Atom <dim>;

				atom->SetInitialPosition(*unit_cell_atoms[i]);
				atom->SetMaterialPosition(*unit_cell_atoms[i]);
				atom->SetSpatialPosition(*unit_cell_atoms[i]);
				atom->SetID(ID_3d);
				atom->SetCellID(-1);

				unrelax_atoms.push_back(atom);
				ID_3d+=1;

			}


		for (int i=0; i<number_x; ++i)
		{
			for (int j=0; j<number_y; ++j)
			{
				for (int k=0; k<number_z; ++k)
				{
					if(i==0 && j==0 && k!=0)
					{
						for (int point=0; point< unit_cell_atoms.size(); ++point)
						{

								Point <dim> *P=new Point <dim>;
								if (unit_cell_atoms[point]->GetZCoord()!=0)
								{

									P->SetXCoord(unit_cell_atoms[point]->GetXCoord());
									P->SetYCoord(unit_cell_atoms[point]->GetYCoord());
									P->SetZCoord(unit_cell_atoms[point]->GetZCoord()+k*lattice_const);

									Atom <dim> *atom=new Atom <dim>;

									atom->SetInitialPosition(*P);
									atom->SetMaterialPosition(*P);
									atom->SetSpatialPosition(*P);
									atom->SetID(ID_3d);
									atom->SetCellID(-1);

									unrelax_atoms.push_back(atom);
									ID_3d+=1;

								}

						}
					}

					if(i==0 && j!=0 && k==0)
					{

						for (int point=0; point< unit_cell_atoms.size(); ++point)
						{

								if (unit_cell_atoms[point]->GetYCoord()!=0)
								{
									Point <dim> *P=new Point <dim>;

									P->SetXCoord(unit_cell_atoms[point]->GetXCoord());
									P->SetYCoord(unit_cell_atoms[point]->GetYCoord()+j*lattice_const);
									P->SetZCoord(unit_cell_atoms[point]->GetZCoord());
									Atom <dim> *atom=new Atom <dim>;

									atom->SetInitialPosition(*P);
									atom->SetMaterialPosition(*P);
									atom->SetSpatialPosition(*P);
									atom->SetID(ID_3d);
									atom->SetCellID(-1);

									unrelax_atoms.push_back(atom);
									ID_3d+=1;

								}



						}

					}

					if(i==0 && j!=0 && k!=0)
					{
						for (int point=0; point< unit_cell_atoms.size(); ++point)
						{

								if (unit_cell_atoms[point]->GetYCoord()!=0 && unit_cell_atoms[point]->GetZCoord()!=0)
								{
									Point <dim> *P=new Point <dim>;

									P->SetXCoord(unit_cell_atoms[point]->GetXCoord());
									P->SetYCoord(unit_cell_atoms[point]->GetYCoord()+j*lattice_const);
									P->SetZCoord(unit_cell_atoms[point]->GetZCoord()+k*lattice_const);
									Atom <dim> *atom=new Atom <dim>;

									atom->SetInitialPosition(*P);
									atom->SetMaterialPosition(*P);
									atom->SetSpatialPosition(*P);
									atom->SetID(ID_3d);
									atom->SetCellID(-1);

									unrelax_atoms.push_back(atom);
									ID_3d+=1;

								}

						}

					}

					if(i!=0 && j==0 && k==0)
					{

						for (int point=0; point< unit_cell_atoms.size(); ++point)
						{

								if (unit_cell_atoms[point]->GetXCoord()!=0)
								{
									Point <dim> *P=new Point <dim>;

									P->SetXCoord(unit_cell_atoms[point]->GetXCoord()+i*lattice_const);
									P->SetYCoord(unit_cell_atoms[point]->GetYCoord());
									P->SetZCoord(unit_cell_atoms[point]->GetZCoord());
									Atom <dim> *atom=new Atom <dim>;

									atom->SetInitialPosition(*P);
									atom->SetMaterialPosition(*P);
									atom->SetSpatialPosition(*P);
									atom->SetID(ID_3d);
									atom->SetCellID(-1);

									unrelax_atoms.push_back(atom);
									ID_3d+=1;

								}
						}

					}

					if(i!=0 && j!=0 && k!=0)
					{

						for (int point=0; point< unit_cell_atoms.size(); ++point)
							{

									if (unit_cell_atoms[point]->GetXCoord()!=0 && unit_cell_atoms[point]->GetYCoord()!=0 && unit_cell_atoms[point]->GetZCoord()!=0)
									{
										Point <dim> *P=new Point <dim>;

										P->SetXCoord(unit_cell_atoms[point]->GetXCoord()+i*lattice_const);
										P->SetYCoord(unit_cell_atoms[point]->GetYCoord()+j*lattice_const);
										P->SetZCoord(unit_cell_atoms[point]->GetZCoord()+k*lattice_const);
										Atom <dim> *atom=new Atom <dim>;

										atom->SetInitialPosition(*P);
										atom->SetMaterialPosition(*P);
										atom->SetSpatialPosition(*P);
										atom->SetID(ID_3d);
										atom->SetCellID(-1);

										unrelax_atoms.push_back(atom);
										ID_3d+=1;

									}
							}

					}


					if(i!=0 && j!=0 && k==0)
					{

						for (int point=0; point< unit_cell_atoms.size(); ++point)
							{

									if (unit_cell_atoms[point]->GetXCoord()!=0 && unit_cell_atoms[point]->GetYCoord()!=0 )
									{
										Point <dim> *P=new Point <dim>;

										P->SetXCoord(unit_cell_atoms[point]->GetXCoord()+i*lattice_const);
										P->SetYCoord(unit_cell_atoms[point]->GetYCoord()+j*lattice_const);
										P->SetZCoord(unit_cell_atoms[point]->GetZCoord());
										Atom <dim> *atom=new Atom <dim>;

										atom->SetInitialPosition(*P);
										atom->SetMaterialPosition(*P);
										atom->SetSpatialPosition(*P);
										atom->SetID(ID_3d);
										atom->SetCellID(-1);

										unrelax_atoms.push_back(atom);
										ID_3d+=1;

									}
							}

					}

					if(i!=0 && j==0 && k!=0)
					{

						for (int point=0; point< unit_cell_atoms.size(); ++point)
							{

									if (unit_cell_atoms[point]->GetXCoord()!=0 && unit_cell_atoms[point]->GetZCoord()!=0 )
									{
										Point <dim> *P=new Point <dim>;

										P->SetXCoord(unit_cell_atoms[point]->GetXCoord()+i*lattice_const);
										P->SetYCoord(unit_cell_atoms[point]->GetYCoord());
										P->SetZCoord(unit_cell_atoms[point]->GetZCoord()+k*lattice_const);
										Atom <dim> *atom=new Atom <dim>;

										atom->SetInitialPosition(*P);
										atom->SetMaterialPosition(*P);
										atom->SetSpatialPosition(*P);
										atom->SetID(ID_3d);
										atom->SetCellID(-1);

										unrelax_atoms.push_back(atom);
										ID_3d+=1;

									}
							}

					}
				}

			}

		}



// uncomment it to generate 3d position file

		ofstream file("dump.unrelaxed_config_silicon");

		int id=0;

    	if (file.is_open())
    	{

    		file <<"ITEM: TIMESTEP" <<endl;
    		file <<"0" <<endl;
    		file <<"ITEM: NUMBER OF ATOMS" << endl;
    		file << unrelax_atoms.size() << endl;
    		file <<"ITEM: BOX BOUNDS mm mm pp"<<endl;

    		file << "0.0" << " " << "1.0" << endl;
    		file << "0.0" << " " << "1.0" << endl;
    		file << "0.0" << " " << "1.0" << endl;

    		file <<"ITEM: ATOMS id type x y z" <<endl;

    		for (int i=0; i< unrelax_atoms.size(); ++i)
    		{

    			    Point <dim> atom_material_position=unrelax_atoms[i]->GetMaterialPosition();

        			double materialp_x=atom_material_position.GetXCoord();
        			double materialp_y=atom_material_position.GetYCoord();
        			double materialp_z=atom_material_position.GetZCoord();

            		file << id<< " " << "4" <<" " << materialp_x<<" "
            		<< setprecision(5)<< materialp_y << " "<< materialp_z <<endl;

        			id+=1;

    		}

    		file << "EOF" <<endl;

    	}

        else
    	{
    	    cerr << "Couldn't open unrelaxed_config.xtl for writing." << std::endl;
    	}



	}

	return(unrelax_atoms);


}

template vector < Atom <2>* > UnrelaxedConfigGenerator<2>(int number_x, int number_y, int number_z, double lattice_const, int dimension, string type);
template vector < Atom <3>* > UnrelaxedConfigGenerator<3>(int number_x, int number_y, int number_z, double lattice_const, int dimension, string type);
