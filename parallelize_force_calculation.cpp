/*
 * parallelize_force_calculation.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. candidate
 *      Chair of Applied Mechanics
 *      University of Erlangen-Nuremberg
 */

#include <iostream>
#include "cell.h"
#include "cell.cpp"
#include "point.h"
#include "force.h"
#include "force.cpp"
#include "point.h"
#include "bond.h"
#include "bond.cpp"
#include <vector>
#include <omp.h>
#include "generate_cells.cpp"
#include "unrelaxed_config.cpp"
#include "atoms_in_cell.cpp"
#include "neighbors_cell_based.cpp"
#include "cell_neighbors.cpp"

int main ()
{

	int num_cells_x=20;
	int num_cells_y=20;

    vector < Cell <2>* > cells = GenerateCells <2> (2.2, 1.1, num_cells_x, num_cells_y, 0);

    vector < Cell <2>* > cells_ndim [10][10];
    int iter_cells=0;

    for (int i=0; i < (num_cells_x)/2 ; ++i)
    {
    	for (int j=0; j < (num_cells_y)/2 ; ++j)
    	{

    		Cell <2>* cell=cells[i+j];

    		cells_ndim[i][j].push_back(cells[iter_cells]);

    		iter_cells+=1;

    	}

    }


//    for (int i=0; i < num_cells_x ; ++i)
//    {
//    	for (int j=0; j < num_cells_y ; ++j)
//    	{
//
//    		int cell_index=cells_ndim[i][j][0]->GetCellIndex();
//    		cout << "cell index: " << cell_index << endl;
//
//    	}
//
//    }

    vector < Atom <2>* > unrelax_atoms;
    unrelax_atoms=UnrelaxedConfigGenerator <2> (num_cells_x*2.2, num_cells_y*2.2, 0, 1.1, 2);

	Force <2> force(1.0,1.0);

	CellsNeighbors <2> (cells,2.2);
	AtominCell (unrelax_atoms, cells, 2.2);
	AtomNeighbors <2> (cells,2.2);

	vector < Cell <2>* > wave;
	vector < vector < Cell <2>* > > waves;
	int num_wave_x=10/3;
	int num_wave_y=10/3;
	int size_waves=num_wave_x*num_wave_y;


	for (int iter_waves=0; iter_waves< size_waves; ++iter_waves)
	{

		int iter_cells_x=0;

		vector < Cell <2>* > wave;

		if (iter_waves> 2 && iter_waves < 6 )
		{

			iter_cells_x=1;

		}

		if (iter_waves> 5 && iter_waves< size_waves )
		{

			iter_cells_x=2;

		}

		cout << iter_cells_x <<endl;

		while (iter_cells_x < num_cells_x/2)
		{

			int iter_cells_y=iter_waves;

			if (iter_waves> 2 && iter_waves < 6 )
			{

				iter_cells_y=iter_waves-3;

			}

			if (iter_waves> 5 && iter_waves< size_waves )
			{

				iter_cells_y=iter_waves-6;

			}

			while (iter_cells_y < num_cells_y/2)
			{

				wave.push_back(cells_ndim[iter_cells_x][iter_cells_y][0]);

				iter_cells_y+=3;

			}

			iter_cells_x+=3;

		}

		waves.push_back(wave);

	}


	int iter=0;
	int iter_2=0;
//	int size_waves=waves.size();
//	int size_wave=wave.size();

	while ( iter < waves.size())

	{
		vector < Cell <2>* > wave=waves[iter];
		int iter_2=0;

		cout << "wave: " << iter << endl;

		while (iter_2 < wave.size())
		{

			Cell <2> * cell= wave[iter_2];
			int cell_index=cell->GetCellIndex();

			cout << " cell index: " << cell_index << endl;
			iter_2+=1;

		}

		iter+=1;

	}

//    vector < Atom <2>* >::iterator At=unrelax_atoms.begin();
//
//    #pragma omp parallel for num_threads(3)
//    for (int i=0; i < (unrelax_atoms.size()); ++i)
//    {
//
//    	int atomID=(*At)->GetID();
//    	int this_thread = omp_get_thread_num();
//
//    	cout <<  this_thread << endl;
//
//    	At++;
//
//    }

	return 0;
}




