/*
 * generate_waves.cpp
 *
 *  Created on: May 8, 2019
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Chair of Applied Mechanics
 *      University of ERlangen-Nuremberg
 */


#include <vector>
#include "math.h"
#include <fstream>
#include <string>
#include <stdio.h>

using namespace std;

int EvenOdd(int cell_number)
{
   int chunk=0;

   if (cell_number % 2 == 0)
      {

	   chunk=(cell_number-3)/2;

      }

   else
      {

	   chunk=(cell_number-2)/2;

      }

   return (chunk);

}

void GenerateWaves (int division_x, int division_y, int division_z, int dim)
{

	ofstream file("/calculate/elmira/Implementations/PyQC_C++/pyQC _C++ version/tests/test.data");

	vector < vector <int> > waves;

	int chunk_x_low= EvenOdd(division_x);
	int chunk_x_up= division_x-EvenOdd(division_x);

	int chunk_y_low=EvenOdd(division_y);
	vector <string> colors;
	colors.push_back("blue");
	colors.push_back("red");
	colors.push_back("green");
	colors.push_back("yellow");
	colors.push_back("orange");
	colors.push_back("violet");
	colors.push_back("pink");
	colors.push_back("black");
	colors.push_back("cyan");
	int color_iter=0;


    if (file.is_open())
    {
    switch (dim)
    {

       case(2):
	   {

    	    vector < vector <int> > waves;

    	    for (int j=0; j<= chunk_y_low; ++j)

    	    {


    		   for (int i=0; i<=chunk_x_low; ++i)
    		   {

    			   vector <int> cells;

    			   int pair1=i+j;
    			   int pair2=(division_x-i)+(division_y-j);
    			   int pair3=i+(division_y-j);
    			   int pair4=(division_x-i)+j;

    			   cells.push_back(pair1);
    			   cells.push_back(pair2);
    			   cells.push_back(pair3);
    			   cells.push_back(pair4);
    			   waves.push_back(cells);

    			   file << i << " " << j << " " <<colors[color_iter]<< "\n";
    			   file << division_x-i << " " << division_y-j << " " << colors[color_iter] <<"\n";
    			   file << i << " " << division_y-j  << " " << colors[color_iter]<<"\n";
    			   file << division_x-i << " " << j << " "<< colors[color_iter]  <<"\n";
    			   color_iter+=1;

    		   }

    	   }

	   }

       break;


       case(3):
	   {

    	    vector < vector <int> > waves;

       	    for (int j=0; j<= chunk_y_low; ++j)

       	    {

       		   for (int i=0; i<=chunk_x_low; ++i)
       		   {

       			   vector <int> cells_even;
       			   for (int z=0; z<=division_z; z+=2 )

       			   {


        			   int pair1=i+j+z;
        			   int pair2=(division_x-i)+(division_y-j)+z;
        			   int pair3=i+(division_y-j)+z;
        			   int pair4=(division_x-i)+j+z;

        			   cells_even.push_back(pair1);
        			   cells_even.push_back(pair2);
        			   cells_even.push_back(pair3);
        			   cells_even.push_back(pair4);

        			   file << i << " " << j << " " << z << " " << colors[color_iter]<<"\n";
        			   file << division_x-i << " " << division_y-j << " " << z <<" " << colors[color_iter]<< "\n";
        			   file << i << " " << division_y-j << " " << z <<" " << colors[color_iter] <<"\n";
        			   file << division_x-i << " " << j << " " << z <<" " << colors[color_iter] <<"\n";



       			   }

       			   color_iter+=1;
       			   waves.push_back(cells_even);

       			   vector <int> cells_odd;
       			   for (int z=1; z<=division_z; z+=2)
       			   {

        			   int pair1=i+j+z;
        			   int pair2=(division_x-i)+(division_y-j)+z;
        			   int pair3=i+(division_y-j)+z;
        			   int pair4=(division_x-i)+j+z;

        			   cells_odd.push_back(pair1);
        			   cells_odd.push_back(pair2);
        			   cells_odd.push_back(pair3);
        			   cells_odd.push_back(pair4);

        			   file << i << " " << j << " " << z << " " << colors[color_iter]<<"\n";
        			   file << division_x-i << " " << division_y-j << " " << z << " " << colors[color_iter]<<"\n";
        			   file << i << " " << division_y-j << " " << z << " " << colors[color_iter]<<"\n";
        			   file << division_x-i << " " << j << " " << z << " " << colors[color_iter]<<"\n";



       			   }

       			   color_iter+=1;
       			   waves.push_back(cells_odd);

       		   }


       	   }

	   }

    }
    }


}




