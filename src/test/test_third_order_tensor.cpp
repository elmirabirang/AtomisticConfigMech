/*
 * test_third_order_tensor.cpp
 *
 *  Created on: Jun 16, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      University of Erlangen-Nürnberg
 */

#include "third_order_tensor.h"
#include <iostream>
#include <vector>

using namespace std;


int main(int argc, char **argv)
{

	ThirdTensor <double> tensor(3,3,3,0.0);

	ThirdTensor <double> tensor1(3,3,3,1.0);

	for (int i=0; i <3 ; i++)
	{
		for(int j=0; j< 3 ; j++)
		{
			cout << tensor(i,j,0) << "\n";
		}

	}


}





