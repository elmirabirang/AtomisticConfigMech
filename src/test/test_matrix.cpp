/*
 * test_matrix.cpp
 *
 *  Created on: Jul 28, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering PhD Candidate
 *      University of Erlangen-Nuremberg
 */

#include "matrix.h"
#include "matrix.cpp"
#include "third_order_tensor.h"
#include "third_order_tensor.cpp"
#include <iostream>
#include <vector>

int main()
{
	SecondTensor<double> obj(3,3,0.);
	SecondTensor<double> mat1(3,3,2.);
	SecondTensor<double> mat2(3,3,-1.);
	std::vector<double> vec(3,1.);


	SecondTensor<double> mat3=mat1*2.0;

	std::vector<double> result(3,2.);


	ThirdTensor <double> maxvec(3,3,3,0.);

	ThirdTensor <double> vecmax(3,3,3,0.);

	maxvec=obj.outer_product_MaxVec(mat1,result);
	vecmax=obj.outer_product_VecMax(result,mat1);

	for (int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{

				std::cout << vecmax(i,j,k) << ",";

			}
			std::cout << std:: endl;
		}
		std::cout << std:: endl;
	}


//	for (int i=0; i<result.size(); i++)
//	{
//
//		std::cout << result[i] << ",";
//
//	}


//	for (int i=0; i<mat3.get_rows(); i++)
//	{
//		for (int j=0; j<mat3.get_cols(); j++)
//		{
//
//			std::cout << mat3(i,j) << ",";
//
//		}
//
//		std::cout << std:: endl;
//
//
//	}



	return 0;
}



