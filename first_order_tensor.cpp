/*
 * first_order_tensor.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This source file provides class of first order tensors and corresponding operators on it.
 */

#include "first_order_tensor.h"
#include "matrix.h"
#include "matrix.cpp"


template <typename T>
FirstTensor<T>::FirstTensor()
{}


template <typename T>
FirstTensor <T>::~FirstTensor()
{}

template <typename T>
FirstTensor <T>::FirstTensor(int number_rows, T& initial_value)
{

	this->nrows=number_rows;

	Vec.resize(nrows,initial_value);


}

template <typename T>
FirstTensor <T> FirstTensor<T>::operator+(FirstTensor<T> &rhs_vector)
{

	FirstTensor result(nrows, 0.0);

	for (int i=0; i<nrows; ++i)
	{

		result(i)=this->Vec[i]+rhs_vector(i);

	}

	return (result);

}

template <typename T>
FirstTensor <T> FirstTensor <T>::operator-(FirstTensor <T> &rhs_vector)
{

	FirstTensor result(nrows, 0.0);

	for (int i=0; i<nrows; ++i)
	{

		result(i)=this->Vec[i]-rhs_vector(i);

	}

	return(result);

}

template <typename T>
FirstTensor <T>& FirstTensor <T>::operator+=(FirstTensor<T> &rhs_vector)
{

	for (int i=0; i< nrows; ++i)
	{

		this->Vec[i]+=rhs_vector(i);

	}

	return (this->Vec);

}

template <typename T>
FirstTensor <T>& FirstTensor <T>::operator-=(FirstTensor <T> &rhs_vector)
{

	for (int i=0; i< nrows; ++i)
	{

	   this->Vec[i]-=rhs_vector(i);

	}

	return (this->Vec);

}

template <typename T>
FirstTensor <T> FirstTensor <T>::operator*(double factor)
{

	FirstTensor result(nrows, 0.0);

	for (int i=0; i< nrows; ++i)
	{

		result(i)=this->Vec[i]*factor;

	}

	return (result);

}

template <typename T>
double FirstTensor <T>::inner_product(FirstTensor<T> &vec1, FirstTensor<T> &vec2)
{

	double result=0.;

	for (int i=0; i<nrows; ++i)
	{

		result+=vec1(i)*vec2(i);

	}

	return (result);

}

template <typename T>
SecondTensor <T> FirstTensor <T>::outer_product(FirstTensor <T> &vec1, FirstTensor <T> &vec2)
{

    //call lhs_tesnor
	SecondTensor <T> lhs_tensor (nrows, nrows, 0.0);

	for (int i=0; i<nrows; ++i)
	{
		for (int j=0; j<nrows; ++j)
		{

			lhs_tensor(i,j)=vec1(i)*vec2(j);

		}

	}

	return (lhs_tensor);

}

template <typename T>
FirstTensor <T> FirstTensor <T>::cross_product (FirstTensor <T> &vec1, FirstTensor <T> &vec2)
{

	FirstTensor result(nrows, 0.0);

	result(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2);
	result(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3);
	result(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1);

	return (result);

}
