/*
 * third_order_tensor.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This source file provides class of third order tensors and corresponding operators on it.
 */

#include "third_order_tensor.h"
#include "matrix.h"
#include "matrix.cpp"
#include <vector>




template <typename T>
ThirdTensor <T>::ThirdTensor(unsigned _number_rows, unsigned _number_columns, unsigned _number_layers, T _initial_val)
{

	third_tensor.resize(_number_rows);


	for (int i=0 ; i <_number_rows ; ++i)
	{

		third_tensor[i].resize(_number_columns);

		for (int j=0; j < _number_columns ; ++j)
		{

			third_tensor[i][j].resize(_number_layers, _initial_val);
		}

	}

	nrows=_number_rows;
	ncols=_number_columns;
	nlayers=_number_layers;

}

template <typename T>
ThirdTensor <T> :: ThirdTensor(const ThirdTensor <T> & rhs)
{
	third_tensor=rhs.third_tensor;
	nrows=rhs.get_rows();
	ncols=rhs.get_cols();
	nlayers=rhs.get_layers();

}


template <typename T>
ThirdTensor <T>::~ThirdTensor()
{}

template <typename T>
ThirdTensor <T>& ThirdTensor <T>::operator =(const ThirdTensor<T>& rhs)
{

	if (&rhs==this)
		return *this;

	unsigned new_rows=rhs.get_rows();
	unsigned new_cols=rhs.get_cols();
	unsigned new_layers=rhs.get_layers();

	third_tensor.resize(new_rows);

	for (unsigned i=0; i<third_tensor.size(); i++)
	{
		third_tensor[i].resize(new_cols);
	}

	for (unsigned i=0; i<new_rows; i++)
	{
		for (unsigned j=0; j <new_cols; j++)
		{

			third_tensor[i][j].resize(new_layers);

		}

	}

	for (unsigned i=0; i<new_rows; i++)
	{
		for (unsigned j=0; j <new_cols; j++)
		{
			for (unsigned k=0; k<new_layers; k++)
			{

			   third_tensor[i][j][k]=rhs(i,j,k);

			}

		}

	}

	nrows=new_rows;
	ncols=new_cols;
	nlayers=new_layers;

	return *this;

}

template <typename T>
ThirdTensor <T> ThirdTensor <T>::operator+(const ThirdTensor <T> &rhs)
{

	unsigned nrows=rhs.get_rows();
	unsigned ncols=rhs.get_cols();
	unsigned nlayers=rhs.get_layers();

	ThirdTensor result(nrows, ncols, nlayers, 0.0);

	for (unsigned i=0; i<nrows; i++)
	{
		for (unsigned j=0; j<ncols; j++)
		{
			for (unsigned k=0; k<nlayers; ++k)
			{

				result(i,j,k)=this->third_tensor[i][j][k]+rhs(i,j,k);

			}

		}

	}

	return result;

}

template <typename T>
ThirdTensor <T>& ThirdTensor <T>::operator+=(const ThirdTensor<T> &rhs)
{

	unsigned nrows=rhs.get_rows();
	unsigned ncols=rhs.get_cols();
	unsigned nlayers=rhs.get_layers();



	for (unsigned i=0; i< nrows; i++)
	{

		for (unsigned j=0; j<ncols; j++)
		{

			for (unsigned k=0; k<nlayers; ++k)
			{
				this->third_tensor[i][j][k]+=rhs(i,j,k);
			}

		}
	}

	return *this;

}


template <typename T>
ThirdTensor <T> ThirdTensor <T>::operator-(const ThirdTensor <T> &rhs)
{

	unsigned nrows=rhs.get_rows();
	unsigned ncols=rhs.get_cols();
	unsigned nlayers=rhs.get_layers();

	ThirdTensor result(nrows, ncols, nlayers, 0.0);

	for (unsigned i=0; i<nrows; i++)
	{
		for (unsigned j=0; j<ncols; j++)
		{
			for (unsigned k=0; k<nlayers; ++k)
			{

				result(i,j,k)=this->third_tensor[i][j][k]-rhs(i,j,k);

			}

		}

	}

	return result;

}

template <typename T>
ThirdTensor <T>& ThirdTensor <T>::operator-=(const ThirdTensor<T> &rhs)
{

	unsigned nrows=rhs.get_rows();
	unsigned ncols=rhs.get_cols();
	unsigned nlayers=rhs.get_layers();

	for (unsigned i=0; i< nrows; i++)
	{

		for (unsigned j=0; j<ncols; j++)
		{

			for (unsigned k=0; k<nlayers; ++k)
			{
				this->third_tensor[i][j][k]-=rhs(i,j,k);
			}

		}
	}

	return *this;

}

template <typename T>
const T& ThirdTensor<T>::operator()(const unsigned& row, const unsigned& col,const unsigned& layer) const
{

	return this->third_tensor[row][col][layer];

}

template<typename T>
T& ThirdTensor<T>::operator()(const unsigned& row, const unsigned& col,const unsigned& layer)
{

  return this->third_tensor[row][col][layer];

}

template <typename T>
unsigned ThirdTensor <T>::get_rows() const
{
	return this->nrows;
}

template <typename T>
unsigned ThirdTensor <T>::get_cols() const
{
	return this->ncols;
}

template <typename T>
unsigned ThirdTensor <T>::get_layers() const
{
	return this->nlayers;
}






