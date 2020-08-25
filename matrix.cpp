/*
 * matrix.cpp
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This source file contains required functions and operations of a second order tensor.
 */
#include "matrix.h"
#include "third_order_tensor.h"

template <typename T>
SecondTensor <T>::SecondTensor()
{}

template <typename T>
SecondTensor <T>::SecondTensor(unsigned _nrows, unsigned _ncols, const T& _initial_values)
{

	matrix.resize(_nrows);

	for (int i=0; i < matrix.size(); ++i)
	{

		matrix[i].resize(_ncols, _initial_values);

	}

	nrows=_nrows;
	ncols=_ncols;


}

template <typename T>
SecondTensor <T>::SecondTensor(const SecondTensor<T>& rhs)
{
	matrix=rhs.matrix;
	nrows=rhs.get_rows();
	ncols=rhs.get_cols();
}

template <typename T>
SecondTensor <T>::~SecondTensor()
{}


template <typename T>
SecondTensor <T>& SecondTensor<T>::operator=(const SecondTensor<T> &rhs_matrix)
{
	if (&rhs_matrix==this)
		return *this;

	unsigned new_rows=rhs_matrix.get_rows();
	unsigned new_cols=rhs_matrix.get_cols();

	matrix.resize(new_rows);

	for (unsigned i=0; i<matrix.size(); i++)
	{
		matrix[i].resize(new_cols);
	}

	for (unsigned i=0; i<new_rows; i++)
	{
		for (unsigned j=0; j <new_cols; j++)
		{

			matrix[i][j]=rhs_matrix(i,j);

		}


	}
	nrows=new_rows;
	ncols=new_cols;

	return *this;

}

template <typename T>
SecondTensor <T> SecondTensor <T>::operator+(const SecondTensor <T> &rhs_matrix)
{

	SecondTensor<T> result(nrows, ncols, 0.0);

	for (unsigned i=0; i<nrows; i++)
	{
		for (unsigned j=0; j<ncols; j++)
		{

			result(i,j)=this->matrix[i][j]+rhs_matrix(i,j);

		}

	}

	return result;

}

template <typename T>
SecondTensor <T>& SecondTensor <T>::operator+=(const SecondTensor<T> &rhs_matrix)
{

	unsigned nrows=rhs_matrix.get_rows();
	unsigned ncols=rhs_matrix.get_cols();



	for (int i=0; i< nrows; i++)
	{

		for (int j=0; j<ncols; j++)
		{

			this->matrix[i][j]+=rhs_matrix(i,j);

		}
	}

	return *this;

}

template <typename T>
SecondTensor <T> SecondTensor <T>::operator-(const SecondTensor <T> &rhs_matrix)
{

	unsigned nrows=rhs_matrix.get_rows();
	unsigned ncols=rhs_matrix.get_cols();

	SecondTensor result(nrows, ncols, 0.0);

	for (unsigned i=0; i<nrows; i++)
	{
		for (unsigned j=0; j<ncols; j++)
		{

			result(i,j)=this->matrix[i][j]-rhs_matrix(i,j);

		}

	}

	return result;

}

template <typename T>
SecondTensor <T>& SecondTensor <T>::operator-=(const SecondTensor<T> &rhs_matrix)
{

	unsigned nrows=rhs_matrix.get_rows();
	unsigned ncols=rhs_matrix.get_cols();



	for (int i=0; i< nrows; i++)
	{

		for (int j=0; j<ncols; j++)
		{

			this->matrix[i][j]-=rhs_matrix(i,j);

		}
	}

	return *this;

}

template <typename T>
SecondTensor <T> SecondTensor <T>::operator*(const SecondTensor& rhs_matrix)
{

	unsigned nrows=rhs_matrix.get_rows();
	unsigned ncols=rhs_matrix.get_cols();

	SecondTensor result(nrows, ncols, 0.0);

	for (unsigned i=0; i<nrows; i++)
	{
		for (unsigned j=0; j<ncols; j++)
		{
			for (unsigned k=0; k<nrows; k++)
			{

				result(i,j)+=this->matrix[i][k]*rhs_matrix(k,j);


			}

		}


	}

	return result;

}

template <typename T>
SecondTensor <T>& SecondTensor<T>::operator*=(const SecondTensor<T>& rhs_matrix)
{

	SecondTensor result=(*this) * rhs_matrix;
	(*this)=result;
	return *this;

}

template<typename T>
SecondTensor<T> SecondTensor<T>::transpose()
{

  SecondTensor result(nrows, ncols, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      result(i,j) = this->matrix[j][i];
    }
  }

  return result;
}


template<typename T>
SecondTensor<T> SecondTensor<T>::operator+(const T& rhs) {
  SecondTensor result(nrows, ncols, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      result(i,j) = this->matrix[i][j] + rhs;
    }
  }

  return result;
}

// Matrix/scalar subtraction
template<typename T>
SecondTensor<T> SecondTensor<T>::operator-(const T& rhs) {
  SecondTensor result(nrows, ncols, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      result(i,j) = this->matrix[i][j] - rhs;
    }
  }

  return result;
}

// Matrix/scalar multiplication
template<typename T>
SecondTensor<T> SecondTensor<T>::operator*(const T& rhs) {
  SecondTensor result(nrows, ncols, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      result(i,j) = this->matrix[i][j] * rhs;
    }
  }

  return result;
}

// Matrix/scalar division
template<typename T>
SecondTensor<T> SecondTensor<T>::operator/(const T& rhs) {
  SecondTensor result(nrows, ncols, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {
      result(i,j) = this->matrix[i][j] / rhs;
    }
  }

  return result;
}

template<typename T>
std::vector<T> SecondTensor<T>::operator*(const std::vector<T>& rhs)
{
  std::vector<T> result(rhs.size(), 0.0);

  for (unsigned i=0; i<nrows; i++) {
    for (unsigned j=0; j<ncols; j++) {

      result[i] += this->matrix[i][j] * rhs[j];
    }
  }

  return result;
}

template<typename T>
std::vector<T> SecondTensor<T>::diag_vec() {
  std::vector<T> result(nrows, 0.0);

  for (unsigned i=0; i<nrows; i++) {
    result[i] = this->matrix[i][i];
  }

  return result;
}

template <typename T>
const T& SecondTensor<T>::operator()(const unsigned& row, const unsigned& col) const
{

	return this->matrix[row][col];

}

template<typename T>
T& SecondTensor<T>::operator()(const unsigned& row, const unsigned& col)
{

  return this->matrix[row][col];

}

template <typename T>
unsigned SecondTensor <T>::get_rows() const
{
	return this->nrows;
}

template <typename T>
unsigned SecondTensor <T>::get_cols() const
{
	return this->ncols;
}



template <typename T>
ThirdTensor <T> SecondTensor <T>::outer_product_MaxVec(SecondTensor <T>& lhs_matrix, std::vector <T>& rhs_vector)
{

	int nlayers=rhs_vector.size();

	ThirdTensor<T> result(nrows,ncols,nlayers,0.0);

	for (int i=0; i<nrows; i++)
	{
		for (int j=0; j<ncols; j++)
		{

			for (int k=0; k<nlayers; k++)
			{

				result(i,j,k)=lhs_matrix(i,j)*rhs_vector[k];

			}

		}

	}

	return result;

}

template <typename T>
ThirdTensor <T> SecondTensor <T>::outer_product_VecMax (std::vector <T>& lhs_vector, SecondTensor <T>& rhs_matrix)
{
	
	int nrows=lhs_vector.size();
	int ncols=rhs_matrix.get_rows();
	int nlayers=rhs_matrix.get_cols();

	ThirdTensor<T> result(nrows, ncols,nlayers, 0.0);

	for (int k=0; k < nlayers; k++)

	{
		for (int i=0; i < nrows; i++)
		{

			for (int j=0; j< ncols; ++j)
			{
				result(i,j,k)=rhs_matrix(j,k)*lhs_vector[i];
			}

		}

	}

	return result;

}

template <typename T>
SecondTensor <T> SecondTensor <T>::dot(const ThirdTensor <T> &lhs_third_order,const std::vector <T> &rhs_vec)
{

	unsigned nrows=lhs_third_order.get_rows();
	unsigned ncols=lhs_third_order.get_cols();
	unsigned nlayers=lhs_third_order.get_layers();

	SecondTensor <T> result(3,3,0.);

	for(int k=0; k<nlayers; k++)
	{
		for (int i=0; i<nrows; i++)
		{
			for (int j=0; j<ncols; j++)
			{

				result(i,j)+=lhs_third_order(i,j,k)*rhs_vec[k];

			}

		}

	}

	return (result);

}










































