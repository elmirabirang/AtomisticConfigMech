/*
 * matrix.h
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This header file contains required functions and operations of a second order tensor.
 */
#pragma once

#include "point.h"
#include "third_order_tensors.h"
#include <vector>

/*I define a template class with double placeholder
 * Which operations I need on matrices
 * inner product
 * outer product with a vector
 * inner product with a vector
 * Define unity tensor
 *
 *
 *
 */

template <typename T>
class SecondTensor
{

private:

   unsigned nrows;
   unsigned ncols;
   std::vector <std::vector <T> > matrix;

  public:

   SecondTensor();
   SecondTensor(unsigned _nrows, unsigned _ncols, const T& _initial_values);
   SecondTensor(const SecondTensor<T>& rhs);
   virtual ~SecondTensor();

   SecondTensor<T>& operator=(const SecondTensor<T>& rhs);

	//Matrix-Matrix Operations

	SecondTensor <T> operator+(const SecondTensor <T> &rhs_matrix);
	SecondTensor <T>& operator+=(const SecondTensor <T> &rhs_matrix);

	SecondTensor <T> operator-(const SecondTensor <T> &rhs_matrix);
	SecondTensor <T>& operator-=(const SecondTensor <T> &rhs_matrix);

	SecondTensor <T> operator*(const SecondTensor <T> &rhs_matrix);
	SecondTensor <T>& operator*=(const SecondTensor <T> &rhs_matrix);

	SecondTensor <T> transpose();

	//Matrix-Scalar Operations

	SecondTensor <T> operator+(const T& rhs_factor);
	SecondTensor <T> operator-(const T& rhs_factor);
	SecondTensor <T> operator*(const T& rhs_factor);
	SecondTensor <T> operator/(const T& rhs_factor);

	//Matrix-Vector Operations

	std::vector<T> operator*(const std::vector<T>& rhs_vec);
	std::vector<T> diag_vec();

	//Access the individual element

	T& operator()(const unsigned& row, const unsigned& col);

	const T& operator()(const unsigned& row, const unsigned& column) const;

	//Access the row and column size
	unsigned get_rows() const;
	unsigned get_cols() const;

	ThirdTensor<T> outer_product_MaxVec(SecondTensor <T>& lhs_matrix, std::vector <T>& rhs_vector);
	ThirdTensor<T> outer_product_VecMax (std::vector <T>& lhs_vector, SecondTensor <T>& rhs_matrix);
	SecondTensor <T> dot(const ThirdTensor <T>& lhs_vec, const std::vector <T>& rhs_vec);

};
