/*
 * first_order_tensor.h
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This header file provides class of first order tensors and corresponding operators on it.
 */

#ifndef TESTS_FIRST_ORDER_TENSOR_H_
#define TESTS_FIRST_ORDER_TENSOR_H_

#include <vector>
#include "matrix.h"

using namespace std;

template <typename T>
class FirstTensor
{

private:

   int nrows;
   vector<T> Vec;

  public:

	FirstTensor();
	~FirstTensor();
	FirstTensor(int nrows, T &initial_values);

	FirstTensor <T> operator+(FirstTensor <T> &rhs_vector);
	FirstTensor <T> &operator+=(FirstTensor <T> &rhs_vector);

	FirstTensor <T> operator-(FirstTensor <T> &rhs_vector);
	FirstTensor <T> &operator-=(FirstTensor <T> &rhs_vector);

	FirstTensor <T> operator*(double rhs_factor);

	double inner_product(FirstTensor <T> &vec1, FirstTensor <T> &vec2);
	SecondTensor <T> outer_product(FirstTensor <T> &vec1, FirstTensor <T> &vec2);
	FirstTensor <T> cross_product(FirstTensor <T> &vec1, FirstTensor <T> &vec2);

};


#endif /* TESTS_FIRST_ORDER_TENSOR_H_ */
