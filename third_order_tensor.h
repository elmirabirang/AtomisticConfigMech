/*
 * third_order_tensor.h
 *
 *  Created on: Jun 15, 2020
 *      Author: S.Elmira Birang.O
 *      Mechanical Engineering Ph.D. Candidate
 *      Central Institute for Scientific Computing
 *      University of Erlangen-Nürnberg
 *      This header file provides class of third order tensors and corresponding operators on it.
 */

#ifndef TESTS_THIRD_ORDER_TENSOR_H_
#define TESTS_THIRD_ORDER_TENSOR_H_

#include <vector>
#include "matrix.h"


template <typename T>
class ThirdTensor
{

	private:

		int nrows;
		int ncols;
		int nlayers;

		std::vector <std::vector <std::vector <T> > > third_tensor;

	public:

		ThirdTensor(unsigned _nrows, unsigned _ncols, unsigned _nlayers, T _initial_values);

		ThirdTensor (const ThirdTensor<T>& rhs);

		virtual ~ThirdTensor();

		ThirdTensor<T>& operator=(const ThirdTensor<T>& rhs);

		ThirdTensor<T> operator+(const ThirdTensor<T>& rhs);
		ThirdTensor<T>& operator+=(const ThirdTensor<T>& rhs);
		ThirdTensor<T> operator-(const ThirdTensor<T>& rhs);
		ThirdTensor<T>& operator-=(const ThirdTensor<T>& rhs);


		const T& operator()(const unsigned& row, const unsigned& column, const unsigned& layer) const;
		T& operator()(const unsigned& row, const unsigned& col, const unsigned& layer);

		unsigned get_rows() const;
		unsigned get_cols() const;
		unsigned get_layers() const;

};



#endif /* TESTS_THIRD_ORDER_TENSOR_H_ */
