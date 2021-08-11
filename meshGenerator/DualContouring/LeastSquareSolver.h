#ifndef LEAST_SQUARE_SOLVER_H
#define LEAST_SQUARE_SOLVER_H

#include "Vector3.h"

namespace LeastSquareSolver
{
	//Various simple low-dimensional linear algebra solvers

	/*
	 Solves for x in  A*x = b.
	 'A' contains the matrix row-wise.
	 'b' and 'x' are column vectors.
	 */
	Vector3 solve3x3(const real* A, const real b[3]);


	/*
	Solves A*x = b for over-determined systems.
	Solves using  At*A*x = At*b   trick where At is the transponate of A.
	Note that each LORD::Vector3 in A is a row in the A-matrix.
	N is the length of A and b.
	*/
	Vector3 leastSquares(size_t N, const Vector3* A, const real* b);
}

#endif