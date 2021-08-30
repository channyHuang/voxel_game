#include "LeastSquareSolver.h"

namespace LeastSquareSolver
{
    float determinant(
        float a, float b, float c,
        float d, float e, float f,
        float g, float h, float i)
	{
		return a * e * i
			+ b * f * g
			+ c * d * h
			- a * f * h
			- b * d * i
			- c * e * g;
	}

    Vector3 solve3x3(const float* A, const float b[3]) {
		auto det = determinant(
			A[0 * 3 + 0], A[0 * 3 + 1], A[0 * 3 + 2],
			A[1 * 3 + 0], A[1 * 3 + 1], A[1 * 3 + 2],
			A[2 * 3 + 0], A[2 * 3 + 1], A[2 * 3 + 2]);
		//Oh-oh - small determinant
		if (abs(det) <= 1e-12) {
			return Vector3(NAN, NAN, NAN);
		}

		return Vector3{
			determinant(
				b[0],    A[0 * 3 + 1], A[0 * 3 + 2],
				b[1],    A[1 * 3 + 1], A[1 * 3 + 2],
				b[2],    A[2 * 3 + 1], A[2 * 3 + 2]),

			determinant(
				A[0 * 3 + 0], b[0],    A[0 * 3 + 2],
				A[1 * 3 + 0], b[1],    A[1 * 3 + 2],
				A[2 * 3 + 0], b[2],    A[2 * 3 + 2]),

			determinant(
				A[0 * 3 + 0], A[0 * 3 + 1], b[0]   ,
				A[1 * 3 + 0], A[1 * 3 + 1], b[1]   ,
				A[2 * 3 + 0], A[2 * 3 + 1], b[2])

		} / det;
	}

    Vector3 leastSquares(size_t N, const Vector3* A, const float* b)
	{
		if (N == 3) {
            const float A_mat[3 * 3] = {
				A[0].x, A[0].y, A[0].z,
				A[1].x, A[1].y, A[1].z,
				A[2].x, A[2].y, A[2].z,
			};
			return solve3x3(A_mat, b);
		}

        float At_A[3][3];
        float At_b[3];

		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
                float sum = 0;
				for (size_t k = 0; k < N; ++k) {
					sum += A[k][i] * A[k][j];
				}
				At_A[i][j] = sum;
			}
		}

		for (int i = 0; i < 3; ++i) {
            float sum = 0;

			for (size_t k = 0; k < N; ++k) {
				sum += A[k][i] * b[k];
			}

			At_b[i] = sum;
		}


		/*
		// Improve conditioning:
        float offset = 0.0001;
		At_A[0][0] += offset;
		At_A[1][1] += offset;
		At_A[2][2] += offset;
		 */

        static_assert(sizeof(At_A) == 9 * sizeof(float), "pack");

		return solve3x3(&At_A[0][0], At_b);
	}
}
