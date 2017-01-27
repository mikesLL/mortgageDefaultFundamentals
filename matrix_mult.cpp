#include "headers.h"

void matrix_mult(double (&c)[3], double (&A)[3][3], double (&v0)[3], double (&v1)[3], int niter ) {

	int i, j, n; //niter;

	double v0_tmp[3];

	for (i = 0; i < 3; i++) {
		v0_tmp[i] = v0[i];
	}

	for (n = 0; n < niter; n++) {

		for (i = 0; i < 3; i++) {

			v1[i] = c[i];

			for (j = 0; j < 3; j++) {
				v1[i] += A[i][j] * v0_tmp[j];
			}
		}

		for (i = 0; i < 3; i++) {
			v0_tmp[i] = v1[i];
		}
				
	}

	return;
}