/*!
 *  \file   dfx_dct_demo.c
 *  \brief  An example testing the fast DCT-II and DCT-IV algorithms of dfx_dct.c.
 *
 *  The program compares the outputs of dct2_rec() and dct4_rec() against the direct evaluation
 *  of the (unnormalized) transform definitions
 *
 *      DCT-II:  X[k] = sum_{n=0}^{N-1} x[n] cos( pi (2n+1) k / (2N) ),
 *      DCT-IV:  X[k] = sum_{n=0}^{N-1} x[n] cos( pi (2n+1)(2k+1) / (4N) ),
 *
 *  computed in double precision, for random inputs. It covers the short-length modules
 *  (N = 2, 3, 4, 5, 7, 9, 15) and composite lengths N = 2^m * N1 up to 3840, including the block
 *  sizes used by image and video codecs (4, 8, 16, 32, 64) and the common image dimensions 480,
 *  768, 1024, 1920, and 3840. (720 = 16 * 45 and 1080 = 8 * 135 are not supported, since 45 and
 *  135 are not among the short-length modules.) It also checks that every length from 2 to 64 is
 *  either rejected or transformed correctly. The program exits with a nonzero code if any error
 *  exceeds a tolerance consistent with single-precision arithmetic.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.00
 *  \date    November 8, 2026
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dfx.h"

#define N_MAX       3840		/* largest transform length tested */
#define N_TRIALS    5			/* random inputs per length */
#define TOLERANCE   1e-5		/* maximum allowed relative error */

/*!
 *  \brief Reference unnormalized DCT-II: X[k] = sum_n x[n] cos(pi (2n+1) k / (2N)).
 */
static void dct2_ref(const float* x, double* X, int N)
{
	int n, k;
	for (k = 0; k < N; k++) {
		double s = 0.;
		for (n = 0; n < N; n++)
			s += x[n] * cos(M_PI * (double)((2 * n + 1) * k) / (2.0 * N));
		X[k] = s;
	}
}

/*!
 *  \brief Reference unnormalized DCT-IV: X[k] = sum_n x[n] cos(pi (2n+1)(2k+1) / (4N)).
 */
static void dct4_ref(const float* x, double* X, int N)
{
	int n, k;
	for (k = 0; k < N; k++) {
		double s = 0.;
		for (n = 0; n < N; n++)
			s += x[n] * cos(M_PI * (double)((2 * n + 1) * (2 * k + 1)) / (4.0 * N));
		X[k] = s;
	}
}

/*!
 *  \brief Fill a vector with uniform random values in [-0.5, 0.5].
 */
static void random_vector(float* x, int n)
{
	int k;
	for (k = 0; k < n; k++) x[k] = (float)rand() / (float)RAND_MAX - 0.5f;
}

/*!
 *  \brief Largest error between a computed transform and the reference, relative to the largest reference magnitude.
 */
static double rel_error(const float* y, const double* Y, int n)
{
	int k; double e = 0, scale = 0;
	for (k = 0; k < n; k++) {
		if (fabs(y[k] - Y[k]) > e) e = fabs(y[k] - Y[k]);
		if (fabs(Y[k]) > scale) scale = fabs(Y[k]);
	}
	return e / (scale > 0 ? scale : 1.0);
}

/*!
 *  \brief Test one transform for one length; returns 1 if the test failed.
 */
static int test_dct(const char* name, int (*transform)(float*, float*, int), void (*reference)(const float*, double*, int), int n,
	float* x, float* y, float* tmp, double* X)
{
	double worst = 0, e;
	int q, err;

	for (q = 0; q < N_TRIALS; q++) {
		random_vector(x, n);
		reference(x, X, n);
		memcpy(y, x, n * sizeof(float));
		if ((err = transform(y, tmp, n)) != DFX_SUCCESS) {
			printf("%-8s n=%5d: returned error code %d\n", name, n, err);
			return 1;
		}
		if ((e = rel_error(y, X, n)) > worst) worst = e;
	}
	printf("%-8s n=%5d: max relative error = %.2e %s\n", name, n, worst, worst < TOLERANCE ? "ok" : "FAILED");
	return worst >= TOLERANCE;
}

/*!
 *  \brief DCT test/demo program.
 */
int main(void)
{
	/* work arrays: */
	static float x[N_MAX], y[N_MAX], tmp[N_MAX];
	static double X[N_MAX];

	/* transform lengths: short modules, then composite lengths 2^m * N1 */
	static const int short_sizes[] = {2, 3, 4, 5, 7, 9, 15};
	static const int composite_sizes[] = {6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 30, 32, 36, 40, 48, 56, 60, 64, 72, 80, 96,
		112, 120, 128, 144, 160, 192, 224, 240, 256, 288, 320, 384, 448, 480, 512, 576, 640, 768, 960, 1024, 1280, 1536, 1920, 2048, 3840};
	int i, fails = 0;

	srand(777);

	printf("Short-length transform modules:\n");
	for (i = 0; i < (int)(sizeof(short_sizes) / sizeof(short_sizes[0])); i++) {
		fails += test_dct("dct2_rec", dct2_rec, dct2_ref, short_sizes[i], x, y, tmp, X);
		fails += test_dct("dct4_rec", dct4_rec, dct4_ref, short_sizes[i], x, y, tmp, X);
	}

	printf("\nComposite-length transforms:\n");
	for (i = 0; i < (int)(sizeof(composite_sizes) / sizeof(composite_sizes[0])); i++) {
		fails += test_dct("dct2_rec", dct2_rec, dct2_ref, composite_sizes[i], x, y, tmp, X);
		fails += test_dct("dct4_rec", dct4_rec, dct4_ref, composite_sizes[i], x, y, tmp, X);
	}

	/* every length must be either rejected or transformed correctly: */
	printf("\nAll lengths 2..64 through dct2_rec() and dct4_rec(); rejected as unsupported: ");
	for (i = 2; i <= 64; i++) {
		int rc2, rc4; double e2 = 0, e4 = 0;
		random_vector(x, i);
		dct2_ref(x, X, i);
		memcpy(y, x, i * sizeof(float));
		if ((rc2 = dct2_rec(y, tmp, i)) == DFX_SUCCESS) e2 = rel_error(y, X, i);
		dct4_ref(x, X, i);
		memcpy(y, x, i * sizeof(float));
		if ((rc4 = dct4_rec(y, tmp, i)) == DFX_SUCCESS) e4 = rel_error(y, X, i);
		if (rc2 == DFX_NOTSUP && rc4 == DFX_NOTSUP) {
			printf("%d ", i);
		} else if (rc2 != rc4 || e2 >= TOLERANCE || e4 >= TOLERANCE) {
			printf("\nn=%d: dct2_rec returned %d (error %.1e), dct4_rec returned %d (error %.1e): FAILED\n", i, rc2, e2, rc4, e4);
			fails++;
		}
	}
	printf("\n");

	printf("\n%s\n", fails ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
	return fails ? 1 : 0;
}

/* dfx_dct_demo.c -- end of file */
