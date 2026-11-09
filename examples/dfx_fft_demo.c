/*!
 *  \file   dfx_fft_demo.c
 *  \brief  An example testing the fast Fourier transform algorithms of dfx_fft.c.
 *
 *  The program compares the outputs of the radix-2 decimation-in-time and decimation-in-frequency
 *  FFTs and of the prime-factor FFT against the direct evaluation of the DFT definition
 *
 *      X[k] = sum_{n=0}^{N-1} x[n] exp(-2 pi i n k / N),
 *
 *  computed in double precision, for random complex inputs and all supported transform lengths
 *  up to 5040. It prints the largest relative error for each length and algorithm, checks that
 *  fft() either rejects or correctly transforms every length from 2 to 100, and exits with a
 *  nonzero code if any error exceeds a tolerance consistent with single-precision arithmetic.
 *  It also times the direct and fast computations for a few lengths.
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
#include <time.h>
#include "dfx.h"

#define N_MAX       5040		/* largest transform length tested */
#define N_TRIALS    3			/* random inputs per length */
#define TOLERANCE   1e-5		/* maximum allowed relative error */

/*!
 *  \brief Reference DFT, evaluated directly from the definition in double precision.
 */
static void dft_ref(const float* xr, const float* xi, double* Xr, double* Xi, int n)
{
	int k, j;
	for (k = 0; k < n; k++) {
		double sr = 0, si = 0;
		for (j = 0; j < n; j++) {
			double a = -2.0 * M_PI * (double)((k * j) % n) / (double)n;   /* k*j reduced modulo n: the twiddle factors are periodic */
			double c = cos(a), s = sin(a);
			sr += xr[j] * c - xi[j] * s;
			si += xr[j] * s + xi[j] * c;
		}
		Xr[k] = sr; Xi[k] = si;
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
 *  \brief Largest error between a computed spectrum and the reference, relative to the largest reference magnitude.
 */
static double rel_error(const float* yr, const float* yi, const double* Yr, const double* Yi, int n)
{
	int k; double e = 0, scale = 0, d, m;
	for (k = 0; k < n; k++) {
		d = hypot(yr[k] - Yr[k], yi[k] - Yi[k]);
		m = hypot(Yr[k], Yi[k]);
		if (d > e) e = d;
		if (m > scale) scale = m;
	}
	return e / (scale > 0 ? scale : 1.0);
}

/*!
 *  \brief Test one algorithm for one transform length; returns 1 if the test failed.
 */
static int test_fft(const char* name, int (*transform)(float*, float*, int), int n,
	float* xr, float* xi, float* yr, float* yi, double* Xr, double* Xi)
{
	double worst = 0, e;
	int q, err;

	for (q = 0; q < N_TRIALS; q++) {
		random_vector(xr, n);
		random_vector(xi, n);
		dft_ref(xr, xi, Xr, Xi, n);
		memcpy(yr, xr, n * sizeof(float));
		memcpy(yi, xi, n * sizeof(float));
		if ((err = transform(yr, yi, n)) != DFX_SUCCESS) {
			printf("%-8s n=%5d: returned error code %d\n", name, n, err);
			return 1;
		}
		if ((e = rel_error(yr, yi, Xr, Xi, n)) > worst) worst = e;
	}
	printf("%-8s n=%5d: max relative error = %.2e %s\n", name, n, worst, worst < TOLERANCE ? "ok" : "FAILED");
	return worst >= TOLERANCE;
}

/*!
 *  \brief Compare running times of the direct DFT and the FFT for a given length.
 */
static void time_fft(int n, float* xr, float* xi, double* Xr, double* Xi)
{
	clock_t t0, t1, t2;
	int reps = 20, q;

	random_vector(xr, n);
	random_vector(xi, n);
	t0 = clock();
	for (q = 0; q < reps; q++) dft_ref(xr, xi, Xr, Xi, n);
	t1 = clock();
	for (q = 0; q < reps; q++) fft(xr, xi, n);   /* repeatedly transforming the output is fine for timing purposes */
	t2 = clock();
	printf("n=%5d: direct DFT %8.3f ms, fft() %8.3f ms per transform\n", n,
		1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC / reps,
		1000.0 * (double)(t2 - t1) / (double)CLOCKS_PER_SEC / reps);
}

/*!
 *  \brief FFT test/demo program.
 */
int main(void)
{
	/* work arrays: */
	static float xr[N_MAX], xi[N_MAX], yr[N_MAX], yi[N_MAX];
	static double Xr[N_MAX], Xi[N_MAX];

	/* lengths for the prime factor algorithm: all products of mutually prime factors from {2,3,4,5,7,8,9,16} */
	static const int twos[] = {1, 2, 4, 8, 16}, threes[] = {1, 3, 9}, fives[] = {1, 5}, sevens[] = {1, 7};
	static const int timing_sizes[] = {256, 720, 1024, 5040};

	int i, j, k, l, n, fails = 0;

	srand(12345);

	printf("Radix-2 FFTs (decimation in time and decimation in frequency):\n");
	for (n = 2; n <= 4096; n *= 2) {
		fails += test_fft("fft_dit", fft_dit, n, xr, xi, yr, yi, Xr, Xi);
		fails += test_fft("fft_dif", fft_dif, n, xr, xi, yr, yi, Xr, Xi);
	}

	printf("\nPrime factor FFT (Winograd modules of lengths 2,3,4,5,7,8,9,16):\n");
	for (i = 0; i < 5; i++) for (j = 0; j < 3; j++) for (k = 0; k < 2; k++) for (l = 0; l < 2; l++) {
		n = twos[i] * threes[j] * fives[k] * sevens[l];
		if (n > 1)
			fails += test_fft("fft_pfa", fft_pfa, n, xr, xi, yr, yi, Xr, Xi);
	}

	/* every length must be either rejected or transformed correctly: */
	printf("\nAll lengths 2..100 through fft(); rejected as unsupported: ");
	for (n = 2; n <= 100; n++) {
		random_vector(xr, n);
		random_vector(xi, n);
		dft_ref(xr, xi, Xr, Xi, n);
		memcpy(yr, xr, n * sizeof(float));
		memcpy(yi, xi, n * sizeof(float));
		if (fft(yr, yi, n) == DFX_NOTSUP) {
			printf("%d ", n);
		} else if (rel_error(yr, yi, Xr, Xi, n) >= TOLERANCE) {
			printf("\nfft(n=%d) accepted the length but produced a wrong result: FAILED\n", n);
			fails++;
		}
	}
	printf("\n");

	printf("\nTiming (direct O(N^2) evaluation vs. fft()):\n");
	for (i = 0; i < (int)(sizeof(timing_sizes) / sizeof(timing_sizes[0])); i++)
		time_fft(timing_sizes[i], xr, xi, Xr, Xi);

	printf("\n%s\n", fails ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
	return fails ? 1 : 0;
}

/* dfx_fft_demo.c -- end of file */
