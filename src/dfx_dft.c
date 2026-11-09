/*!
 *  \file   dfx_dft.c
 *  \brief  DFX library: 2D Discrete Fourier Transform (DFT).
 * 
 *  This module provides a very basic, separable, matrix multiply-based implementation of the 2D DFT transform:
 *
 *      F(u,v) = sum_{x=0}^{W-1} sum_{y=0}^{H-1} f(x,y) exp(-2 pi i (u x / W + v y / H)),   u = 0..W-1, v = 0..H-1,
 *
 *  computed as a 1D DFT of every row followed by a 1D DFT of every column of the result.
 *  Each 1D DFT is evaluated directly from its definition using precomputed tables of
 *  cos(2 pi k n / N) and sin(2 pi k n / N), so the cost is O(W H (W + H)).
 *  Faster methods for DFT computation can be found in module dfx_fft.c.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    November 8, 2026
 */

#include <math.h>
#include <stdlib.h>
#include "dfx.h"

/******************
 *
 * Internal functions supporting computation of 2D DFT:
 *
 *   gen_dft_coeffs()     - generate DFT coefficients
 *   free_dft_coeffs()    - frees memory used by DFT coefficients
 * 
 *   dft_row()            - 1D DFT for row (real input)
 *   dft_column()         - 1D DFT for column (complex input)
 */

/*! 
 *  \brief Initialize arrays of transform coefficients for 2D DFT computation.
 * 
 *  \param[out]    pC, pS - pointers to variables to receive the N x N tables C[k][n] = cos(2 pi k n / N), S[k][n] = sin(2 pi k n / N)
 *  \param[in]     N      - transform length
 * 
 *  \returns       DFX_X error code
 */
static int gen_dft_coeffs(float** pC, float** pS, int N)
{
	float* C, * S;
	int k, n;

	/* check args: */
	if (pC == NULL || pS == NULL || N < 0) return DFX_INVARG;

	/* set both pointers to NULL first: */
	*pC = NULL;
	*pS = NULL;

	/* allocate arrays */
	if ((C = (float*)malloc((size_t)N * N * sizeof(float) + 1)) == NULL) return DFX_NOMEM;
	if ((S = (float*)malloc((size_t)N * N * sizeof(float) + 1)) == NULL) { free(C); return DFX_NOMEM; }

	/* generate transform coefficients; the product k*n is reduced modulo N first (the
	 * twiddle factors are periodic), which keeps the argument of cos/sin small and accurate: */
	for (k = 0; k < N; k++) for (n = 0; n < N; n++) {
		double a = 2 * M_PI * (double)((k * n) % N) / (double)N;
		C[k * N + n] = (float)cos(a);
		S[k * N + n] = (float)sin(a);
	}

	/* set pointers to initialized arrays: */
	*pC = C;
	*pS = S;
	return DFX_SUCCESS;
}

/*!
 *  Frees transform coefficients
 */
static int free_dft_coeffs(float* C, float* S)
{
	if (C == NULL || S == NULL) return DFX_INVARG;
	free(C);
	free(S);
	return DFX_SUCCESS;
}

/*!
 *  \brief Apply 1D DFT across a row in an image.
 * 
 *  \param[in]  f         - pointer to a row in an input image (real) 
 *  \param[out] reZ, imZ  - pointers to a row in a temp image to contain real/imaginary parts of the row-processed spectrum
 *  \param[in]  width     - image width 
 *  \param[in]  C, S      - precomputed transform coefficients
 *
 *  \returns    DFX_X error code
 */
static int dft_row(const float* f, float* reZ, float* imZ, int width, const float *C, const float *S)
{
	int k, n;

	/* check parameters: */
	if (f == NULL || reZ == NULL || imZ == NULL || width < 0 || C == NULL || S == NULL)
		return DFX_INVARG;

	/* compute the transform: */
	for (k = 0; k < width; k++)
	{
		/* compute k-th Fourier coefficient: */
		double re = 0, im = 0;
		for (n = 0; n < width; n++) {
			re += f[n] * C[k * width + n]; /* cos(2 * M_PI * k * n / width); */
			im -= f[n] * S[k * width + n]; /* sin(2 * M_PI * k * n / width); */
		}
		reZ[k] = (float)re;
		imZ[k] = (float)im;
	}

	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Apply 1D DFT across a column in an image.
 * 
 *  \param[in]  reZ, imZ  - pointers to a column in a temp image containing real/imaginary parts of row-processed spectrum
 *  \param[out] reF, imF  - pointers to a column in an image to contain final 2D DFT spectrum
 *  \param[in]  width     - full width of the padded planes = distance between consecutive samples of a column
 *  \param[in]  height    - image height
 *  \param[in]  C, S      - precomputed transform coefficients
 *
 *  \returns    DFX_X error code
 */
static int dft_column(const float* reZ, const float* imZ, float* reF, float* imF, int width, int height, const float* C, const float* S)
{
	int k, n;

	/* check parameters: */
	if (reZ == NULL || imZ == NULL || reF == NULL || imF == NULL || width < 0 || height < 0 || C == NULL || S == NULL)
		return DFX_INVARG;

	/* compute the transform: */
	for (k = 0; k < height; k++)
	{
		/* compute k-th Fourier coefficient: */
		double re = 0, im = 0;
		for (n = 0; n < height; n++) 
		{
			re += reZ[n * width] * C[k * height + n]; /* cos(2 * M_PI * k * n / height); */
			re += imZ[n * width] * S[k * height + n]; /* sin(2 * M_PI * k * n / height); */ 
			im += imZ[n * width] * C[k * height + n]; /* cos(2 * M_PI * k * n / height); */
			im -= reZ[n * width] * S[k * height + n]; /* sin(2 * M_PI * k * n / height); */
		}
		reF[k * width] = (float)re;
		imF[k * width] = (float)im;
	}

	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Image pre-processing shifting 2D DFT's DC to the center.
 *
 *  Multiplying f(x,y) by (-1)^(x+y) = exp(i pi (x + y)) shifts its spectrum by half the
 *  image size in each direction, so that the DC term lands at (width/2, height/2).
 *  The shift is exact when width and height are even.
 *
 *  \param[in]      f       - input image (real)
 *  \param[out]     g       - output image (real)
 *  \param[in]      width   - image width
 *  \param[in]      height  - image height
 *  \param[in]      p       - padding parameter
 *
 *  \returns        DFX_X error code
 */
static int dft_preprocess(const float* f, float* g, int width, int height, int p)
{
	int x, y;

	/* check parameters: */
	if (f == NULL || g == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* generate pre-processed image: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) 
	{
		if ((x + y) & 1)
			g[(y + p) * (width + 2 * p) + x + p] = -f[(y + p) * (width + 2 * p) + x + p];
		else 
			g[(y + p) * (width + 2 * p) + x + p] = f[(y + p) * (width + 2 * p) + x + p];
	}

	/* success */
	return DFX_SUCCESS;
}

/*!
 *  \brief Two-dimensional Discrete Fourier Transform over an image plane
 *
 *  \param[in]  f          - input image plane (real)
 *  \param[out] reF, imF   - real/imaginary parts of Fourier spectrum of input image (planes of the same geometry as f)
 *  \param[in]  width      - image width 
 *  \param[in]  height     - image height 
 *  \param[in]  p          - padding parameter
 *  \param[in]  center_dc  - if nonzero, shift DC to the center of the 2D spectral image 
 *
 *  \returns    DFX_X error code
 */
int dft_plane(const float* f, float* reF, float* imF, int width, int height, int p, int center_dc)
{
	float *reZ=NULL, *imZ=NULL;				/* temporary planes holding the results of row transforms */
	float *C_row=NULL, *S_row=NULL;			/* transform coefficients for rows */
	float *C_col = NULL, * S_col = NULL;	/* transform coefficients for columns */
	const float *f_in;                      /* pointer to an input image to be used for transformation */
	
	int x, y;

	/* check parameters: */
	if (f == NULL || reF == NULL || imF == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* allocate temporary image planes & arrays of filter coeffs: */
	if (alloc_plane(&reZ, width, height, p) != DFX_SUCCESS
	 || alloc_plane(&imZ, width, height, p) != DFX_SUCCESS
	 || gen_dft_coeffs(&C_row, &S_row, width) != DFX_SUCCESS
	 || gen_dft_coeffs(&C_col, &S_col, height) != DFX_SUCCESS) 
	{
		/* free buffers & exit with error */
		if (reZ != NULL) free_plane(reZ); 
		if (imZ != NULL) free_plane(imZ); 
		if (C_row != NULL && S_row != NULL) free_dft_coeffs(C_row, S_row);
		if (C_col != NULL && S_col != NULL) free_dft_coeffs(C_col, S_col);
		return DFX_NOMEM;
	}

	/* check if we need to shift DC to center: */
	f_in = f;
	if (center_dc) {
		/* use output array as temporary buffer for preprocessed image: */
		dft_preprocess(f, reF, width, height, p);
		f_in = reF;
	}

	/* transforms rows: */
	for (y = 0; y < height; y++) {
		/* call DFT function on each row (skipping padding): */
		dft_row(&f_in[(y + p) * (width + 2 * p) + p], 
				&reZ[(y + p) * (width + 2 * p) + p], 
				&imZ[(y + p) * (width + 2 * p) + p], 
				width, C_row, S_row);
	}

	/* transforms columns: */
	for (x = 0; x < width; x++) {
		/* call DFT function on each column (skipping padding): */
		dft_column(
			&reZ[p * (width + 2 * p) + x + p],
			&imZ[p * (width + 2 * p) + x + p],
			&reF[p * (width + 2 * p) + x + p],
			&imF[p * (width + 2 * p) + x + p],
			width+2*p, height, C_col, S_col);
	}

	/* free intermediate buffers & exit: */
	free_plane(reZ);
	free_plane(imZ);
	free_dft_coeffs(C_row, S_row);
	free_dft_coeffs(C_col, S_col);

	return DFX_SUCCESS;
}

/*!
 *  \brief Compute magnitude of a 2D DFT spectrum for an image. 
 *
 *  \param[in]  reF, imF  - real/imaginary parts of Fourier spectrum of input image
 *  \param[out] magF      - magnitude of the Fourier spectrum
 *  \param[in]  width     - image width
 *  \param[in]  height    - image height
 *  \param[in]  p         - padding parameter
 *
 *  \returns    DFX_X error code
 */
int dft_magnitude(const float* reF, const float* imF, float* magF, int width, int height, int p)
{
	float re, im, mag;
	int x, y;

	/* check parameters: */
	if (reF == NULL || imF == NULL || magF == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* compute magnitude: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		re = reF[(y + p) * (width + 2 * p) + x + p];
		im = imF[(y + p) * (width + 2 * p) + x + p];
		mag = sqrtf(re*re + im*im);
		magF[(y + p) * (width + 2 * p) + x + p] = mag;
	}

	/* success */
	return DFX_SUCCESS;
}

/*!
 *  \brief Compute phase of a 2D DFT spectrum for an image.
 *
 *  \param[in]  reF, imF  - real/imaginary parts of Fourier spectrum of input image
 *  \param[out] phaseF    - phase of the Fourier spectrum, in the range (-pi, pi]
 *  \param[in]  width     - image width
 *  \param[in]  height    - image height
 *  \param[in]  p         - padding parameter
 *
 *  \returns    DFX_X error code
 */
int dft_phase(const float* reF, const float* imF, float* phaseF, int width, int height, int p)
{
	float re, im, phase;
	int x, y;

	/* check parameters: */
	if (reF == NULL || imF == NULL || phaseF == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* compute phase: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		re = reF[(y + p) * (width + 2 * p) + x + p];
		im = imF[(y + p) * (width + 2 * p) + x + p];
		phase = atan2f(im, re);
		phaseF[(y + p) * (width + 2 * p) + x + p] = phase;
	}

	/* success */
	return DFX_SUCCESS;
}

/* dfx_dft.c -- end of file */
