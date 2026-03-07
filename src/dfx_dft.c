/*!
 *  \file   dfx_dft.c
 *  \brief  DFX library: 2D Discrete Fourier Transform (DFT).
 * 
 *  This module provodes a very basic, separable, matrix multiply-based implementation of the 2D DFT transform. 
 *  Father methods for DFT computation can be found in module dfx_fft.c.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "dfx.h"

/******************
 *
 * Internal functions supporting computation of 2D DFT:
 *
 *   gen_dft_coeffs()     - generate DTF coefficients
 *   free_dft_coeffs()    - frees memory used by DFT coefficients
 * 
 *   dft_row()            - 1D DFT for row (real input)
 *   dft_column()         - 1D DFT for column (complex input)
 */

/*! 
 *  \brief Initialize arrays of transform coefficients for 2D DFT computation.
 * 
 *  \param[in,out] pC, pS - pointers to arrays of precomputed transform coefficients
 * 
 *  \returns       DFX_X error codes
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

	/* allocate arays */
	if ((C = (float*)malloc(N * N * sizeof(float))) == NULL) return DFX_NOMEM;
	if ((S = (float*)malloc(N * N * sizeof(float))) == NULL) { free(C); return DFX_NOMEM; }

	/* generate transform coefficients: */
	for (k = 0; k < N; k++) for (n = 0; n < N; n++) {
		C[k * N + n] = cosf((float)(2 * M_PI * k * n)/(float)N);
		S[k * N + n] = sinf((float)(2 * M_PI * k * n)/(float)N);
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
 *  \param[out] reZ, imZ  - pointer to a row in a temp image to contain real/imaginnary parts of the row-processed spectrum
 *  \param[in]  width     - image width 
 *  \param[in]  C, S      - precomputed transform coefficients
 *
 *  \returns    DFX_X error codes
 */
static int dft_row(float* f, float* reZ, float* imZ, int width, float *C, float *S)
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
 *  \param[in]  reZ, imZ  - pointer to a column in a temp image containing real/imaginnary parts of row-processed spectrum
 *  \param[out] reF, imF  - pointer to a column in an image to contain final 2D DFT spectrum
 *  \param[in]  width     - image width
 *  \param[in]  height    - image width
 *  \param[in]  C, S      - precomputed transform coefficients
 *
 *  \returns    DFX_X error codes
 */
static int dft_column(float* reZ, float* imZ, float* reF, float* imF, int width, int height, float* C, float* S)
{
	int k, n;  /* full width of a padded image */

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
 *  \param[in]      f       - input image (real)
 *  \param[out]     g       - output image (real)
 *  \param[in]      height  - image height
 *  \param[in]      width   - image width
 *  \param[in]      p       - padding parameter
 *
 *  \returns        DFX_X error codes
 */
static int dft_preprocess(float* f, float* g, int width, int height, int p)
{
	int x, y;

	/* check parameters: */
	if (f == NULL || width < 0 || height < 0 || p < 0)
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
 *  \brief Two-dimenstional Discrete Fourier Transform over an image plane
 *
 *  \param[in]  f          - input image plane (real)
 *  \param[out] reF, imF   - real/imaginnary parts of Fourier spectrum of input image
 *  \param[in]  height     - image height 
 *  \param[in]  width      - image width 
 *  \param[in]  p          - padding parameter
 *  \param[in]  center_dc  - shift DC to the center of the 2D spectral image 
 *
 *  \returns    DFX_X error codes
 */
int dft_plane(float* f, float* reF, float* imF, int width, int height, int p, int center_dc)
{
	float *reZ=NULL, *imZ=NULL;				/* temporary plane holding the results of row-filtering */
	float *C_row=NULL, *S_row=NULL;			/* transform coefficients for rows */
	float *C_col = NULL, * S_col = NULL;	/* transform coefficients for columns */
	float *f_in;                            /* pointer to an input image to be used for transformation */
	
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

/*
 *  \brief Compute magnitude of a 2D DFT spectrum for an image. 
 *
 *  \param[in]  reF, imF  - real/imaginary parts of Fourier spectrum of input image
 *  \param[out] magF      - magnitude of the Fourier spectrum
 *  \param[in]  height    - image height
 *  \param[in]  width     - image width
 *  \param[in]  p         - padding parameter
 *
 *  \returns    DFX_X error codes
*/
int dft_magnitude(float* reF, float* imF, float* magF, int width, int height, int p)
{
	float re, im, mag;
	int x, y;

	/* check parameters: */
	if (reF == NULL || imF == NULL || magF == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* compute magnitude: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		re = reF[(y + p) * (width + 2 * p) + x + p];
		im = reF[(y + p) * (width + 2 * p) + x + p];
		mag = sqrtf(re*re + im*im);
		magF[(y + p) * (width + 2 * p) + x + p] = mag;
	}

	/* success */
	return DFX_SUCCESS;
}

/*
 *  \brief Compute phase of a 2D DFT spectrum for an image.
 *
 *  \param[in]  reF, imF  - real/imaginary parts of Fourier spectrum of input image
 *  \param[out] phaseF    - phase of the Fourier spectrum
 *  \param[in]  height    - image height
 *  \param[in]  width     - image width
 *  \param[in]  p         - padding parameter
 *
 *  \returns    DFX_X error codes
*/
int dft_phase(float* reF, float* imF, float* phaseF, int width, int height, int p)
{
	float re, im, phase;
	int x, y;

	/* check parameters: */
	if (reF == NULL || imF == NULL || phaseF == NULL || width < 0 || height < 0 || p < 0)
		return DFX_INVARG;

	/* compute phase: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		re = reF[(y + p) * (width + 2 * p) + x + p];
		im = reF[(y + p) * (width + 2 * p) + x + p];
		phase = atan2f(im, re);
		phaseF[(y + p) * (width + 2 * p) + x + p] = phase;
	}

	/* success */
	return DFX_SUCCESS;
}

/* dfx_dft.c -- end of file */
