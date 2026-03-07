/*!
 *  \file   dfx_dft_demo.c
 *  \brief  An example program demonstrating how to use 2D DFT transform.
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
#include <stdio.h>
#include "dfx.h"

/*!
 *  \brief An example program, showing how to use 2D DFT transform and related functions.
 */
int main(int argc, char* argv[])
{
	/* image parameters: */
	int height = 0, width = 0, ppm_x = 0, ppm_y = 0, x, y;

	/* images and planes: */
	unsigned char* sRGB_in = NULL, *sRGB_out = NULL;
	float *R = NULL, *G = NULL, *B = NULL, *Y = NULL;
	float *reDFT = NULL, *imDFT = NULL, *magDFT = NULL, *phaseDFT = NULL;
	float magDC;

	/* check parameters and read sRGB image: */
	if (argc < 2) {
		printf("Usage: dfx_dft_demo <bitmap image>\n");
		return 1;
	}
	if (read_bitmap(argv[1], &sRGB_in, &width, &height, &ppm_x, &ppm_y) != DFX_SUCCESS) {
		printf("Cannot open input file: %s\n", argv[1]);
		return 1;
	}

	/* allocate working images and planes: */
	if (alloc_srgb_image(&sRGB_out, width, height) != DFX_SUCCESS
		|| alloc_image(&R, &G, &B, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&Y, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&reDFT, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&imDFT, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&magDFT, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&phaseDFT, width, height, 0) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free allocated memory & exit: */
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB_out != NULL) free_srgb_image(sRGB_out);
		if (R != NULL && G != NULL && B != NULL) free_image(R, G, B);
		if (Y != NULL) free_plane(Y);
		if (reDFT != NULL) free_plane(reDFT);
		if (imDFT != NULL) free_plane(imDFT);
		if (magDFT != NULL) free_plane(magDFT);
		if (phaseDFT != NULL) free_plane(phaseDFT);
		return 1;
	}

	/* Extract luminance: */
	srgb_to_linear(sRGB_in, R, G, B, width, height, 0);
	linear_to_luminance(Y, R, G, B, width, height, 0);

	/* Visualize luminance: */
	linear_to_srgb(sRGB_out, Y, Y, Y, width, height, 0);
	write_bitmap("dft_demo_luminance.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* compute DFT: */
	dft_plane(Y, reDFT, imDFT, width, height, 0, 1);
	dft_magnitude(reDFT, imDFT, magDFT, width, height, 0);
	dft_phase(reDFT, imDFT, phaseDFT, width, height, 0);

	/* Visualize DFT magnitude: */
	magDC = magDFT[(height/2) * width + width / 2];
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) magDFT[y * width + x] /= magDC;  /* normalize to fit in [0,1] range */
	linear_to_srgb(sRGB_out, magDFT, magDFT, magDFT, width, height, 0);
	write_bitmap("dft_demo_magnitude.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* Visualize DFT log-magnitude: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) magDFT[y * width + x] = logf(1+magDFT[y * width + x] * magDC) / logf(magDC);  /* normalize to fit in [0,1] range */
	linear_to_srgb(sRGB_out, magDFT, magDFT, magDFT, width, height, 0);
	write_bitmap("dft_demo_log_magnitude.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* Visualize DFT phase: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) phaseDFT[y * width + x] = phaseDFT[y * width + x] / ((float)(2.0*M_PI)) + 0.5f; /* normalize to fit in [0,1] range */
	linear_to_srgb(sRGB_out, phaseDFT, phaseDFT, phaseDFT, width, height, 0);
	write_bitmap("dft_demo_phase.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R, G, B);
	free_plane(Y);
	free_plane(reDFT);
	free_plane(imDFT);
	free_plane(magDFT);
	free_plane(phaseDFT);

	return 0;
}
