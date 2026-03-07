/*!
 *  \file   dfx_filter_demo.c
 *  \brief  An example program, showing how to use 2D image filtering functions.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#include <stdio.h>
#include "dfx.h"

/*!
 *  \brief An example program, showing how to use 2D DFT transform and related functions.
 */
int main(int argc, char* argv[])
{
	/* image parameters: */
	int height = 0, width = 0, ppm_x = 0, ppm_y = 0;

	/* images and planes: */
	unsigned char* sRGB_in = NULL, *sRGB_out = NULL;
	float *R_in = NULL, *G_in = NULL, *B_in = NULL;
	float* R_out = NULL, * G_out = NULL, * B_out = NULL;

	/* padding in linear form: */
	int p = 10;

	/* filter width parameters: */
	float fc = 0.25;								
	int n = 5;

	/* some basic checks: */
	if (n > p) {
		printf("Error: filter width cannot exceed padding of an image\n");
		return 1;
	}
	if (fc <= 0 || fc > 1.0) {
		printf("Error: incorrect value of cutoff frequency parameter\n");
		return 1;
	}

	/* check program parameters and read sRGB image: */
	if (argc < 2) {
		printf("Usage: dfx_filter_image <bitmap image>\n");
		return 1;
	}
	if (read_bitmap(argv[1], &sRGB_in, &width, &height, &ppm_x, &ppm_y) != DFX_SUCCESS) {
		printf("Cannot open input file: %s\n", argv[1]);
		return 1;
	}

	/* allocate working images and planes: */
	if (alloc_srgb_image(&sRGB_out, width, height) != DFX_SUCCESS
		|| alloc_image(&R_in, &G_in, &B_in, width, height, p) != DFX_SUCCESS
		|| alloc_image(&R_out, &G_out, &B_out, width, height, p) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free memory and exit: */
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB_out != NULL) free_srgb_image(sRGB_out);
		if (R_in != NULL && G_in != NULL && B_in != NULL) free_image(R_in, G_in, B_in);
		if (R_out != NULL && G_out != NULL && B_out != NULL) free_image(R_out, G_out, B_out);
		return 1;
	}

	/* convert to linear RGB: */
	srgb_to_linear(sRGB_in, R_in, G_in, B_in, width, height, p);
	linear_to_srgb(sRGB_out, R_in, G_in, B_in, width, height, p);
	write_bitmap("dfx_filter_original.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply Gaussian filter: */
	zero_image(R_out, G_out, B_out, width, height, p);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_GAUSS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_filter_gauss.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply sinc filter: */
	zero_image(R_out, G_out, B_out, width, height, p);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_SINC, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_filter_sinc.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply Lanczos filter: */
	zero_image(R_out, G_out, B_out, width, height, p);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_LANCZOS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_filter_lanczos.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R_in, G_in, B_in);
	free_image(R_out, G_out, B_out);

	return 0;
}
