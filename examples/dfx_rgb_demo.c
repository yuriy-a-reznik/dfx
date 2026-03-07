/*!
 *  \file   dfx_rgb_demo.c
 *  \brief  An example program demonstrating basic operations with RGB images, and bitmap files. 
 *
 *  This example shows how to split an sRGB image into component images containing its R, G, B, and luminance channels.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#include <stdio.h>
#include <math.h>
#include "dfx.h"

 /*!
  *  \brief An example generating gray-scale gradient images, explaining the conversion between linear and sRGB spaces.
  */
int main(int argc, char *argv[])
{
	/* image parameters: */
	int height = 0, width = 0, ppm_x = 0, ppm_y = 0;

	/* linear and srgb images/planes: */
	unsigned char *sRGB_in = NULL, *sRGB_out = NULL;
	float *R = NULL, *G = NULL, *B = NULL, *Y = NULL, *Z = NULL;

	/* check parameters and read sRGB image: */
	if (argc < 2) {
		printf("Usage: dfx_split_rgb <bitmap image>\n");
		return 1;
	}
	if (read_bitmap(argv[1], &sRGB_in, &width, &height, &ppm_x, &ppm_y) != DFX_SUCCESS) {
		printf("Cannot open input file: %s\n", argv[1]);
		printf("This program only works with uncompressed 24-bit RGB bitmap files\n");
		return 1;
	}

	/* allocate working images and planes: */
	if (alloc_srgb_image(&sRGB_out, width, height) != DFX_SUCCESS
		|| alloc_image(&R, &G, &B, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&Y, width, height, 0) != DFX_SUCCESS 
		|| alloc_plane(&Z, width, height, 0) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free allocated memory & exit: */
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB_out != NULL) free_srgb_image(sRGB_out);
		if (R != NULL && G != NULL && B != NULL) free_image(R, G, B);
		if (Y != NULL) free_plane(Y);
		if (Z != NULL) free_plane(Z);
		return 1;
	}

	/* Extract linear RGB channels: */
	srgb_to_linear(sRGB_in, R, G, B, width, height, 0);
	
	/* Extract luminance: */
	linear_to_luminance(Y, R, G, B, width, height, 0);

	/* Produce zero-plane: */
	zero_plane(Z, width, height, 0);

	/*
	 * Produce R-channel image:
	 */
	linear_to_srgb(sRGB_out, R, Z, Z, width, height, 0);
	write_bitmap("red.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/*
	 * Produce G-channel image:
	 */
	linear_to_srgb(sRGB_out, Z, G, Z, width, height, 0);
	write_bitmap("green.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/*
	 * Produce B-channel image:
	 */
	linear_to_srgb(sRGB_out, Z, Z, B, width, height, 0);
	write_bitmap("blue.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/*
	 * Produce gray-scale image:
	 */
	linear_to_srgb(sRGB_out, Y, Y, Y, width, height, 0);
	write_bitmap("grayscale.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R, G, B);
	free_plane(Y);
	free_plane(Z);

	return 0;
}
