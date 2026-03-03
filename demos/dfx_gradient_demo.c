/*!
 *  \file   dfx_gradient_demo.c
 *  \brief  An example generating gray-scale gradient images, explaining the conversion between linear and sRGB spaces.
 * 
 *  This example also shows that 8-bit sRGB format has limitations. The straightforward rounding to 8-bit valies produces 
 *  "banding" effects visible in this image. Dithered conversion to sRGB helps to reduce the appearance of banding. 
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
int main()
{
	/* assumed image geometry: */
	int height = 1024;
	int width = 1024;
	int ppm_x = 1024;	// 96 DPI pixel density
	int ppm_y = 1024;

	/* linear and srgb images/planes: */
	unsigned char *sRGB = NULL;
	float *R=NULL, *G=NULL, *B = NULL, *Y=NULL;
	
	/* variables: */
	int x, y;
	float L;

	/* allocate images */
	if (alloc_srgb_image(&sRGB, width, height) != DFX_SUCCESS
		|| alloc_image(&R, &G, &B, width, height, 0) != DFX_SUCCESS
		|| alloc_plane(&Y, width, height, 0) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free allocated memory & exit: */
		if (sRGB != NULL) free_srgb_image(sRGB);
		if (R != NULL && G != NULL && B != NULL) free_image(R, G, B);
		if (Y != NULL) free_plane(Y);
		return 1;
	}

	/*
	 * 1st experiment: linear gradient:
	 */

	/* generate linear gradient luminance image: */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			/* produce luminance as a linear function of x and y */
			L = ((float)x / (width - 1)) + ((float)y / (height - 1)) / 2.0f;
			/* store it in luminance image: */
			Y[y * width + x] = L;
		}
	}

	/* convert it to grayscale image: */
	luminance_to_grayscale_image(Y, R, G, B, width, height, 0);

	/* convert to srgb and store it as a file: */
	linear_to_srgb(sRGB, R, G, B, width, height, 0);
	write_bitmap("linear_gradient.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* convert to srgb by additionally applying error diffusion: */
	linear_to_srgb_dithered(sRGB, R, G, B, width, height, 0);
	write_bitmap("linear_gradient_diffused.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* 
	 * 2nd experiment: power-2.2 gradient: 
	 */

	/* generate power-2.2 gradient luminance image: */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			/* produce luminance as a power-2.2 function of x and y */
			L = powf((float)x / (width - 1),2.2f) + powf((float)y / (height - 1), 2.2f) / 2.0f;
			/* store it in luminance image: */
			Y[y * width + x] = L;
		}
	}

	/* convert it to grayscale image: */
	luminance_to_grayscale_image(Y, R, G, B, width, height, 0);

	/* convert to srgb and store it as a file: */
	linear_to_srgb(sRGB, R, G, B, width, height, 0);
	write_bitmap("pow22_gradient.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* convert to srgb by additionally applying error diffusion: */
	linear_to_srgb_dithered(sRGB, R, G, B, width, height, 0);
	write_bitmap("pow22_gradient_diffused.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB);
	free_image(R,G,B);
	free_plane(Y);

	return 0;
}
