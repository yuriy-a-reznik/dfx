/*!
 *  \file   dfx_gradient_demo.c
 *  \brief  An example generating gray-scale gradient images, illustrating the differences between the linear RGB and gamma-companded spaces.
 * 
 *  This example generates test images, convertes them to sRGB space an saves them as bitmap files.
 *  This example also shows that 8-bit sRGB format has some limitations. The straightforward rounding to 8-bit valies produces 
 *  visible "banding" effects in these test images. Dithered conversion to sRGB helps to reduce the appearance of banding. 
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.01
 *  \date    March 7, 2026
 */

#include <stdio.h>
#include <math.h>
#include "dfx.h"

 /*!
  *  \brief Generate grayscale gradient image.
  *
  *  \param[in,out]  y       - image where to place the pattern;
  *  \param[in]      width   - image width in pixels;
  *  \param[in]      height  - image height in pixels;
  *  \param[in]      p       - padding parameter
  *  \param[in]		 gamma   - gamma indicator (0 - linear, 1 - power-2.2-law)
  *
  *  \return         DFX_SUCCESS if success, or an error code otherwise.
  */
static int generate_gradient(float* Y, int width, int height, int p, int gamma)
{
	float dx, dy, L;
	int x, y;

	/* check args: */
	if (Y == NULL || width <= 0 || height <= 0 || p < 0)
		return DFX_INVARG;

	/* generate linear gradient luminance image: */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			/* produce luminance as a linear function of x and y */
			dx = (float)x / (width - 1);
			dy = (float)y / (height - 1);
			L = sqrtf((dx*dx + dy*dy)/2);   /* distance from 0,0, normalized to [0,1] range */
			if (gamma) {
				/* apply pow-2.2 gamma: */
				L = powf(L, 2.2f);
			}
			/* store it in luminance image: */
			Y[(p + y) * (width + 2 * p) + p +  x] = L;
		}
	}

	/* pad image and return: */
	if (p > 0) pad_plane(Y, width, height, p, PAD_ZERO);
	return DFX_SUCCESS;
}

/*!
 *  \brief DFX: Gray-scale gradient demo program.
 */
int main()
{
	/* assumed image geometry: */
	int height = 4098;
	int width = 4098;
	int ppm_x = 1024;	// 96 DPI pixel density
	int ppm_y = 1024;

	/* linear and srgb images/planes: */
	unsigned char *sRGB = NULL;
	float *R=NULL, *G=NULL, *B = NULL, *Y=NULL;

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
	 * Generate linear gradient:
	 */

	/* generate linear gradient luminance image: */
	generate_gradient(Y, width, height, 0, 0);

	/* convert it to grayscale image: */
	luminance_to_grayscale_image(Y, R, G, B, width, height, 0);

	/* convert to srgb and store it as a file: */
	linear_to_srgb(sRGB, R, G, B, width, height, 0);
	write_bitmap("linear_gradient.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* convert to srgb by additionally applying error diffusion: */
	linear_to_srgb_dithered(sRGB, R, G, B, width, height, 0);
	write_bitmap("linear_gradient_diffused.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* 
	 * Generate power-2.2-law-companded gradient: 
	 */

	/* generate power-2.2 gradient luminance image: */
	generate_gradient(Y, width, height, 0, 1);

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
