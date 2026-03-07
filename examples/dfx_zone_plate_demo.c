/*!
 *  \file   dfx_zone_plate_demo.c
 *  \brief  An example generating classic zone-plate test pattern in linear space and storing it as an sRGB bitmap image.
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
#include <stdlib.h>
#include "dfx.h"

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

/*!
 *  \brief Generate zone-plate test pattern. 
 * 
 *  Using formula:
 * 
 *     g(x) = 0.5 + 0.5 * sin( km * r^2 / (2*rm) ) * (0.5 + 0.5 * tanh( (rm - r) / w))
 * 
 *   where
 *     r   - distance from the center of the image
 *     rm  - maximum radius of the zone plate pattern
 *     km  - maximum spatial frequency in the zone plate pattern
 *     w   - controls the width of the transition band
 * 
 *  \param[in,out]  y       - image where to place the pattern;
 *  \param[in]      width   - image width in pixels;
 *  \param[in]      height  - image height in pixels;
 *  \param[in]      p       - padding parameter
 * 
 *  \return         DFX_SUCCESS if success, or an error code otherwise.
 */
static int generate_zone_plate(float *y, int width, int height, int p)
{
	float r, rm, km, w, g;
	int i, j, di, dj;

	/* check args: */
	if (y == NULL || width <= 0 || height <= 0 || p < 0) 
		return DFX_INVARG;

	/* set parameters: */
	km = (float)M_PI;									/* touch Nyquist at boundary*/
	rm = (float)( min(width, height) / 2.0f);
	w = rm / 10.0f;						

	/* generate zone plate pattern: */
	for (j = 0; j < height; j++) for (i = 0; i < width; i++) {
		/* compute radius: */
		di = i - width / 2;
		dj = j - height / 2;
		r = sqrtf((float)(di * di + dj * dj));
		/* compute zone plate pattern intensity at (i,j): */
		g = 0.735f * (0.5f + 0.5f * sinf(km * r * r / (2 * rm)) * (0.5f + 0.5f * tanhf((rm - r) / w)));
		/* store the result: */
		y[(p+j) * (width + 2*p) + p + i] = g;
	}

	/* pad image and return: */
	if (p > 0) pad_plane(y, width, height, p, PAD_REFLECT);
	return DFX_SUCCESS;
}

/*!
 *  \brief An example generating gray-scale gradient images, explaining the conversion between linear and sRGB spaces.
 */
int main()
{
	/* assumed image geometry: */
	int height = 480;
	int width = 480;
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

	/* generate zone-plate image: */
	generate_zone_plate(Y, width, height, 0);

	/* convert it to grayscale image: */
	luminance_to_grayscale_image(Y, R, G, B, width, height, 0);

	/* convert to srgb and store it as a file: */
	linear_to_srgb(sRGB, R, G, B, width, height, 0);
	write_bitmap("dfx-zone-plate.bmp", sRGB, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB);
	free_image(R,G,B);
	free_plane(Y);

	return 0;
}
