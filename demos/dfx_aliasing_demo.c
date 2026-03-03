/*!
 *  \file   dfx_aliasing_demo.c
 *  \brief  An example program, illustrating the effects of aliasing when subsampling is not followed by a low-pass filter.
 *
 *  For most obvious visual effects please use zone-plate.bmp image. 
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

 /**************
  *
  *  Subsamling operations.
  *
  *  These functions have the effect of subsampling an image by a factor of 2, 
  *  and then upsamling it (with zero-fill of missing pixels) prior to applying an interpolation filter. 
  */

 /*!
  *  \brief Fill odd-rows & odd-columns in RGB image plane with zeros.
  */
static int subsample_plane(float* X_in, float* X_out, int width, int height, int p)
{
	int w_lin = width + 2 * p;
	int x, y;

	/* check params: */
	if (X_in == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* zero odd-rows & odd-columns: */
	for (y = -p; y < height + p; y++) {
		if (y & 1) {
			/* whipe out entire row: */
			for (x = -p; x < width + p; x++)
				X_out[(p + y) * w_lin + p + x] = 0;
		}
		else {
			/* whipe out every odd pixel: */
			for (x = -p; x < width + p; x++) {
				if (x & 1)
					X_out[(p + y) * w_lin + p + x] = 0;
				else 
					X_out[(p + y) * w_lin + p + x] = X_in[(p + y) * w_lin + p + x];
			}
		}
	}
	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill odd-rows & odd-columns in RGB image with zeros.
 */
static int subsample_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p)
{
	if (R_in == NULL || G_in == NULL || B_in == NULL || R_out == NULL || G_out == NULL || B_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	subsample_plane(R_in, R_out, width, height, p);
	subsample_plane(G_in, G_out, width, height, p);
	subsample_plane(B_in, B_out, width, height, p);
	return DFX_SUCCESS;
}

/*!
 *  \brief An example program illustrating the effects of aliasing, and their prevention by using proper low-pass filters.
 */
int main(int argc, char* argv[])
{
	/* image parameters: */
	int height = 0, width = 0, ppm_x = 0, ppm_y = 0;

	/* images and planes: */
	unsigned char* sRGB_in = NULL, *sRGB_out = NULL;
	float *R_in = NULL, *G_in = NULL, *B_in = NULL;
	float *R_temp = NULL, * G_temp = NULL, * B_temp = NULL;
	float* R_out = NULL, * G_out = NULL, * B_out = NULL;

	/* padding & filter parameters: */
	int p = 10;
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
		|| alloc_image(&R_temp, &G_temp, &B_temp, width, height, p) != DFX_SUCCESS
		|| alloc_image(&R_out, &G_out, &B_out, width, height, p) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free memory and exit: */
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB_out != NULL) free_srgb_image(sRGB_out);
		if (R_in != NULL && G_in != NULL && B_in != NULL) free_image(R_in, G_in, B_in);
		if (R_temp != NULL && G_temp != NULL && B_temp != NULL) free_image(R_temp, G_temp, B_temp);
		if (R_out != NULL && G_out != NULL && B_out != NULL) free_image(R_out, G_out, B_out);
		return 1;
	}

	/* convert to linear RGB: */
	srgb_to_linear(sRGB_in, R_in, G_in, B_in, width, height, p);
	linear_to_srgb(sRGB_out, R_in, G_in, B_in, width, height, p);
	write_bitmap("dfx_aliasing_demo_original.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* tests without low-pass filter: */

	/* subsample without pre-filtering: */
	subsample_image(R_in, G_in, B_in, R_temp, G_temp, B_temp, width, height, p);
	linear_to_srgb(sRGB_out, R_temp, G_temp, B_temp, width, height, p);
	write_bitmap("dfx_aliasing_demo_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply Gaussian filter: */
	filter_image(R_temp, G_temp, B_temp, R_out, G_out, B_out, width, height, p, n, FILT_GAUSS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_subsampled_gauss_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply sinc filter: */
	filter_image(R_temp, G_temp, B_temp, R_out, G_out, B_out, width, height, p, n, FILT_SINC, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_subsampled_sinc_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply Lanczos filter: */
	filter_image(R_temp, G_temp, B_temp, R_out, G_out, B_out, width, height, p, n, FILT_LANCZOS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_subsampled_lanczos_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* tests with low-pass filter: */

	/* using Gaussian filter: */
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_GAUSS, fc);
	subsample_image(R_out, G_out, B_out, R_temp, G_temp, B_temp, width, height, p);
	linear_to_srgb(sRGB_out, R_temp, G_temp, B_temp, width, height, p);
	write_bitmap("dfx_aliasing_demo_gauss_lp_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_GAUSS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_gauss_lp_subsampled_gauss_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply sinc filter: */
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_SINC, fc);
	subsample_image(R_out, G_out, B_out, R_temp, G_temp, B_temp, width, height, p);
	linear_to_srgb(sRGB_out, R_temp, G_temp, B_temp, width, height, p);
	write_bitmap("dfx_aliasing_demo_sinc_lp_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_SINC, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_sinc_lp_subsampled_sinc_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* apply Lanczos filter: */
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_LANCZOS, fc);
	subsample_image(R_out, G_out, B_out, R_temp, G_temp, B_temp, width, height, p);
	linear_to_srgb(sRGB_out, R_temp, G_temp, B_temp, width, height, p);
	write_bitmap("dfx_aliasing_demo_lanczos_lp_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, FILT_LANCZOS, fc);
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap("dfx_aliasing_demo_lanczos_lp_subsampled_lanczos_rec.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R_in, G_in, B_in);
	free_image(R_temp, G_temp, B_temp);
	free_image(R_out, G_out, B_out);
	return 0;
}
