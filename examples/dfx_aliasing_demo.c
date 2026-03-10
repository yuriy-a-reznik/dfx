/*!
 *  \file   dfx_aliasing_demo.c
 *  \brief  An example program illustrating the effects of aliasing with different filters.
 *
 *  An example program illustrating the effects of aliasing when the subsampling is done without prior application of a low-pass filter.
 *  For most obvious visual results please use dfx-zone-plate.bmp image. 
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.01
 *  \date    March 7, 2026
 */

#include <stdio.h>
#include "dfx.h"

 /**************
  *
  *  Subsamling operations.
  *
  *  These functions have the effect of subsampling an a plane (or an image) by a factor of 2, 
  *  and then zero-fill upsampling it to the original resolution.
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
 *  \brief Apply reconstruction filter, produce visualizations, compute SNR/MSE, and print results:
 * 
 *  \param[in] n      - filter order
 *  \param[in] filter - reconstruciton filter type (FILTER_SINC or FILTER_LANCZOS)
 *  \param[in] fc     - cutoff
 */
static int test_reconstruction(float* R_ref, float* G_ref, float* B_ref, float* R_in, float *G_in, float *B_in, float* R_out, float* G_out, float* B_out, unsigned char *sRGB_out,
	int width, int height, int p, int n, int filter, float fc, char *bmp_fn, int ppm_x, int ppm_y)
{
	float snr, mse;

	/* check parameters */
	if (R_ref == NULL || G_ref == NULL || B_ref == NULL) return DFX_INVARG;
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL || sRGB_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (bmp_fn == NULL) return DFX_INVARG;

	/* reconstruct image: */
	filter_image(R_in, G_in, B_in, R_out, G_out, B_out, width, height, p, n, filter, fc);
	scale_image(R_out, G_out, B_out, R_out, G_out, B_out, width, height, p, 4.0);

	/* produce bitmap file: */
	linear_to_srgb(sRGB_out, R_out, G_out, B_out, width, height, p);
	write_bitmap(bmp_fn, sRGB_out, width, height, ppm_x, ppm_y);

	/* compute MSE, SNR: */
	snr_image(R_ref, G_ref, B_ref, R_out, G_out, B_out, width, height, p, &mse, &snr);
	printf("original vs %s: mse=%g, snr=%g[db]\n", bmp_fn, mse, snr);

	/* success: */
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
	float* R_lp = NULL, * G_lp = NULL, * B_lp = NULL;
	float *R_ss = NULL, * G_ss = NULL, * B_ss = NULL;
	float* R_out = NULL, * G_out = NULL, * B_out = NULL;

	/* padding & filter parameters: */
	int p = 50;
	float fc = 0.25;								
	int n = 3;

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
		printf("Only 24-bit RGB bitmap files are supported.\n");
		return 1;
	}

	/* allocate working images and planes: */
	if (alloc_srgb_image(&sRGB_out, width, height) != DFX_SUCCESS
		|| alloc_image(&R_in, &G_in, &B_in, width, height, p) != DFX_SUCCESS
		|| alloc_image(&R_lp, &G_lp, &B_lp, width, height, p) != DFX_SUCCESS
		|| alloc_image(&R_ss, &G_ss, &B_ss, width, height, p) != DFX_SUCCESS
		|| alloc_image(&R_out, &G_out, &B_out, width, height, p) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory for images\n");
		/* free memory and exit: */
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB_out != NULL) free_srgb_image(sRGB_out);
		if (R_in != NULL && G_in != NULL && B_in != NULL) free_image(R_in, G_in, B_in);
		if (R_lp != NULL && G_lp != NULL && B_lp != NULL) free_image(R_lp, G_lp, B_lp);
		if (R_ss != NULL && G_ss != NULL && B_ss != NULL) free_image(R_ss, G_ss, B_ss);
		if (R_out != NULL && G_out != NULL && B_out != NULL) free_image(R_out, G_out, B_out);
		return 1;
	}

	/* convert to linear: */
	srgb_to_linear(sRGB_in, R_in, G_in, B_in, width, height, p);
	linear_to_srgb(sRGB_out, R_in, G_in, B_in, width, height, p);
	write_bitmap("original.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* produce a low-pass version (20-th order sinc filter): */
	filter_image(R_in, G_in, B_in, R_lp, G_lp, B_lp, width, height, p, 20, FILT_SINC, fc);
	pad_image(R_lp, G_lp, B_lp, width, height, p, PAD_REFLECT);
	linear_to_srgb(sRGB_out, R_lp, G_lp, B_lp, width, height, p);
	write_bitmap("lowpass.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* subsample original image: */
	subsample_image(R_in, G_in, B_in, R_ss, G_ss, B_ss, width, height, p);
	linear_to_srgb(sRGB_out, R_ss, G_ss, B_ss, width, height, p);
	write_bitmap("original_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* run reconstruction tests: */
	test_reconstruction(R_in, G_in, B_in, R_ss, G_ss, B_ss, R_out, G_out, B_out, sRGB_out, width, height, p, n, FILT_SINC, fc, "sync_reconstructed.bmp", ppm_x, ppm_y);
	test_reconstruction(R_in, G_in, B_in, R_ss, G_ss, B_ss, R_out, G_out, B_out, sRGB_out, width, height, p, n, FILT_LANCZOS, fc, "lanczos_reconstructed.bmp", ppm_x, ppm_y);

	/* subsample lowpass image: */
	subsample_image(R_lp, G_lp, B_lp, R_ss, G_ss, B_ss, width, height, p);
	linear_to_srgb(sRGB_out, R_ss, G_ss, B_ss, width, height, p);
	write_bitmap("lowpass_subsampled.bmp", sRGB_out, width, height, ppm_x, ppm_y);

	/* run reconstruction tests: */
	test_reconstruction(R_in, G_in, B_in, R_ss, G_ss, B_ss, R_out, G_out, B_out, sRGB_out, width, height, p, n, FILT_SINC, fc, "lowpass_sync_reconstructed.bmp", ppm_x, ppm_y);
	test_reconstruction(R_in, G_in, B_in, R_ss, G_ss, B_ss, R_out, G_out, B_out, sRGB_out, width, height, p, n, FILT_LANCZOS, fc, "lowpass_lanczos_reconstructed.bmp", ppm_x, ppm_y);

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R_in, G_in, B_in);
	free_image(R_lp, G_lp, B_lp);
	free_image(R_ss, G_ss, B_ss);
	free_image(R_out, G_out, B_out);
	return 0;
}
