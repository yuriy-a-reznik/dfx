/*!
 *  \file   dfx_resampler_demo.c
 *  \brief  An example demonstrating how to use image resampling functions.
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

 /*!
  *  \brief Test resampler operation.
  *
  *  \param[in] fn_in      - name of input bitmap file to use for scaling
  *  \param[in] fn_out     - name of output bitmap to create to store the results
  *  \param[in] M,N        - rational scaling factor (M/N) to apply 
  *  \param[in] n          - filter order to use
  */
static void test_resampler(char* fn_in, char* fn_out, int M, int N, int n)
{
	/* image parameters: */
	int width_in = 0, height_in = 0, ppm_x = 0, ppm_y = 0;
	int width_out = 0, height_out = 0;

	/* internal buffers: */
	unsigned char* sRGB_in = NULL, * sRGB_out = NULL;
	float* R_in = NULL, * G_in = NULL, * B_in = NULL;
	float* R_out = NULL, * G_out = NULL, * B_out = NULL;

	/* padding: */
	int p = 10;
	if (p < n) p = n;

	/* print test type: */
	printf("Testing conversion: %s -> %s, scale factor: %d/%d, filter order: %d\n", fn_in, fn_out, M, N, n);

	/* read input file: */
	if (read_bitmap(fn_in, &sRGB_in, &width_in, &height_in, &ppm_x, &ppm_y) != DFX_SUCCESS) {
		printf("Cannot open input file: %s\n", fn_in);
		printf("Only 24-bit RGB bitmap files are supported.\n");
		return;
	}

	/* check if input dimensions are meaninfull: */
	if (width_in < 64 || height_in < 48) {
		printf("Input image is too small (%dx%d).", width_in, height_in);
		return;
	}

	/* check if scaling factor makes any sense: */
	if (N < 1 || M < 1) {
		printf("Scaling fraction parameters (%d/%d) must be positive.", M, N);
		return;
	}

	/* check if scaling by requested factor is possible: */
	if (width_in % N != 0 || height_in % N != 0) {
		/* skip test if it is not possible: */
		printf("Conversion of %dx%d image by a factor %d/%d is not possible. Input width/height must be multiple of %d\n", width_in, height_in, M, N, N);
		free_srgb_image(sRGB_in);
		return;
	}

	/* compute scaled dimensions: */
	width_out = (width_in * M) / N;
	height_out = (height_in * M) / N;

	/* allocate working images and planes: */
	if (alloc_image(&R_in, &G_in, &B_in, width_in, height_in, p) == DFX_SUCCESS
		&& alloc_image(&R_out, &G_out, &B_out, width_out, height_out, p) == DFX_SUCCESS
		&& alloc_srgb_image(&sRGB_out, width_out, height_out) == DFX_SUCCESS)
	{
		/* perform conversions: */

		/* input RGB -> linear: */
		srgb_to_linear(sRGB_in, R_in, G_in, B_in, width_in, height_in, p);

		/* apply resampling: */
		resample_image(R_in, G_in, B_in, R_out, G_out, B_out, width_in, height_in, width_out, height_out, p, n);

		/* linear to sRGB: */
		linear_to_srgb(sRGB_out, R_out, G_out, B_out, width_out, height_out, p);

		/* write output bitmap file: */
		write_bitmap(fn_out, sRGB_out, width_out, height_out, ppm_x, ppm_y);
	}
	else
	{
		/* report memory issue: */
		printf("Cannot allocate memory\n");
	}

	/* free memory and exit: */
	free_srgb_image(sRGB_in);
	free_srgb_image(sRGB_out);
	free_image(R_in, G_in, B_in);
	free_image(R_out, G_out, B_out);
	return;
}

/*!
 *  \brief An example program, showing how to use image resamping functions.
 */
int main(int argc, char* argv[])
{
	/* check program parameters: */
	if (argc < 2) {
		printf("Usage: dfx_resample_demo <bitmap image>\n");
		return 1;
	}

	/* test resamping by filter with n=3: */
	test_resampler(argv[1], "dfx_double_3.bmp",          2, 1, 3);
	test_resampler(argv[1], "dfx_triple_3.bmp",          3, 1, 3);
	test_resampler(argv[1], "dfx_quadruple_3.bmp",       4, 1, 3);
	test_resampler(argv[1], "dfx_half_3.bmp",            1, 2, 3);
	test_resampler(argv[1], "dfx_third_3.bmp",           1, 3, 3);
	test_resampler(argv[1], "dfx_quarter_3.bmp",         1, 4, 3);
	test_resampler(argv[1], "dfx_two_thirds_3.bmp",      2, 3, 3);
	test_resampler(argv[1], "dfx_four_thirds_3.bmp",     4, 3, 3);
	test_resampler(argv[1], "dfx_three_quarters_3.bmp",  3, 4, 3);

	/* test resamping by filter with n=5: */
	test_resampler(argv[1], "dfx_double_5.bmp",          2, 1, 5);
	test_resampler(argv[1], "dfx_triple_5.bmp",          3, 1, 5);
	test_resampler(argv[1], "dfx_quadruple_5.bmp",       4, 1, 5);
	test_resampler(argv[1], "dfx_half_5.bmp",            1, 2, 5);
	test_resampler(argv[1], "dfx_third_5.bmp",           1, 3, 5);
	test_resampler(argv[1], "dfx_quarter_5.bmp",         1, 4, 5);
	test_resampler(argv[1], "dfx_two_thirds_5.bmp",      2, 3, 5);
	test_resampler(argv[1], "dfx_four_thirds_5.bmp",     4, 3, 5);
	test_resampler(argv[1], "dfx_three_quarters_5.bmp",  3, 4, 5);

	/* test resamping by filter with n=7: */
	test_resampler(argv[1], "dfx_double_7.bmp",          2, 1, 7);
	test_resampler(argv[1], "dfx_triple_7.bmp",          3, 1, 7);
	test_resampler(argv[1], "dfx_quadruple_7.bmp",       4, 1, 7);
	test_resampler(argv[1], "dfx_half_7.bmp",            1, 2, 7);
	test_resampler(argv[1], "dfx_third_7.bmp",           1, 3, 7);
	test_resampler(argv[1], "dfx_quarter_7.bmp",         1, 4, 7);
	test_resampler(argv[1], "dfx_two_thirds_7.bmp",      2, 3, 7);
	test_resampler(argv[1], "dfx_four_thirds_7.bmp",     4, 3, 7);
	test_resampler(argv[1], "dfx_three_quarters_7.bmp",  3, 4, 7);

	/* all done */
	return 0;
}

/* dfx_reampler_demo.c - end of file */
