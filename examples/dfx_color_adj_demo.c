/*!
 *  \file   dfx_color_adj_demo.c
 *  \brief  An example demonstrating color adjustment operations. 
 *
 *  This program generates color gradient test patterns and then applies adjustments of hue, saturation, and brightness. 
 *  It can also perform these operations with a user-supplied test file as an input.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.03
 *  \date    November 8, 2026
 */

#include <math.h>
#include <stdio.h>
#include "dfx.h"

 /*!
  *  \brief Generate color test image 1: changing brightness for primary and mixed colors.
  *
  *  \param[in,out]  R,G,B   - image where to place the pattern;
  *  \param[in]      width   - image width in pixels;
  *  \param[in]      height  - image height in pixels;
  *  \param[in]      p       - padding parameter
  *
  *  \return         DFX_SUCCESS if success, or an error code otherwise.
  */
static int generate_color_brightness_gradient_image(float* R, float* G, float* B, int width, int height, int p)
{
	float gamma, dx, L;
	int stripe_height, x, y;

	/* check args: */
	if (R == NULL || G == NULL || B == NULL || width <= 100 || height <= 7*10 || p < 0)
		return DFX_INVARG;

	/* test pattern parameters: */
	gamma = 3.0f;						/* gamma to use to produce color gradients */
	stripe_height = height / 7;			/* height of each stripe */

	/* produce test pattern: */
	for (x = 0; x < width; x++) 
	{
		/* compute luminance: */
		dx = (float)x / (width - 1);	/* delta in [0,1] */
		L = powf(dx, gamma);            /* converted to intended luminance, scaled to [0,1] */
		/* paint red stripe: */
		for (y = 0; y < stripe_height; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = L;
			G[(p + y) * (width + 2 * p) + p + x] = 0.;
			B[(p + y) * (width + 2 * p) + p + x] = 0.;
		}
		/* paint yellow stripe: */
		for (; y < stripe_height * 2; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = L;
			G[(p + y) * (width + 2 * p) + p + x] = L;
			B[(p + y) * (width + 2 * p) + p + x] = 0.;
		}
		/* paint green stripe: */
		for (; y < stripe_height * 3; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = 0.;
			G[(p + y) * (width + 2 * p) + p + x] = L;
			B[(p + y) * (width + 2 * p) + p + x] = 0.;
		}
		/* paint cyan stripe: */
		for (; y < stripe_height * 4; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = 0.;
			G[(p + y) * (width + 2 * p) + p + x] = L;
			B[(p + y) * (width + 2 * p) + p + x] = L;
		}
		/* paint blue stripe: */
		for (; y < stripe_height * 5; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = 0.;
			G[(p + y) * (width + 2 * p) + p + x] = 0.;
			B[(p + y) * (width + 2 * p) + p + x] = L;
		}
		/* paint magenta stripe: */
		for (; y < stripe_height * 6; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = L;
			G[(p + y) * (width + 2 * p) + p + x] = 0.;
			B[(p + y) * (width + 2 * p) + p + x] = L;
		}
		/* paint grayscale stripe: */
		for (; y < height; y++) {
			R[(p + y) * (width + 2 * p) + p + x] = L;
			G[(p + y) * (width + 2 * p) + p + x] = L;
			B[(p + y) * (width + 2 * p) + p + x] = L;
		}
	}

	/* pad image and return: */
	if (p > 0) pad_image(R, G, B, width, height, p, PAD_ZERO);
	return DFX_SUCCESS;
}

/*!
   *  \brief Generate color test image 2: changing saturation for primary and mixed colors.
   *
   *  \param[in,out]  R,G,B   - image where to place the pattern;
   *  \param[in]      width   - image width in pixels;
   *  \param[in]      height  - image height in pixels;
   *  \param[in]      p       - padding parameter
   *
   *  \return         DFX_SUCCESS if success, or an error code otherwise.
   */
static int generate_color_saturation_gradient_image(float* R, float* G, float* B, int width, int height, int p)
{
	float gamma, dy, S, L;
	int stripe_height, x, y;

	/* check args: */
	if (R == NULL || G == NULL || B == NULL || width <= 100 || height <= 7 * 10 || p < 0)
		return DFX_INVARG;

	/* test pattern parameters: */
	gamma = 3.0f;						/* gamma to use to produce brightness gradients */
	stripe_height = height / 6;			/* height of each stripe */

	/* produce test pattern: */
	for (x = 0; x < width; x++)
	{
		/* compute saturation: */
		S = (float)x / (width - 1);											/* S = distance from 0, normalized to [0,1] */
		/* paint red stripe: */
		for (y = 0; y < stripe_height; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);										/* L = function of vertical position within a stripe */
			R[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
		}
		/* paint yellow stripe: */
		for (; y < stripe_height * 2; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);
			R[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
		}
		/* paint green stripe: */
		for (; y < stripe_height * 3; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);
			R[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
		}
		/* paint cyan stripe: */
		for (; y < stripe_height * 4; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);
			R[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
		}
		/* paint blue stripe: */
		for (; y < stripe_height * 5; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);
			R[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
		}
		/* paint magenta stripe: */
		for (; y < height; y++) {
			dy = (float)(y % stripe_height) / (stripe_height - 1);
			L = powf(1.f - dy, gamma);
			R[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
			G[(p + y) * (width + 2 * p) + p + x] = 0 * S + L * (1 - S);
			B[(p + y) * (width + 2 * p) + p + x] = L * S + L * (1 - S);
		}
	}

	/* pad image and return: */
	if (p > 0) pad_image(R, G, B, width, height, p, PAD_ZERO);
	return DFX_SUCCESS;
}

/*!
 *  \brief Report the outcome of a color adjustment.
 *
 *  DFX_INVRANGE is a legitimate outcome: it means that some adjusted colors fell outside of the
 *  sRGB gamut (RGB values outside of [0,1]); such values are clipped by linear_to_srgb().
 */
static void report(int err, const char* what)
{
	if (err == DFX_INVRANGE)
		printf("%s: some colors fell outside of the sRGB gamut (clipped on output)\n", what);
	else if (err != DFX_SUCCESS)
		printf("%s: failed with error code %d\n", what, err);
}

/*!
 *  \brief Apply brightness and saturation adjustments with a given scale factor, store results in bitmap files.
 */
static void test_adjustment(int (*adjust)(const float*, const float*, const float*, float*, float*, float*, int, int, int, float, int),
	const char* name, const float* R, const float* G, const float* B, float* R2, float* G2, float* B2, unsigned char* sRGB,
	int width, int height, int p, int ppm_x, int ppm_y, float scale, const char* scale_name)
{
	char fn[128], what[128];

	/* hard clipping: */
	snprintf(what, sizeof(what), "%s(%g, hard)", name, scale);
	report(adjust(R, G, B, R2, G2, B2, width, height, p, scale, 0), what);
	linear_to_srgb(R2, G2, B2, sRGB, width, height, p);
	snprintf(fn, sizeof(fn), "test_%s_%s_hard.bmp", scale_name, name);
	write_bitmap(fn, sRGB, width, height, ppm_x, ppm_y);

	/* soft transition: */
	snprintf(what, sizeof(what), "%s(%g, soft)", name, scale);
	report(adjust(R, G, B, R2, G2, B2, width, height, p, scale, 2), what);
	linear_to_srgb(R2, G2, B2, sRGB, width, height, p);
	snprintf(fn, sizeof(fn), "test_%s_%s_soft.bmp", scale_name, name);
	write_bitmap(fn, sRGB, width, height, ppm_x, ppm_y);
}

/*!
 *  \brief Color operations / adjustments demo program.
 */
int main(int argc, char* argv[])
{
	/* image parameters (geometry of the generated test patterns; replaced by those of the input file, if given): */
	int width = 1120, height = 630, ppm_x = 2835, ppm_y = 2835; /* 72 dpi */
	int p = 16; /* padding parameter */

	/* working images / color planes: */
	unsigned char *sRGB = NULL, *sRGB_in = NULL;
	float *R = NULL, *G = NULL, *B = NULL;
	float *R2 = NULL, *G2 = NULL, *B2 = NULL;

	/* hue rotation angles to test: */
	static const int degrees[] = {15, 30, 45, 60, 75, 90, 135, 180};
	char fn[128], what[128];
	int i;

	/* 
	 * Read the input file first, since it determines the image geometry:
	 */
	if (argc >= 2) {
		if (read_bitmap(argv[1], &sRGB_in, &width, &height, &ppm_x, &ppm_y) != DFX_SUCCESS) {
			printf("Cannot open input file: %s\n", argv[1]);
			printf("Only 24-bit RGB bitmap files are supported.\n");
			return 1;
		}
	}

	/* allocate memory: */
	if (alloc_srgb_image(&sRGB, width, height) != DFX_SUCCESS
	 || alloc_image(&R, &G, &B, width, height, p) != DFX_SUCCESS
	 || alloc_image(&R2, &G2, &B2, width, height, p) != DFX_SUCCESS) {
		/* print error, free allocated memory & exit: */
		printf("Error: cannot allocate memory for images\n");
		if (sRGB_in != NULL) free_srgb_image(sRGB_in);
		if (sRGB != NULL) free_srgb_image(sRGB);
		if (R != NULL || G != NULL || B != NULL) free_image(R, G, B);
		if (R2 != NULL || G2 != NULL || B2 != NULL) free_image(R2, G2, B2);
		return 1;
	}

	/* 
	 * Load or generate test pattern to work with: 
	 */
	if (sRGB_in != NULL) {
		/* convert the input image to linear & visualize: */
		srgb_to_linear(sRGB_in, R, G, B, width, height, p);
		linear_to_srgb(R, G, B, sRGB, width, height, p);
		write_bitmap("test_image.bmp", sRGB, width, height, ppm_x, ppm_y);
	}
	else {
		/* generate test pattern: */
		generate_color_brightness_gradient_image(R, G, B, width, height, p);
		linear_to_srgb(R, G, B, sRGB, width, height, p);
		write_bitmap("color_gradient_test_image_1.bmp", sRGB, width, height, ppm_x, ppm_y);
	}

	/*
	 * Test brightness adjustments by factors 0.5, 2, and 4, with hard and soft clipping:
	 */
	test_adjustment(adjust_brightness, "brightness", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 0.5f, "05x");
	test_adjustment(adjust_brightness, "brightness", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 2.0f, "2x");
	test_adjustment(adjust_brightness, "brightness", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 4.0f, "4x");

	/*
	 * Test saturation adjustments:
	 */
	if (sRGB_in == NULL) {
		/* generate saturation test pattern; */
		generate_color_saturation_gradient_image(R, G, B, width, height, p);
		linear_to_srgb(R, G, B, sRGB, width, height, p);
		write_bitmap("color_gradient_test_image_2.bmp", sRGB, width, height, ppm_x, ppm_y);
	}
	test_adjustment(adjust_saturation, "saturation", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 0.5f, "05x");
	test_adjustment(adjust_saturation, "saturation", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 2.0f, "2x");
	test_adjustment(adjust_saturation, "saturation", R, G, B, R2, G2, B2, sRGB, width, height, p, ppm_x, ppm_y, 4.0f, "4x");

	/*
	 * Test hue adjustments:
	 */
	for (i = 0; i < (int)(sizeof(degrees) / sizeof(degrees[0])); i++) {
		snprintf(what, sizeof(what), "adjust_hue(+%d degrees)", degrees[i]);
		report(adjust_hue(R, G, B, R2, G2, B2, width, height, p, (float)(M_PI * degrees[i] / 180.0)), what);
		linear_to_srgb(R2, G2, B2, sRGB, width, height, p);
		snprintf(fn, sizeof(fn), "test_hue_plus_%ddeg.bmp", degrees[i]);
		write_bitmap(fn, sRGB, width, height, ppm_x, ppm_y);
	}

	/* free memory and exit: */
	if (sRGB_in) free_srgb_image(sRGB_in);
	free_srgb_image(sRGB);
	free_image(R, G, B);
	free_image(R2, G2, B2);

	return 0;
}

/* dfx_color_adj_demo.c -- end of file */
