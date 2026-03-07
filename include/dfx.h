/*!
 *  \file   dfx.h
 *  \brief  DFX: Digital image processing Functions and eXamples.
 *
 *  This library provides a collection of basic image‑processing primitives
 *  implemented entirely in C, without relying on external libraries.
 *
 *  •  read_bitmap() and write_bitmap() load and store 24‑bit packed sRGB images.
 *     In memory, these images are represented as arrays of unsigned chars.
 *     The parameters “width” and “height” specify the image dimensions.
 *
 *  •  srgb_to_linear() and linear_to_srgb() convert between sRGB and linear RGB.
 *     In linear RGB form, images are stored as three float arrays: R (red),
 *     G (green), and B (blue). Linear representation is fundamental. All color
 *     mixing and image‑processing operations in this library are implemented 
 *     in linear space.
 *
 *  •  When converting from sRGB to linear RGB, images may be padded.
 *     The parameter “p” defines the number of pixels reserved around the
 *     boundary of the original image. Padding is useful for filtering and
 *     resampling operations.
 *
 *  •  linear_to_luminance() and luminance_to_grayscale_image() extract image
 *     luminance (the Y channel in CIE 1931 XYZ space) and convert luminance
 *     channel into a grayscale RGB image.
 *
 *  •  filter_plane() and filter_image() implement a variety of image‑filtering
 *     operations.
 *
 *  •  resample_plane() and resample_image() implement image‑resampling
 *     operations.
 *
 *  •  dft_plane(), dft_magnitude(), and dft_phase() compute the Discrete
 *     Fourier Transform of an image plane and extract magnitude and phase.
 *
 *  Examples demonstrating how to use these operations can be found in
 *  the dfx/examples directory.
 *  
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 * 
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#ifndef _DFX_H_
#define _DFX_H_  1

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Filter kenels: */
enum {
	FILT_SINC = 0,		/* truncated sinc filter */
	FILT_LANCZOS = 1,	/* Lanczos filter */
	FILT_GAUSS = 2,		/* Gaussian filter */
	/* ... */
	N_FILTERS			/* number of supported filters */
};

/* Image padding operations: */
enum {
	PAD_ZERO = 0,		/* padding by zeros */
	PAD_REPLICATE = 1,	/* padding by replicating boundary pixels */
	PAD_REFLECT = 2		/* padding by reflecting image along boundary */
};

/* Error codes returned by DFX library functions: */
enum {
	DFX_SUCCESS = 0,	/* success */
	DFX_INVARG = 1,		/* invalid arguments */
	DFX_CANTOPEN = 2,	/* file cannot be opened/created */
	DFX_IOERR = 3,		/* file reading/writing error */
	DFX_NOTSUP = 4,		/* not supported input format */
	DFX_NOMEM = 5		/* out of memory */
};

/* Reading/writing bitmap files (dfx_bmp.c): */
extern int read_bitmap(char *filename, unsigned char **p_srgb, int* width, int* height, int* ppm_x, int* ppm_y);
extern int write_bitmap(char* filename, unsigned char* srgb, int width, int height, int ppm_x, int ppm_y);

/* Allocation/deallocation of sRGB images (dfx_bmp.c): */
extern int alloc_srgb_image(unsigned char** p_srgb, int width, int height);
extern int free_srgb_image(unsigned char* srgb);

/* Alloc/dealloc of linear images (dfx_image.c): */
extern int alloc_plane(float** pX, int wip_xdth, int height, int p);
extern int free_plane(float* X);
extern int alloc_image(float** pR, float** pG, float** pB, int width, int height, int p);
extern int free_image(float* R, float* G, float* B);

/* Zero & copy functions: */
extern int zero_plane(float* X, int width, int height, int p);
extern int zero_image(float* R, float* G, float* B, int width, int height, int p);
extern int copy_plane(float* X_in, float* X_out, int width, int height, int p);
extern int copy_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p);

/* Padding functions: */
extern int pad_plane(float* X, int width, int height, int p, int padding_type);
extern int pad_image(float* R, float* G, float* B, int width, int height, int p, int padding_type);

/* Conversions between sRGB and linear RGB: */
extern int srgb_to_linear(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p);
extern int linear_to_srgb(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p);
extern int linear_to_srgb_dithered(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p);

/* Luminance extraction and visualization: */
extern int linear_to_luminance(float* Y, float* R, float* G, float* B, int width, int height, int p);
extern int luminance_to_grayscale_image(float* Y, float* R, float* G, float* B, int width, int height, int p);

/* Filtering operations (dfx_filter.c): */
extern int filter_plane(float* X, float* X_out, int width, int height, int p, int n, int t, float fc);
extern int filter_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, int n, int t, float fc);

/* Resampling operations: */
extern int resample_plane(float* X, float* X_out, int width_in, int height_in, int width_out, int height_out, int p, int n);
extern int resample_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width_in, int height_in, int width_out, int height_out, int p, int n);

/* DFT-related funcitons (dfx_dft.c): */
extern int dft_plane(float* x, float* reX, float* imX, int width, int height, int p, int center_dc);
extern int dft_magnitude(float* reX, float* imX, float* magX, int width, int height, int p);
extern int dft_phase(float* reX, float* imX, float* phaseX, int width, int height, int p);

#ifdef __cplusplus
}
#endif

#endif // _DFX_H
