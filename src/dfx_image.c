/*!
 *  \file   dfx_image.c
 *  \brief  Basic operations with linear RGB images.
 *
 *  We use planar format for linear RGB images, where each channel is stored in a separate plane.
 *  Each plane is allocated as a single block of memory, containing pixel values in row-major order,
 *  with padding around the boundary: a plane for a width x height image with padding p holds
 *  (width + 2p) x (height + 2p) values, and the image pixel (x,y) is stored at index
 *  (p + y) * (width + 2p) + (p + x).
 *
 *  This is a preferred format for image filtering, blending, and resampling operations.
 *  By Grassmann's laws of color mixing, all such operations on linear RGB images can be performed
 *  independently in each channel.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    November 8, 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dfx.h"

/**************
 *
 *  Allocation/free operations:
 *
 *   plane_size()
 *   alloc_plane()
 *   alloc_image()
 *   free_plane()
 *   free_image()
 */

/*!
 * \brief Compute size of a padded plane of linear image (in pixels):
 */
static size_t plane_size(int width, int height, int p)
{
	return (size_t)(2 * p + height) * (size_t)(2 * p + width);	        /* size of padded plane, in pixels */
}

/*!
 *  \brief Allocate a channel of linear RGB image.
 *
 *  \param[out]     pX     - pointer to a variable to receive the pointer to the allocated plane
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *
 *  \returns        DFX_X error code
 */
int alloc_plane(float** pX, int width, int height, int p)
{
	size_t size;
	if (pX == NULL || height < 0 || width < 0 || p < 0)  return DFX_INVARG;
	size = plane_size(width, height, p);
	if ((*pX = (float*)malloc((size > 0 ? size : 1) * sizeof(float))) == NULL) return DFX_NOMEM;	/* malloc(0) may return NULL */
	return DFX_SUCCESS;
}

/*!
 *  \brief Free image plane.
 */
int free_plane(float* X)
{
	if (X != NULL) free(X);
	return DFX_SUCCESS;
}

/*!
 *  \brief Allocate linear RGB image.
 *
 *  \param[out]     pR, pG, pB  - pointers to variables to receive the pointers to allocated R,G,B planes
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *
 *  \returns        DFX_X error code
 */
int alloc_image(float** pR, float** pG, float** pB, int width, int height, int p)
{
	float* R = NULL, * G = NULL, * B = NULL;

	/* check arguments */
	if (pR == NULL || pG == NULL || pB == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* allocate planes: */
	if (alloc_plane(&R, width, height, p) != DFX_SUCCESS
		|| alloc_plane(&G, width, height, p) != DFX_SUCCESS
		|| alloc_plane(&B, width, height, p) != DFX_SUCCESS) {
		/* free memory & exit: */
		free_image(R, G, B);
		return DFX_NOMEM;
	}

	/* pass pointers & exit: */
	*pR = R; *pG = G; *pB = B;
	return DFX_SUCCESS;
}

/*!
 *  \brief Free planes of linear RGB image.
 */
int free_image(float* R, float* G, float* B)
{
	free_plane(R);
	free_plane(G);
	free_plane(B);
	return DFX_SUCCESS;
}

/**************
 *
 *  Zero & unit initializations:
 *
 *   zero_plane()
 *   zero_image()
 *
 *   unit_plane()
 *   unit_image()
 */

 /*!
  *  \brief Fill channel of linear RGB image with zeros.
  */
int zero_plane(float* X, int width, int height, int p)
{
	if (X == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	/* set everything to 0s */
	memset((void*)X, 0, plane_size(width, height, p) * sizeof(float));
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill linear RGB image with zeros.
 */
int zero_image(float* R, float* G, float* B, int width, int height, int p)
{
	int err;
	if (R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if ((err = zero_plane(R, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = zero_plane(G, width, height, p)) != DFX_SUCCESS) return err;
	return zero_plane(B, width, height, p);
}

/*!
 *  \brief Fill channel of linear RGB image with ones.
 */
int unit_plane(float* X, int width, int height, int p)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;

	/* check parameters: */
	if (X == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* set everything to 1s: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++)
		X[y * width_p + x] = 1.0f;

	return DFX_SUCCESS;
}

/*!
 *  \brief Fill linear RGB image with ones.
 */
int unit_image(float* R, float* G, float* B, int width, int height, int p)
{
	int err;

	/* check parameters: */
	if (R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* initialize planes: */
	if ((err = unit_plane(R, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = unit_plane(G, width, height, p)) != DFX_SUCCESS) return err;
	return unit_plane(B, width, height, p);
}

/**************
 *
 *  Copy, scale, add, subtract, blend, and clip images:
 *
 *   copy_plane()
 *   copy_image()
 *   scale_plane()
 *   scale_image()
 *   add_planes()
 *   add_images()
 *   subtract_planes()
 *   subtract_images()
 *   blend_planes()
 *   blend_images()
 *   clip_plane()
 *   clip_image()
 *
 *  All these operations process the padded area as well, so the padding of their results stays valid.
 */

/*!
 *  \brief Copy channel/plane of a linear RGB image.
 */
int copy_plane(const float* X_in, float* X_out, int width, int height, int p)
{
	/* check parameters: */
	if (X_in == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* copy plane (memmove allows the two planes to overlap): */
	if (X_in != X_out)
		memmove((void*)X_out, (const void*)X_in, plane_size(width, height, p) * sizeof(float));
	return DFX_SUCCESS;
}

/*!
 *  \brief Copy linear RGB image.
 */
int copy_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p)
{
	int err;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* copy planes: */
	if ((err = copy_plane(R_in, R_out, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = copy_plane(G_in, G_out, width, height, p)) != DFX_SUCCESS) return err;
	return copy_plane(B_in, B_out, width, height, p);
}

/*!
 *  \brief Scale channel of linear RGB image by a constant.
 */
int scale_plane(const float* X_in, float* X_out, int width, int height, int p, float scale)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;

	/* check parameters: */
	if (X_in == NULL || X_out == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* scale image: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++)
		X_out[y * width_p + x] = X_in[y * width_p + x] * scale;

	return DFX_SUCCESS;
}

/*!
 *  \brief Scale linear RGB image by a constant.
 */
int scale_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float scale)
{
	int err;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* scale planes: */
	if ((err = scale_plane(R_in, R_out, width, height, p, scale)) != DFX_SUCCESS) return err;
	if ((err = scale_plane(G_in, G_out, width, height, p, scale)) != DFX_SUCCESS) return err;
	return scale_plane(B_in, B_out, width, height, p, scale);
}

/*!
 *  \brief Add two linear RGB image planes.
 */
int add_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;

	/* check parameters: */
	if (X_1 == NULL || X_2 == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* add two images: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++)
		X_out[y * width_p + x] = X_1[y * width_p + x] + X_2[y * width_p + x];

	return DFX_SUCCESS;
}

/*!
 *  \brief Add two linear RGB images.
 */
int add_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p)
{
	int err;

	/* check parameters: */
	if (R_1 == NULL || G_1 == NULL || B_1 == NULL) return DFX_INVARG;
	if (R_2 == NULL || G_2 == NULL || B_2 == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* add planes: */
	if ((err = add_planes(R_1, R_2, R_out, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = add_planes(G_1, G_2, G_out, width, height, p)) != DFX_SUCCESS) return err;
	return add_planes(B_1, B_2, B_out, width, height, p);
}

/*!
 *  \brief Subtract two linear RGB image planes.
 */
int subtract_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;

	/* check parameters: */
	if (X_1 == NULL || X_2 == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* subtract two images: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++)
		X_out[y * width_p + x] = X_1[y * width_p + x] - X_2[y * width_p + x];

	return DFX_SUCCESS;
}

/*!
 *  \brief Subtract two linear RGB images.
 */
int subtract_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p)
{
	int err;

	/* check parameters: */
	if (R_1 == NULL || G_1 == NULL || B_1 == NULL) return DFX_INVARG;
	if (R_2 == NULL || G_2 == NULL || B_2 == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* subtract planes: */
	if ((err = subtract_planes(R_1, R_2, R_out, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = subtract_planes(G_1, G_2, G_out, width, height, p)) != DFX_SUCCESS) return err;
	return subtract_planes(B_1, B_2, B_out, width, height, p);
}

/*!
 *  \brief Blend two linear RGB image planes: X_out = alpha * X_1 + (1 - alpha) * X_2, with alpha in [0,1].
 */
int blend_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p, float alpha)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;

	/* check parameters: */
	if (X_1 == NULL || X_2 == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (alpha < 0.f || alpha > 1.0f) return DFX_INVARG;

	/* blend two images: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++)
		X_out[y * width_p + x] = alpha * X_1[y * width_p + x] + (1.0f - alpha) * X_2[y * width_p + x];

	return DFX_SUCCESS;
}

/*!
 *  \brief Blend two linear RGB images.
 */
int blend_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p, float alpha)
{
	int err;

	/* check parameters: */
	if (R_1 == NULL || G_1 == NULL || B_1 == NULL) return DFX_INVARG;
	if (R_2 == NULL || G_2 == NULL || B_2 == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (alpha < 0.f || alpha > 1.0f) return DFX_INVARG;

	/* blend planes: */
	if ((err = blend_planes(R_1, R_2, R_out, width, height, p, alpha)) != DFX_SUCCESS) return err;
	if ((err = blend_planes(G_1, G_2, G_out, width, height, p, alpha)) != DFX_SUCCESS) return err;
	return blend_planes(B_1, B_2, B_out, width, height, p, alpha);
}

/*!
 *  \brief Clip values of a linear RGB image plane to the range [lo, hi].
 *
 *  Filtering with kernels having negative taps (sinc, Lanczos) can produce values slightly
 *  outside of the [0,1] range of linear RGB. This function brings them back into range,
 *  which is required before conversions to XYZ and other color spaces.
 */
int clip_plane(const float* X_in, float* X_out, int width, int height, int p, float lo, float hi)
{
	int width_p = width + 2 * p;	/* width of padded plane */
	int height_p = height + 2 * p;	/* height of padded plane */
	int x, y;
	float v;

	/* check parameters: */
	if (X_in == NULL || X_out == NULL || height < 0 || width < 0 || p < 0 || lo > hi) return DFX_INVARG;

	/* clip: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++) {
		v = X_in[y * width_p + x];
		X_out[y * width_p + x] = (v < lo) ? lo : (v > hi) ? hi : v;
	}

	return DFX_SUCCESS;
}

/*!
 *  \brief Clip values of a linear RGB image to the range [lo, hi].
 */
int clip_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float lo, float hi)
{
	int err;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0 || lo > hi) return DFX_INVARG;

	/* clip planes: */
	if ((err = clip_plane(R_in, R_out, width, height, p, lo, hi)) != DFX_SUCCESS) return err;
	if ((err = clip_plane(G_in, G_out, width, height, p, lo, hi)) != DFX_SUCCESS) return err;
	return clip_plane(B_in, B_out, width, height, p, lo, hi);
}

/*******************
 *
 *  Image padding functions
 *
 *   pad_plane()
 *   pad_image()
 */

/*!
 *  \brief Fill the padding border of a given image plane/channel.
 *
 *  Three rules are supported:
 *    PAD_ZERO       - border samples are set to 0;
 *    PAD_REPLICATE  - border samples repeat the nearest boundary pixel;
 *    PAD_REFLECT    - the image is mirrored about its boundary, repeating the boundary pixel
 *                     (..., x[2], x[1], x[0] | x[0], x[1], x[2], ...); requires p <= width and p <= height.
 *
 *  \param[in,out]  X      - pointer to a plane to process
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *  \param[in]      t      - padding type (see PAD_X constants)
 *
 *  \returns        DFX_X error code
 */
int pad_plane(float* X, int width, int height, int p, int t)
{
	int w_lin = width + 2 * p;			    /* w_lin = width of padded linear RGB image */
	int x, y;

	/* check parameters: */
	if (X == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (t < 0 || t >= N_PADDINGS) return DFX_INVARG;
	if (t == PAD_REFLECT && (p > width || p > height)) return DFX_INVARG;

	/* nothing to do if there is no padding: */
	if (p == 0) return DFX_SUCCESS;

	/* check padding type: */
	if (t == PAD_ZERO)
	{
		/* zero-pad left and right edges: */
		for (y = 0; y < height; y++) for (x = 0; x < p; x++) {
			X[(p + y) * w_lin + x] = 0;
			X[(p + y) * w_lin + p + width + x] = 0;
		}
		/* zero-pad top and bottom edges: */
		for (y = 0; y < p; y++) for (x = 0; x < w_lin; x++) {
			X[y * w_lin + x] = 0;
			X[(p + height + y) * w_lin + x] = 0;
		}
	}
	else if (t == PAD_REPLICATE)
	{
		/* replicate-pad left and right edges: */
		for (y = 0; y < height; y++) for (x = 0; x < p; x++) {
			X[(p + y) * w_lin + x] = X[(p + y) * w_lin + p];
			X[(p + y) * w_lin + p + width + x] = X[(p + y) * w_lin + p + width - 1];
		}
		/* replicate-pad top and bottom edges: */
		for (y = 0; y < p; y++) for (x = 0; x < w_lin; x++) {
			X[y * w_lin + x] = X[p * w_lin + x];
			X[(p + height + y) * w_lin + x] = X[(p + height - 1) * w_lin + x];
		}
	}
	else /* PAD_REFLECT */
	{
		/* reflect-pad left and right edges: */
		for (y = 0; y < height; y++) for (x = 0; x < p; x++) {
			X[(p + y) * w_lin + p - 1 - x] = X[(p + y) * w_lin + p + x];                  /* flip horizontal */
			X[(p + y) * w_lin + p + width + x] = X[(p + y) * w_lin + p + width - 1 - x];
		}
		/* reflect-pad top and bottom edges: */
		for (y = 0; y < p; y++) for (x = 0; x < w_lin; x++) {
			X[(p - 1 - y) * w_lin + x] = X[(p + y) * w_lin + x];                          /* flip vertical */
			X[(p + height + y) * w_lin + x] = X[(p + height - 1 - y) * w_lin + x];
		}
	}

	return DFX_SUCCESS;
}

/*!
 *  \brief Fill the padding borders of an image.
 */
int pad_image(float* R, float* G, float* B, int width, int height, int p, int t)
{
	int err;
	if (R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if ((err = pad_plane(R, width, height, p, t)) != DFX_SUCCESS) return err;
	if ((err = pad_plane(G, width, height, p, t)) != DFX_SUCCESS) return err;
	return pad_plane(B, width, height, p, t);
}

/* dfx_image.c -- end of file */
