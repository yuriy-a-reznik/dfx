/*!
 *  \file   dfx_srgb.c
 *  \brief  DFX library: conversions between sRGB and linear RGB representations.
 *
 *  The sRGB transfer function (IEC 61966-2-1) maps linear light x in [0,1] to the
 *  companded value y in [0,1]:
 *
 *      y = 12.92 x                      for x <= 0.0031308,
 *      y = 1.055 x^(1/2.4) - 0.055      for x >  0.0031308,
 *
 *  and the 8-bit code value is round(255 y). This module implements this mapping,
 *  its inverse, and an error-diffusion (dithered) variant of the forward mapping.
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
 *  sRGB-related operations:
 *
 *   quant_8bit()
 *   rec_8bit()
 *   to_srgb()
 *   to_linear()
 */

/*!
 *  \brief 8-bit uniform quantizer function (rounding to the nearest code, with clipping to [0,255]):
 */
static unsigned char quant_8bit(float x)
{
	int y = (int)floor(x * 255 + 0.5);
	if (y < 0) y = 0; else if (y > 255) y = 255;
	return (unsigned char)y;
}

/*!
 *  \brief 8-bit uniform quantizer - reconstruction function:
 */
static float rec_8bit(unsigned char y)
{
	float x = (float)y / 255.f;
	return x;
}

/*!
 *  \brief Linear to sRGB gamma space conversion:
 */
static float to_srgb(float x)
{
	/* sRGB gamma: */
	float y = (float)((x > 0.0031308) ? 1.055 * pow(x, 1 / 2.4) - 0.055 : 12.92 * x);
	return y;
}

/*!
 *  \brief sRGB gamma to linear space conversion:
 */
static float to_linear(float y)
{
	/* inverse sRGB gamma (the threshold 0.04045 = 12.92 * 0.0031308 matches the forward function): */
	float x = (float)((y > 0.0031308 * 12.92) ? pow((y + 0.055) / 1.055, 2.4) : y / 12.92);
	return x;
}

/**************
 *
 *  Frame-level sRGB/linear conversion functions:
 *
 *   srgb_to_linear()
 *   linear_to_srgb()
 *   linear_to_srgb_dithered()
 */

 /*!
  *  \brief Converts an 8-bit-per channel sRGB image to padded linear-space RGB image.
  *
  *  Inverse sRGB gamma is applied. RGB primary colors stay the same as in sRGB.
  *  The sRGB image is assumed to be stored bottom-up (as per Microsoft's conventions in bitmap files);
  *  the linear planes are stored top-down.
  *  Output linear planes are padded (PAD_REPLICATE) to support subsequent filtering operations.
  *
  *  \param[in]  sRGB     - pointer to an sRGB image
  *  \param[out] R, G, B  - pointers to linear R, G, and B planes
  *  \param[in]  width    - image width
  *  \param[in]  height   - image height
  *  \param[in]  p        - padding parameter
  *
  *  \returns    DFX_X error code
  */
int srgb_to_linear(const unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);	/* w_srgb = offset to the next line in srgb bitmap image */
	int w_lin = width + 2 * p;			    /* w_lin = width of padded linear RGB image */
	int x, y;

	/* check parameters: */
	if (sRGB == NULL || R == NULL || G == NULL || B == NULL
	 || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* extract R,G,B channels: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++)
	{
		B[(p + y) * w_lin + p + x] = to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 0]));
		G[(p + y) * w_lin + p + x] = to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 1]));
		R[(p + y) * w_lin + p + x] = to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 2]));
	}

	/* add padding: */
	return pad_image(R, G, B, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation.
 *
 *  sRGB gamma is applied, and the resulting values quantized to 8-bit outputs (values outside [0,1] are clipped).
 *  Padding is removed. The output consists of packed 24-bit sRGB values, suitable for writing in a bitmap file.
 *
 *  \param[in]  R, G, B  - pointers to linear R, G, and B planes
 *  \param[out] sRGB     - pointer to an sRGB image
 *  \param[in]  width    - image width
 *  \param[in]  height   - image height
 *  \param[in]  p        - padding parameter
 *
 *  \returns    DFX_X error code
 */
int linear_to_srgb(const float* R, const float* G, const float* B, unsigned char* sRGB, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);		/* w_srgb = offset to the next line in bitmap image */
	int w_lin = 2 * p + width;					/* w_lin = width of padded linear RGB image */
	int x, y;

	/* check parameters: */
	if (sRGB == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* produce sRGB pixel values: */
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++)
		{
			sRGB[(height-1-y) * w_srgb + x * 3 + 0] = quant_8bit(to_srgb(B[(p + y) * w_lin + p + x]));
			sRGB[(height-1-y) * w_srgb + x * 3 + 1] = quant_8bit(to_srgb(G[(p + y) * w_lin + p + x]));
			sRGB[(height-1-y) * w_srgb + x * 3 + 2] = quant_8bit(to_srgb(R[(p + y) * w_lin + p + x]));
		}
		/* clear the row padding bytes, so that the buffer is fully initialized: */
		for (x = width * 3; x < w_srgb; x++)
			sRGB[(height-1-y) * w_srgb + x] = 0;
	}

	return DFX_SUCCESS;
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation with error diffusion.
 *
 *  sRGB gamma is applied, and the resulting values are quantized to 8-bit outputs.
 *  Banding is minimized by Floyd-Steinberg diffusion of the quantization errors: the error made at
 *  pixel (x,y), measured in linear space, is added to the not-yet-quantized neighbors with weights
 *
 *                 (x,y)   7/16
 *      3/16       5/16    1/16      (next row)
 *
 *  Terms that would fall outside of the image are dropped. The input planes are not modified:
 *  the accumulated errors are kept in two row buffers (the current and the next row).
 *  Padding is removed. The output consists of packed 24-bit sRGB values, suitable for writing in a bitmap file.
 *
 *  \param[in]  R, G, B  - pointers to linear R, G, and B planes
 *  \param[out] sRGB     - pointer to an sRGB image
 *  \param[in]  width    - image width
 *  \param[in]  height   - image height
 *  \param[in]  p        - padding parameter
 *
 *  \returns    DFX_X error code
 */
int linear_to_srgb_dithered(const float* R, const float* G, const float* B, unsigned char* sRGB, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);	/* w_srgb = offset to the next line in bitmap image */
	int w_lin = 2 * p + width;			    /* w_lin = width of padded linear RGB image */
	float* err = NULL;						/* accumulated errors: 2 rows x width pixels x 3 channels */
	float* err_cur, * err_next, * swap;
	const float* in[3];
	float v, d;
	int x, y, c;

	/* check parameters: */
	if (sRGB == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* allocate and clear error buffers: */
	if ((err = (float*)calloc((size_t)width * 6 + 6, sizeof(float))) == NULL)
		return DFX_NOMEM;
	err_cur = err;
	err_next = err + width * 3;

	/* channel order matches the B,G,R byte order of the sRGB image: */
	in[0] = B; in[1] = G; in[2] = R;

	/* produce sRGB pixel values: */
	for (y = 0; y < height; y++)
	{
		for (x = 0; x < width; x++) for (c = 0; c < 3; c++)
		{
			/* input value plus the error diffused into this pixel: */
			v = in[c][(p + y) * w_lin + p + x] + err_cur[x * 3 + c];

			/* quantize: */
			sRGB[(height-1-y) * w_srgb + x * 3 + c] = quant_8bit(to_srgb(v));

			/* compute the quantization error in linear space: */
			d = v - to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + c]));

			/* distribute the error using Floyd-Steinberg weights: */
			if (x + 1 < width) err_cur[(x + 1) * 3 + c] += d * 7.f / 16.f;		/* FS (0,+1)  */
			if (y + 1 < height) {
				if (x > 0) err_next[(x - 1) * 3 + c] += d * 3.f / 16.f;			/* FS (+1,-1) */
				err_next[x * 3 + c] += d * 5.f / 16.f;								/* FS (+1,0)  */
				if (x + 1 < width) err_next[(x + 1) * 3 + c] += d * 1.f / 16.f;	/* FS (+1,+1) */
			}
		}

		/* clear the row padding bytes, so that the buffer is fully initialized: */
		for (x = width * 3; x < w_srgb; x++)
			sRGB[(height-1-y) * w_srgb + x] = 0;

		/* move to the next row: the next row's errors become current, and the new next row starts empty */
		swap = err_cur; err_cur = err_next; err_next = swap;
		memset(err_next, 0, (size_t)width * 3 * sizeof(float));
	}

	free(err);
	return DFX_SUCCESS;
}

/* dfx_srgb.c -- end of file */
