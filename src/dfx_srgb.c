/*!
 *  \file   dfx_srgb.c
 *  \brief  DFX library: gamma-conversions between sRGB and linear RGB formats
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.01
 *  \date    March 7, 2026
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
 *  \brief 8-bit uniform quantizer function:
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
	/* inverse sRGB gamma: */
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
  *  sRGB image is assumed to be vertically flipped (as per Microsoft's conventions in bitmap files)
  *  Output linear planes are padded to support subsequent filtering operations. 
  *
  *  \param[in]  sRGB     - pointer to an sRGB image
  *  \param[in]  R, G, B  - pointers to linear R, G, and B channels
  *  \param[in]  height   - image height
  *  \param[in]  width    - image width
  *  \param[in]  p        - padding parameter
  * 
  *  \returns    DFX_X error code
  */
int srgb_to_linear(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);	/* w_srgb = offset to the next line in srgb bitmap image */
	int w_lin = width + 2 * p;			    /* w_lin = width of padded liner RGB image */
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
	pad_image(R, G, B, width, height, p, PAD_REPLICATE);

	return DFX_SUCCESS;
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation.
 *
 *  sRGB gamma is applied, and the resulting values quantized to 8-bit outputs.
 *  Padding is removed. The output consists of packed 24bit sRGB values, suitable for writing in a bitmap file.
 *
 *  \param[in]  R, G, B  - pointers to linear R, G, and B channels
 *  \param[out] sRGB     - pointer to an sRGB image
 *  \param[in]  height   - image height
 *  \param[in]  width    - image width
 *  \param[in]  p        - padding parameter
 *
 *  \returns    DFX_X error code
 */
int linear_to_srgb(float* R, float* G, float* B, unsigned char* sRGB, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);		/* w_srgb = offset to the next line in bitmap image */
	int w_lin = 2 * p + width;					/* w_lin = width of padded liner RGB image */
	int x, y;

	/* check parameters: */
	if (sRGB == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* produce sRGB pixel values: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++)
	{
		sRGB[(height-1-y) * w_srgb + x * 3 + 0] = quant_8bit(to_srgb(B[(p + y) * w_lin + p + x]));
		sRGB[(height-1-y) * w_srgb + x * 3 + 1] = quant_8bit(to_srgb(G[(p + y) * w_lin + p + x]));
		sRGB[(height-1-y) * w_srgb + x * 3 + 2] = quant_8bit(to_srgb(R[(p + y) * w_lin + p + x]));
	}

	return DFX_SUCCESS;
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation.
 *
 *  sRGB gamma is applied, and the resulting values quantized to 8-bit outputs.
 *  Banding is minimized by using Floyd-Steinberg-style diffusion of the quantization errors.
 *  Padding is removed. The output consists of packed 24bit sRGB values, suitable for writing in bitmap file.
 *
 *  \param[in]  R, G, B  - pointers to linear R, G, and B channels
 *  \param[out] sRGB     - pointer to an sRGB image
 *  \param[in]  height - image height
 *  \param[in]  width  - image width
 *  \param[in]  p      - padding parameter
 *
 *  \returns    DFX_X error code
 */
int linear_to_srgb_dithered(float* R, float* G, float* B, unsigned char* sRGB, int width, int height, int p)
{
	int w_srgb = (width * 3 + 3) & (~3);	/* w_srgb = offset to the next line in bitmap image */
	int w_lin = 2 * p + width;			    /* w_lin = width of padded liner RGB image */
	float dR, dG, dB;
	int x, y;

	/* check parameters: */
	if (sRGB == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0) 
		return DFX_INVARG;

	/* produce sRGB pixel values: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++)
	{
		/* compute current sRGB pixel: */
		sRGB[(height-1-y) * w_srgb + x * 3 + 0] = quant_8bit(to_srgb(B[(p + y) * w_lin + p + x]));
		sRGB[(height-1-y) * w_srgb + x * 3 + 1] = quant_8bit(to_srgb(G[(p + y) * w_lin + p + x]));
		sRGB[(height-1-y) * w_srgb + x * 3 + 2] = quant_8bit(to_srgb(R[(p + y) * w_lin + p + x]));

		/* distriibute noise for all rows and columns, except the last ones: */
		if (y < height - 1 && x < width - 1) 
		{
			/* compute quant_8bit error: */
			dB = B[(p + y) * w_lin + p + x] - to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 0]));
			dG = G[(p + y) * w_lin + p + x] - to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 1]));
			dR = R[(p + y) * w_lin + p + x] - to_linear(rec_8bit(sRGB[(height-1-y) * w_srgb + x * 3 + 2]));

			/* distribute noise using Floyd-Steinberg diffusion filter: */
			B[(p + y) * w_lin + p + x + 1] += dB * 7.f / 16.f;     /* FS (0,+1) */
			G[(p + y) * w_lin + p + x + 1] += dG * 7.f / 16.f;
			R[(p + y) * w_lin + p + x + 1] += dR * 7.f / 16.f;
			B[(p + y + 1) * w_lin + p + x] += dB * 5.f / 16.f;     /* FS (+1,0) */
			G[(p + y + 1) * w_lin + p + x] += dG * 5.f / 16.f;
			R[(p + y + 1) * w_lin + p + x] += dR * 5.f / 16.f;
			B[(p + y + 1) * w_lin + p + x - 1] += dB * 3.f / 16.f; /* FS (+1,-1) */
			G[(p + y + 1) * w_lin + p + x - 1] += dG * 3.f / 16.f;
			R[(p + y + 1) * w_lin + p + x - 1] += dR * 3.f / 16.f;
			B[(p + y + 1) * w_lin + p + x + 1] += dB * 1.f / 16.f; /* FS (+1,+1) */
			G[(p + y + 1) * w_lin + p + x + 1] += dG * 1.f / 16.f;
			R[(p + y + 1) * w_lin + p + x + 1] += dR * 1.f / 16.f;
		}
	}

	return DFX_SUCCESS;
}

/* dfx_srgb.c -- end of file */
