/*!
 *  \file   dfx_image.c
 *  \brief  DFX library: basic operations with linear RGB images
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
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
static unsigned int plane_size(int width, int height, int p)
{
	unsigned int size = (2 * p + height) * (2 * p + width);	        /* size of padded plane, in pixels */
	return size;
}

/*!
 *  \brief Allocate a channel of linear RGB image.
 *
 *  \param[in,out]  pX     - pointer to a pointer to allocated plane
 *  \param[in]      height - image height
 *  \param[in]      width  - image width
 *  \param[in]      p      - padding parameter
 *
 *  \returns        DFX_X error codes
 */
int alloc_plane(float** pX, int width, int height, int p)
{
	unsigned int size = plane_size(width, height, p);
	if (pX == NULL || height < 0 || width < 0 || p < 0)  return DFX_INVARG;
	if ((*pX = (float*)malloc(size * sizeof(float))) == NULL) return DFX_NOMEM;
	return DFX_SUCCESS;
}

/*!
 *  \brief Allocate linear RGB image.
 *
 *  \param[in,out]  pR, pG, pB  - pointers to a pointers to allocated R,G,B planes
 *  \param[in]      height - image height
 *  \param[in]      width  - image width
 *  \param[in]      p      - padding parameter
 *
 *  \returns        DFX_X error codes
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

	/* success: */
	*pR = R; *pG = G; *pB = B;
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
 *  Zero & copy operations:
 *
 *   zero_plane()
 *   zero_image()
 *   copy_plane()
 *   copy_image()
 */

/*!
 *  \brief Fill channel of linear RGB image with zeros.
 */
int zero_plane(float* X, int width, int height, int p)
{
	unsigned int size = plane_size(width, height, p); 
	if (X == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	/* setting all bytes to 0 is equivalent to assigning positive floating-point 0s */
	memset((void*)X, 0, size * sizeof(float));
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill linear RGB image with zeros.
 */
int zero_image(float* R, float* G, float* B, int width, int height, int p)
{
	if (R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	zero_plane(R, width, height, p);
	zero_plane(G, width, height, p);
	zero_plane(B, width, height, p);
	return DFX_SUCCESS;
}

/*!
 *  \brief Copy channel/plain of a linear RGB image.
 */
int copy_plane(float* X_in, float* X_out, int width, int height, int p)
{
	unsigned int size = plane_size(width, height, p);
	if (X_in == NULL || X_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	memcpy((void*)X_out, (void*)X_in, size * sizeof(float));
	return DFX_SUCCESS;
}

/*!
 *  \brief Copy linear RGB image.
 */
int copy_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p)
{
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	copy_plane(R_in, R_out, width, height, p);
	copy_plane(G_in, G_out, width, height, p);
	copy_plane(B_in, B_out, width, height, p);
	return DFX_SUCCESS;
}

/*******************
 * 
 *  Image padding functions 
 * 
 *   pad_plane()
 *   pad_image()
 */

/*!
 *  \brief Fill boundary for a given image plane/channel.
 * 
 *  \param[in,out]  X      - pointer to a channel to process
 *  \param[in]      height - image height
 *  \param[in]      width  - image width
 *  \param[in]      p      - padding parameter
 *  \param[in]      t      - padding type
 *
 *  \returns        DFX_X error codes
 */
int pad_plane(float* X, int width, int height, int p, int t)
{
	int w_lin = width + 2 * p;			    /* w_lin = width of padded liner RGB image */
	int x, y;

	/* check parameters: */
	if (X == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	
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
			X[(p + y) * w_lin + p - 1 - x] = X[(p + y) * w_lin + p + x];                  /* flip horisontal */
			X[(p + y) * w_lin + p + width + x] = X[(p + y) * w_lin + p + width - 1 - x];
		}
		/* repflect-pad top and bottom edges: */
		for (y = 0; y < p; y++) for (x = 0; x < w_lin; x++) {
			X[(p - 1 - y) * w_lin + x] = X[(p + y) * w_lin + x];                          /* flip vertical */
			X[(p + height + y) * w_lin + x] = X[(p + height - 1 - y) * w_lin + x];
		}
	}

	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill boundary for an image.
 */
int pad_image(float* R, float* G, float* B, int width, int height, int p, int t)
{
	if (R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0 || t < 0) return DFX_INVARG;
	pad_plane(R, width, height, p, t);
	pad_plane(G, width, height, p, t);
	pad_plane(B, width, height, p, t);
	return DFX_SUCCESS;
}

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
  *  sRGB image is considered vertically flipped (as per Microsoft's convention in bitmap files)
  *  Output linear planes are padded to support subsequent filtering operations. 
  *  The effective size of a padded image is (2*p + height) x (2*p + width), where p is a padding parameter.
  *
  *  \param[in]  sRGB   - pointer to an sRGB image
  *  \param[out] R      - pointer to linear R channel
  *  \param[out] G      - pointer to linear G channel
  *  \param[out] B      - pointer to linear B channel
  *  \param[in]  height - image height
  *  \param[in]  width  - image width
  *  \param[in]  p      - padding parameter
  * 
  *  \returns    DFX_X error codes
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

	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation.
 *
 *  sRGB gamma is applied, and the resulting values quantized to 8-bit outputs.
 *  Padding is removed. The output consists of packed 24bit sRGB values, suitable for writing in a bitmap file.
 *
 *  \param[out] sRGB  - pointer to an sRGB image
 *  \param[in]  R      - pointer to linear R channel
 *  \param[in]  G      - pointer to linear G channel
 *  \param[in]  B      - pointer to linear B channel
 *  \param[in]  height - image height
 *  \param[in]  width  - image width
 *  \param[in]  p      - padding parameter
 *
 *  \returns    DFX_X error codes
 */
int linear_to_srgb(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p)
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

	/* success: */
	return DFX_SUCCESS;
}

/*!
 *  \brief Converts linear-space RGB image to an 8-bit-per channel sRGB representation.
 *
 *  sRGB gamma is applied, and the resulting values quantized to 8-bit outputs.
 *  Banding is minimized by using Floyd-Steinberg-style diffusion of the quantization errors.
 *  Padding is removed. The output consists of packed 24bit sRGB values, suitable for writing in bitmap file.
 *
 *  \param[out] sRGB  - pointer to an sRGB image
 *  \param[in]  R      - pointer to linear R channel
 *  \param[in]  G      - pointer to linear G channel
 *  \param[in]  B      - pointer to linear B channel
 *  \param[in]  height - image height
 *  \param[in]  width  - image width
 *  \param[in]  p      - padding parameter
 *
 *  \returns    DFX_X error codes
 */
int linear_to_srgb_dithered(unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p)
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
			B[(p + y) * w_lin + p + x + 1] += dR * 7.f / 16.f;     /* FS (0,+1) */
			G[(p + y) * w_lin + p + x + 1] += dG * 7.f / 16.f;
			R[(p + y) * w_lin + p + x + 1] += dR * 7.f / 16.f;
			B[(p + y + 1) * w_lin + p + x] += dR * 5.f / 16.f;     /* FS (+1,0) */
			G[(p + y + 1) * w_lin + p + x] += dG * 5.f / 16.f;
			R[(p + y + 1) * w_lin + p + x] += dR * 5.f / 16.f;
			B[(p + y + 1) * w_lin + p + x - 1] += dR * 3.f / 16.f; /* FS (+1,-1) */
			G[(p + y + 1) * w_lin + p + x - 1] += dG * 3.f / 16.f;
			R[(p + y + 1) * w_lin + p + x - 1] += dR * 3.f / 16.f;
			B[(p + y + 1) * w_lin + p + x + 1] += dR * 1.f / 16.f; /* FS (+1,+1) */
			G[(p + y + 1) * w_lin + p + x + 1] += dG * 1.f / 16.f;
			R[(p + y + 1) * w_lin + p + x + 1] += dR * 1.f / 16.f;
		}
	}
	/* success: */
	return DFX_SUCCESS;
}

/**************
 *
 *  Limunance-related functions:
 *
 *   linear_to_luminance()
 *	 luminance_to_grayscale_rgb()
 */

/*!
  *  \brief Extracts luminance (XYZ's Y channel) from linear-space RGB image.
  *
  *  Conversion to luminance is accomplished by assuming that RGB space uses BT.709 primaries with D65 white.
  *
  *  \param[out] Y      - pointer to an luminance channel
  *  \param[in]  R      - pointer to linear R channel
  *  \param[in]  G      - pointer to linear G channel
  *  \param[in]  B      - pointer to linear B channel
  *  \param[in]  height - image height
  *  \param[in]  width  - image width
  *  \param[in]  p      - padding boundary to add to the linear image
  *
  *  \returns    DFX_X error codes
  */
int linear_to_luminance(float* Y, float* R, float* G, float* B, int width, int height, int p)
{
	int h_lin = 2 * p + height;			    /* h_lin = height of padded liner RGB image */
	int w_lin = 2 * p + width;			    /* w_lin = width of padded liner RGB image */
	int x, y;

	/* check parameters: */
	if (Y == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* compute luminance image: */
	for (y = 0; y < h_lin; y++) for (x = 0; x < w_lin; x++)
	{
		/* downmix to luminance based on BT.709 primaries: */
		Y[y * w_lin + x] = (float)(0.2126 * R[y * w_lin + x] + 0.7152 * G[y * w_lin + x] + 0.0722 * B[y * w_lin + x]);
	}

	/* success: */
	return DFX_SUCCESS;
}

/*!
  *  \brief Converts luminance to linear-space gray-scale RGB image.
  *
  *  \param[in]  Y      - pointer to an luminance channel
  *  \param[out] R      - pointer to linear R channel
  *  \param[out] G      - pointer to linear G channel
  *  \param[out] B      - pointer to linear B channel
  *  \param[in]  height - image height
  *  \param[in]  width  - image width
  *  \param[in]  p      - padding boundary to add to the linear image
  */
int luminance_to_grayscale_image(float* Y, float* R, float* G, float* B, int width, int height, int p)
{
	/* check parameters: */
	if (Y == NULL || R == NULL || G == NULL || B == NULL
		|| height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* replicate luminance values to all 3 channels: */
	copy_plane(Y, R, width, height, p);
	copy_plane(Y, G, width, height, p);
	copy_plane(Y, B, width, height, p);

	/* success: */
	return DFX_SUCCESS;
}

/* dfx_linear.c -- end of file */
