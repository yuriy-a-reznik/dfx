/*!
 *  \file   dfx_colors.c
 *  \brief  Color-related functions. 
 * 
 *  This module implements conversions between linear RGB and CIE XYZ, L*a*b*, and L*u*v* spaces.
 *  It also implements several basic color adjustment operations, such as hue rotation, saturation 
 *  adjustment, and brightness adjustment.
 * 
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.03
 *  \date    November 8, 2026
 */

#include <math.h>
#include <float.h>	
#include <stdlib.h>
#include "dfx.h"

 /*********************
  *
  *  Primary colors of practical RGB color systems:
  *
  *    xy_prim[]          - table of xy primary colors of RGB color systems
  *    uv_prim[]          - table of uv primary colors of RGB color systems
  *    gen_uv_prim()      - function that generates the uv_prim[] table
  */

/* Tolerance for the [0,1] range check of linear RGB values (about 1/30 of the smallest 8-bit sRGB step in linear light): */
#define RGB_RANGE_TOLERANCE  1e-5f

/* CIE 1931 (x,y) coordinates */
typedef struct { float x, y; } xy_t;

/* CIE 1976 (u',v') coordinates: */
typedef struct { float u, v; } uv_t;

/* Primary colors of standard RGB color systems (indexed by COLORS_* values)  */
static xy_t xy_prim[N_COLORS][4] = {
	/*** red: ******/  /*** green: ***/  /*** blue: ****/  /*** white: *****/
	{{0.640f, 0.330f}, {0.300f, 0.600f}, {0.150f, 0.060f}, {0.3127f, 0.3290f}}, /* BT.709, sRGB */
	{{0.670f, 0.330f}, {0.210f, 0.710f}, {0.140f, 0.080f}, {0.3100f, 0.3160f}}, /* 1953 NTSC legacy colors */
	{{0.630f, 0.340f}, {0.310f, 0.595f}, {0.155f, 0.070f}, {0.3127f, 0.3290f}}, /* SMPTE PR 145 (SMPTE-C) */
	{{0.640f, 0.330f}, {0.290f, 0.600f}, {0.150f, 0.060f}, {0.3127f, 0.3290f}}, /* PAL/SECAM, EBU Tech 3213 */
	{{0.708f, 0.292f}, {0.170f, 0.797f}, {0.131f, 0.046f}, {0.3127f, 0.3290f}}, /* BT.2020 */
	{{0.680f, 0.320f}, {0.265f, 0.690f}, {0.150f, 0.060f}, {0.3127f, 0.3290f}}, /* DCI-P3 D65 */
	{{0.680f, 0.320f}, {0.265f, 0.690f}, {0.150f, 0.060f}, {0.3140f, 0.3510f}}  /* DCI-P3 */
};

/* the same table, converted to (u'v') coordinates */
static uv_t uv_prim[N_COLORS][4];
static int uv_prim_generated = 0;

/*!
 *  \brief Generate a table of uv primary colors.
 */
static void gen_uv_prim(void)
{
	float x, y;
	int i, j;

	/* check if the table was already generated: */
	if (!uv_prim_generated) 
	{
		/* generate the uv_prim[] table: */
		for (i = 0; i < N_COLORS; i++) for (j = 0; j <= 3; j++) {
			/* pull (x,y) pair: */
			x = xy_prim[i][j].x;
			y = xy_prim[i][j].y;
			/* convert (x,y) to (u',v') (CIE 1976 UCS: u' = 4X/(X+15Y+3Z), v' = 9Y/(X+15Y+3Z), with X=x/y, Y=1, Z=(1-x-y)/y): */
			uv_prim[i][j].u = 4.f * x / (-2.f * x + 12.f * y + 3.f);
			uv_prim[i][j].v = 9.f * y / (-2.f * x + 12.f * y + 3.f);
		}
		uv_prim_generated = 1;
	}
}

/********************  
 * 
 * 3x3 matrix operations:
 * 
 *   matrix_vector_multiply()
 *   matrix_inverse()
 */

/*! 3x3 matrix-vector multiply: c = A b */
static void matrix_vector_multiply(float A[][3], const float b[], float c[])
{
	c[0] = A[0][0] * b[0] + A[0][1] * b[1] + A[0][2] * b[2];
	c[1] = A[1][0] * b[0] + A[1][1] * b[1] + A[1][2] * b[2];
	c[2] = A[2][0] * b[0] + A[2][1] * b[1] + A[2][2] * b[2];
}

/*! 3x3 matrix inversion by cofactors: W = M^-1 */
static int matrix_inverse(float M[][3], float W[][3])
{
	float a, b, c, d, e, f, g, h, i;
	double A, B, C, D, E, F, G, H, I;
	double det, det_inv;

	/* load matrix */
	a = M[0][0];  b = M[0][1];  c = M[0][2];
	d = M[1][0];  e = M[1][1];  f = M[1][2];
	g = M[2][0];  h = M[2][1];  i = M[2][2];

	/* compute cofactors */
	A = e * i - f * h;  D = -b * i + c * h;  G = b * f - c * e;
	B = -d * i + f * g;  E = a * i - c * g;  H = -a * f + c * d;
	C = d * h - e * g;  F = -a * h + b * g;  I = a * e - b * d;

	/* compute determinant */
	det = a * A + b * B + c * C;
	if (fabs(det) < FLT_MIN * 1000) return DFX_INVARG; 
	det_inv = 1.0 / det;

	/* store the inverse matrix */
	W[0][0] = (float)(A * det_inv);  W[0][1] = (float)(D * det_inv);  W[0][2] = (float)(G * det_inv);
	W[1][0] = (float)(B * det_inv);  W[1][1] = (float)(E * det_inv);  W[1][2] = (float)(H * det_inv);
	W[2][0] = (float)(C * det_inv);  W[2][1] = (float)(F * det_inv);  W[2][2] = (float)(I * det_inv);

	return DFX_SUCCESS;
}

/*********************
 *
 *  Functions supporting conversions between linear RGB systems and XYZ:
 * 
 *    compute_Ys()                 - compute Y values of primary colors given Y of white
 *    compute_rgb_to_xyz_matrix()  - computes RGB to XYZ conversion matrix given RGB primaries and Y-level of white
 */

/*!
 *  \brief Compute luminance values for primary colors.
 *
 *  The RGB->XYZ matrix of a color system has columns proportional to the XYZ coordinates
 *  of its primaries. The proportionality factors Y_r, Y_g, Y_b are found from the
 *  requirement that RGB = (1,1,1) maps to the white point with luminance Y_white.
 *
 *  \param[in]  p        - xy primaries + white point 
 *  \param[in]  Y_white  - luminance of white point
 *  \param[out] Y        - luminance values of primary colors in XYZ space
 *
 *  \return DFX_X error code
 */
static int compute_primary_Ys(const xy_t* p, float Y_white, float* Y)
{
	float M[3][3] = {0}, M_inv[3][3] = { 0 }, W[3] = { 0 };
	int i;

	/* check parameters: */
	if (p == NULL || Y == NULL || Y_white <= 0.f) return DFX_INVARG;

	/* XYZ coordinates of reference white (W): */
	W[0] = Y_white * (p[3].x / p[3].y);						/* X = Y*x/y  */
	W[1] = Y_white;											/* Y = Y      */
	W[2] = Y_white * ((1.0f - p[3].x - p[3].y) / p[3].y);	/* Z = Y*z/y  */

	/* normalized RGB to XYZ conversion matrix (M): */
	for (i = 0; i < 3; i++) {
		M[0][i] = p[i].x / p[i].y;                          /* x/y */
		M[1][i] = 1.0f;                                     /* y/y */
		M[2][i] = (1.0f - p[i].x - p[i].y) / p[i].y;        /* z/y */
	}

	/* solve: M * Y = W; */
	if ((i = matrix_inverse(M, M_inv)) != DFX_SUCCESS) return i;
	matrix_vector_multiply(M_inv, W, Y);

	return DFX_SUCCESS;
}

/*!
 *  \brief Compute mapping of an RGB space to XYZ.
 *
 *  \param[in]  p          - xy primaries + white point 
 *  \param[in]  Y_white    - luminance of white point
 *  \param[out] rgb_to_xyz - RGB to XYZ conversion matrix to generate
 *
 *  \returns    DFX_X error code
 */
static int compute_rgb_to_xyz_matrix(const xy_t* p, float Y_white, float rgb_to_xyz[][3])
{
	float Y[3] = {0};
	int i;

	/* check parameters: */
	if (p == NULL || rgb_to_xyz == NULL || Y_white <= 0.f) return DFX_INVARG;

	/* compute Y-values for primary colors */
	if ((i = compute_primary_Ys(p, Y_white, Y)) != DFX_SUCCESS) return i;

	/* define full RGB to XYZ conversion matrix: */
	for (i = 0; i < 3; i++) {
		rgb_to_xyz[0][i] = Y[i] * p[i].x / p[i].y;
		rgb_to_xyz[1][i] = Y[i];
		rgb_to_xyz[2][i] = Y[i] * (1.0f - p[i].x - p[i].y) / p[i].y;
	}

	return DFX_SUCCESS;
}

/*****************
 * 
 *  Conversions between linear RGB to XYZ color spaces:
 * 
 *   range()               - computes range of rgb values subject to conversion
 *   check_rgb_range()     - checks that RGB values are in [0,1] range
 *   linear_rgb_to_xyz()   - RGB->XYZ conversion
 *   xyz_to_linear_rgb()   - XYZ->RGB conversion
 */

 /*!
   *  \brief Compute range of sample values in an image plane.
   *
   *  \param[in]      X        - image plane
   *  \param[in]      width    - image width in pixels
   *  \param[in]      height   - image height in pixels
   *  \param[in]      p        - padding parameter
   *  \param[out]     p_X_min  - pointer to a variable to store minimum sample value
   *  \param[out]     p_X_max  - pointer to a variable to store maximum sample value
   *
   *  \returns        DFX_X error code
   */
static int range(const float* X, int width, int height, int p, float* p_X_min, float* p_X_max)
{
	int x, y;
	float v, v_min, v_max;

	/* check parameters: */
	if (X == NULL || p_X_min == NULL || p_X_max == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* compute max sample value: */
	v_min = FLT_MAX;
	v_max = -FLT_MAX;
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		v = X[(p + y) * (width + 2 * p) + p + x];
		if (v < v_min) v_min = v;
		if (v > v_max) v_max = v;
	}

	/* store the results & exit: */
	*p_X_min = v_min;
	*p_X_max = v_max;
	return DFX_SUCCESS;
}

/*!
 *  \brief Check if all samples in an RGB image are in the [0,1] range.
 * 
 *  RGB values are essentially the gains that must be applied to the primary colors to reproduce image. 
 *  Such gains are always normalized to the maximum range that can be stored in a particular format. 
 *  E.g. in 8-bit/channel sRGB bitmaps, the maximum value is 255. 
 *
 *  When we convert to linear RGB, we map them to [0,1] range. 
 *  We expect that proper filtering and adjustment operations will also keep them in the same range
 *  (filters with negative taps may overshoot slightly; use clip_image() in that case). 
 * 
 *  This function checks this range for cases when we are about to convert RGB images to 
 *  the XYZ space. If RGB values are outside of [0,1] range, such conversion may not be correct, 
 *  and the results may be meaningless.
 *
 *  A small tolerance (RGB_RANGE_TOLERANCE, well below one 8-bit quantization step) is allowed,
 *  so that rounding errors of round trips through XYZ and L*u*v* are not reported as range errors.
 *
 *  \param[in]   R, G, B  - image planes
 *  \param[in]   width    - image width in pixels
 *  \param[in]   height   - image height in pixels
 *  \param[in]   p        - padding parameter
 *
 *  \returns     DFX_SUCCESS if success, 
 *               DFX_INVRANGE if any of the RGB planes contains values outside of [0,1] range, or 
 *               DFX_INVARG if any of the parameters is invalid.
 */
static int check_rgb_range(const float* R, const float* G, const float* B, int width, int height, int p)
{ 
	float R_min, R_max, G_min, G_max, B_min, B_max;

	/* check parameters: */
	if (R == NULL || G == NULL || B == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;

	/* compute ranges of RGB values in each plane: */
	if (range(R, width, height, p, &R_min, &R_max) != DFX_SUCCESS) return DFX_INVARG;
	if (range(G, width, height, p, &G_min, &G_max) != DFX_SUCCESS) return DFX_INVARG;
	if (range(B, width, height, p, &B_min, &B_max) != DFX_SUCCESS) return DFX_INVARG;

	/* check ranges: */
	if (R_min < -RGB_RANGE_TOLERANCE || R_max > 1.0f + RGB_RANGE_TOLERANCE) return DFX_INVRANGE;
	if (G_min < -RGB_RANGE_TOLERANCE || G_max > 1.0f + RGB_RANGE_TOLERANCE) return DFX_INVRANGE;
	if (B_min < -RGB_RANGE_TOLERANCE || B_max > 1.0f + RGB_RANGE_TOLERANCE) return DFX_INVRANGE;

	return DFX_SUCCESS;
}

/*!
 *  \brief Linear RGB to XYZ conversion.
 * 
 *  Input RGB values are assumed to be limited to [0,1] range.
 *  Parameter Y_white defines the peak luminance to assign to the image in the XYZ space.
 *
 *  \param[in]  R, G, B    - R,G,B channels in input linear RGB image
 *  \param[out] X, Y, Z    - X,Y,Z channels in XYZ image
 *  \param[in]  width      - image width
 *  \param[in]  height     - image height
 *  \param[in]  p          - padding parameter
 *  \param[in]  colors     - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
 *  \param[in]  Y_white    - luminance of white point
 * 
 *  \returns    DFX_X error code
 */
int linear_rgb_to_xyz(const float* R, const float* G, const float* B, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white)
{ 
	float rgb_to_xyz[3][3] = {0}, rgb[3] = { 0 }, xyz[3] = { 0 };
	int x, y;

	/* check parameters: */
	if (R == NULL || G == NULL || B == NULL) return DFX_INVARG;
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* check the range of RGB values: */
	if ((x = check_rgb_range(R, G, B, width, height, p)) != DFX_SUCCESS) return x;

	/* compute RGB->XYZ conversion matrix: */
	if ((x = compute_rgb_to_xyz_matrix(xy_prim[colors], Y_white, rgb_to_xyz)) != DFX_SUCCESS) return x;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* retrieve RGB values: */
		rgb[0] = R[(p + y) * (width + 2 * p) + p + x];
		rgb[1] = G[(p + y) * (width + 2 * p) + p + x];
		rgb[2] = B[(p + y) * (width + 2 * p) + p + x];
		/* perform conversion: */
		matrix_vector_multiply(rgb_to_xyz, rgb, xyz);
		/* store the results: */
		X[(p + y) * (width + 2 * p) + p + x] = xyz[0];
		Y[(p + y) * (width + 2 * p) + p + x] = xyz[1];
		Z[(p + y) * (width + 2 * p) + p + x] = xyz[2];
	}
 
	/* fill the padding of the output planes & exit: */
	return pad_image(X, Y, Z, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief XYZ to linear RGB conversion.
 *
 *  The output is produced for all pixels; DFX_INVRANGE is returned if any of the resulting
 *  RGB values are outside of the [0,1] range, i.e. the color is outside of the RGB gamut.
 */
int xyz_to_linear_rgb(const float* X, const float* Y, const float* Z, float* R, float* G, float* B, int width, int height, int p, int colors, float Y_white)
{
	float rgb_to_xyz[3][3] = {0}, xyz_to_rgb[3][3] = { 0 }, rgb[3] = { 0 }, xyz[3] = { 0 };
	int x, y;

	/* check parameters: */
	if (R == NULL || G == NULL || B == NULL) return DFX_INVARG;
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* compute RGB->XYZ conversion matrix: */
	if ((x = compute_rgb_to_xyz_matrix(xy_prim[colors], Y_white, rgb_to_xyz)) != DFX_SUCCESS) return x;

	/* compute its inverse: */
	if ((x = matrix_inverse(rgb_to_xyz, xyz_to_rgb)) != DFX_SUCCESS) return x;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* retrieve XYZ values: */
		xyz[0] = X[(p + y) * (width + 2 * p) + p + x];
		xyz[1] = Y[(p + y) * (width + 2 * p) + p + x];
		xyz[2] = Z[(p + y) * (width + 2 * p) + p + x];
		/* perform conversion: */
		matrix_vector_multiply(xyz_to_rgb, xyz, rgb);
		/* store the results: */
		R[(p + y) * (width + 2 * p) + p + x] = rgb[0];
		G[(p + y) * (width + 2 * p) + p + x] = rgb[1];
		B[(p + y) * (width + 2 * p) + p + x] = rgb[2];
	}

	/* fill the padding of the output planes: */
	if ((x = pad_image(R, G, B, width, height, p, PAD_REPLICATE)) != DFX_SUCCESS) return x;

	/* check the range of RGB values & exit: */
	return check_rgb_range(R, G, B, width, height, p);
}

/**************
 *
 *  Luminance (XYZ Y-channel)-related functions:
 *
 *   linear_rgb_to_luminance()
 *   linear_to_luminance()
 *   luminance_to_grayscale_image()
 */

 /*!
  *  \brief Extracts luminance (XYZ's Y channel) from linear-space RGB image.
  *
  *  The luminance is the middle row of the RGB->XYZ matrix applied to (R,G,B); for sRGB primaries
  *  and D65 white this is Y = 0.2126 R + 0.7152 G + 0.0722 B (times Y_white).
  *
  *  \param[in]  R, G, B  - pointers to linear R, G, and B channels
  *  \param[out] Y        - pointer to a luminance plane
  *  \param[in]  width    - image width
  *  \param[in]  height   - image height
  *  \param[in]  p        - padding parameter
  *  \param[in]  colors   - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
  *  \param[in]  Y_white  - luminance of white point
  *
  *  \returns    DFX_X error code
  */
int linear_rgb_to_luminance(const float* R, const float* G, const float* B, float* Y, int width, int height, int p, int colors, float Y_white)
{
	float rgb_to_xyz[3][3] = { 0 }, r, g, b, xyz_y;
	int x, y;

	/* check parameters: */
	if (R == NULL || G == NULL || B == NULL || Y == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* check the range of RGB values: */
	if ((x = check_rgb_range(R, G, B, width, height, p)) != DFX_SUCCESS) return x;

	/* compute RGB->XYZ conversion matrix: */
	if ((x = compute_rgb_to_xyz_matrix(xy_prim[colors], Y_white, rgb_to_xyz)) != DFX_SUCCESS) return x;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* retrieve RGB values: */
		r = R[(p + y) * (width + 2 * p) + p + x];
		g = G[(p + y) * (width + 2 * p) + p + x];
		b = B[(p + y) * (width + 2 * p) + p + x];
		/* downmix to Y: */
		xyz_y = rgb_to_xyz[1][0] * r + rgb_to_xyz[1][1] * g + rgb_to_xyz[1][2] * b;
		/* store the results: */
		Y[(p + y) * (width + 2 * p) + p + x] = xyz_y;
	}

	/* fill the padding of the output plane & exit: */
	return pad_plane(Y, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Simplified version of linear_rgb_to_luminance(), assuming that the input is an sRGB-converted image, and white Y=1.
 */
int linear_to_luminance(const float* R, const float* G, const float* B, float* Y, int width, int height, int p)
{
	return linear_rgb_to_luminance(R, G, B, Y, width, height, p, COLORS_SRGB, 1.0f);
}

/*!
  *  \brief Converts luminance to linear-space gray-scale RGB image.
  *
  *  A gray color has R = G = B, and for sRGB primaries its luminance equals that common value
  *  (the middle row of the RGB->XYZ matrix sums to 1). The luminance plane is therefore expected
  *  to be normalized to Y_white = 1.
  *
  *  \param[in]  Y        - pointer to a luminance plane
  *  \param[out] R, G, B  - pointers to linear R, G, and B channels
  *  \param[in]  width    - image width
  *  \param[in]  height   - image height
  *  \param[in]  p        - padding parameter
  */
int luminance_to_grayscale_image(const float* Y, float* R, float* G, float* B, int width, int height, int p)
{
	int err;

	/* check parameters: */
	if (Y == NULL || R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* replicate luminance values to all 3 channels: */
	if ((err = copy_plane(Y, R, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = copy_plane(Y, G, width, height, p)) != DFX_SUCCESS) return err;
	return copy_plane(Y, B, width, height, p);
}

/*****************
 *
 *  Conversions between XYZ and CIE L*u*v* and CIE L*a*b* color spaces:
 *
 *   cie_f()      - compander function 
 *   cie_f_inv()  - inverse compander function, parameterized by L*
 *   cie_f_inv_t()- inverse compander function
 *   xyz_to_luv() - XYZ to L*u*v* conversion
 *   luv_to_xyz() - L*u*v* to XYZ conversion
 *   xyz_to_lab() - XYZ to L*a*b* conversion
 *   lab_to_xyz() - L*a*b* to XYZ conversion
 *
 *  Notes: These 1976 CIE standards have been recently revised, and their formulas have been updated.
 *  The old formulae listed in the books and many online sources are no longer correct.
 *  The current standards are:
 *    ISO/CIE 11664-4:2019 - CIE L*a*b* color space
 *    ISO/CIE 11664-5:2024 - CIE L*u*v* color space
 */

/*!
 *  \brief Compander used in computation of both L*a*b* and L*u*v* color spaces:
 *
 *      f(t) = t^(1/3)                  for t > (6/29)^3,
 *      f(t) = (841/108) t + 4/29       otherwise,
 *
 *  so that f is continuous with a continuous derivative at t = (6/29)^3 = 0.008856.
 */
static float cie_f(float y)
{
	const double y_lim = pow(6. / 29., 3);
	const double scale = 841. / 108., offs = 4./29.; 
	double f;

	/* conversion formula (ISO/CIE 11664-5:2024): */
	if (y > y_lim)  
		f = pow(y, 1. / 3.);    /* cubic root branch */
	else
		f = y * scale + offs;	/* linear branch */

	return (float) f;
}

/*!
 *  \brief Inverse compander process.
 *  
 *  This function finds y, such that cie_f(y) = (L_star + 16) / 116.
 *  In other words, is not pure inverse of cie_f(y) = x, but an inverse with a substitute 
 *  of a parameter (L_star instead of x). 
 * 
 *  This modification avoids the rounding errors introduced by direct inversion 
 *  and offset computations. 
 */
static float cie_f_inv(float L_star)
{
	const double L_star_lim = 8.;
	const double scale_inv = 108. / 841.;
	double y;

	/* inverse conversion formula: */
	if (L_star > L_star_lim)
		y = pow((L_star + 16.) / 116., 3.);
	else
		y = scale_inv * (L_star / 116);	   /* 16/116 - 4/29 cancel out !!! */

	return (float)y;
}

/*!
 *  \brief Inverse compander: finds t such that cie_f(t) = f.
 *
 *      t = f^3                         for f > 6/29,
 *      t = (108/841) (f - 4/29)        otherwise.
 */
static float cie_f_inv_t(float f)
{
	const double f_lim = 6. / 29.;
	const double scale_inv = 108. / 841., offs = 4. / 29.;
	double t;

	if (f > f_lim)
		t = pow(f, 3.);
	else
		t = scale_inv * (f - offs);

	return (float)t;
}

/*!
 *  \brief XYZ to CIE L*u*v* color space conversion.
 * 
 *   This implementation follows the algorithm defined by ISO/CIE 11664-5:2024 standard. 
 *   The division by zero in X,Y,Z -> 0 region is avoided by forcing (x,y) and (u',v') coordinates 
 *   of such colors to white. 
 * 
 *  \param[in]  X, Y, Z    - X,Y,Z channels in XYZ image
 *  \param[out] L, u, v    - L*, u*, and v* channels in CIE L*u*v* space
 *  \param[in]  height     - image height
 *  \param[in]  width      - image width
 *  \param[in]  p          - padding parameter
 *  \param[in]  colors     - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
 *  \param[in]  Y_white    - luminance of white point
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int xyz_to_luv(const float* X, const float* Y, const float* Z, float* L, float* u, float* v, int width, int height, int p, int colors, float Y_white)
{
	float X_in, Y_in, Z_in, u_norm, u_prime, v_prime, L_star, u_star, v_star;
	float u_prime_white, v_prime_white;
	int x, y;

	/* check parameters: */
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (L == NULL || u == NULL || v == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* get (u',v') coordinates of white: */
	gen_uv_prim();
	u_prime_white = uv_prim[colors][3].u;
	v_prime_white = uv_prim[colors][3].v;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* retrieve XYZ values: */
		X_in = X[(p + y) * (width + 2 * p) + p + x];
		Y_in = Y[(p + y) * (width + 2 * p) + p + x];
		Z_in = Z[(p + y) * (width + 2 * p) + p + x];
		/* compute uv chromaticity coordinates: */
		u_norm = X_in + 15.f * Y_in + 3.f * Z_in;
		if (u_norm > 1e-6) {
			/* normal u', v' conversion: */
			u_prime = 4.f * X_in / u_norm;
			v_prime = 9.f * Y_in / u_norm;
		} else { 
			/* force black to be white: */
			u_prime = u_prime_white;
			v_prime = v_prime_white;
		}
		/* convert to L*u*v*: */
		L_star = 116.f * cie_f(Y_in / Y_white) - 16.f;
		u_star = 13.f * L_star * (u_prime - u_prime_white);
		v_star = 13.f * L_star * (v_prime - v_prime_white);
		/* store the results: */
		L[(p + y) * (width + 2 * p) + p + x] = L_star;
		u[(p + y) * (width + 2 * p) + p + x] = u_star;
		v[(p + y) * (width + 2 * p) + p + x] = v_star;
	}

	/* fill the padding of the output planes & exit: */
	return pad_image(L, u, v, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief CIE L*u*v* to XYZ conversion.
 *  
 *   This function implements the inverse conversion (ISO/CIE 11664-5:2024, Annex A).
 *   The division by zero in X,Y,Z -> 0 region is avoided by forcing (x,y) and (u',v') coordinates 
 *   of such colors to white. 
 * 
 *   The function fails with DFX_INVRANGE error if (x,y)s of the reconstruction point to non-existent colors.
 *
 *  \param[in]  L, u, v    - L*, u*, and v* channels in CIE L*u*v* image
 *  \param[out] X, Y, Z    - X,Y,Z channels in XYZ image
 *  \param[in]  width      - image width
 *  \param[in]  height     - image height
 *  \param[in]  p          - padding parameter
 *  \param[in]  colors     - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
 *  \param[in]  Y_white    - luminance of white point
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int luv_to_xyz(const float* L, const float* u, const float* v, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white)
{
	float L_star, u_star, v_star;
	float X_out, Y_out, Z_out, x_out, y_out, u_prime, v_prime;
	float u_prime_white, v_prime_white;
	int x, y;

	/* check parameters: */
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (L == NULL || u == NULL || v == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* get (u',v') coordinates of white: */
	gen_uv_prim();
	u_prime_white = uv_prim[colors][3].u;
	v_prime_white = uv_prim[colors][3].v;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* retrieve L*u*v* values: */
		L_star = L[(p + y) * (width + 2 * p) + p + x];
		u_star = u[(p + y) * (width + 2 * p) + p + x];
		v_star = v[(p + y) * (width + 2 * p) + p + x];
		/* recover Y, u_prime, v_prime: */
		Y_out = Y_white * cie_f_inv(L_star);
		if (L_star > 1e-4) {
			/* normal u',v' conversion: */
			u_prime = u_star / (13.f * L_star) + u_prime_white;
			v_prime = v_star / (13.f * L_star) + v_prime_white;
		} else { 
			/* force black to be white: */
			u_prime = u_prime_white;
			v_prime = v_prime_white;
		}
		/* convert to (x,y): */
		x_out = 9.f * u_prime / (6.f * u_prime - 16.f * v_prime + 12.f);
		y_out = 4.f * v_prime / (6.f * u_prime - 16.f * v_prime + 12.f);   
		/* sanity check (0.0050 is y of 380nm blue light): */
		if (y_out < 0.0049999) 
			return DFX_INVRANGE;
		/* recover X,Y: */
		X_out = Y_out * x_out / y_out;
		Z_out = Y_out * (1.f - x_out - y_out) / y_out;
		/* store the results: */
		X[(p + y) * (width + 2 * p) + p + x] = X_out;
		Y[(p + y) * (width + 2 * p) + p + x] = Y_out;
		Z[(p + y) * (width + 2 * p) + p + x] = Z_out;
	}

	/* fill the padding of the output planes & exit: */
	return pad_image(X, Y, Z, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief XYZ to CIE L*a*b* color space conversion.
 *
 *   This implementation follows ISO/CIE 11664-4:2019:
 *
 *      L* = 116 f(Y/Y_n) - 16,
 *      a* = 500 [f(X/X_n) - f(Y/Y_n)],
 *      b* = 200 [f(Y/Y_n) - f(Z/Z_n)],
 *
 *   where (X_n, Y_n, Z_n) are the tristimulus values of the white point of the selected color
 *   system scaled to Y_n = Y_white, and f is the compander cie_f().
 *
 *  \param[in]  X, Y, Z    - X,Y,Z channels in XYZ image
 *  \param[out] L, a, b    - L*, a*, and b* channels in CIE L*a*b* image
 *  \param[in]  width      - image width
 *  \param[in]  height     - image height
 *  \param[in]  p          - padding parameter
 *  \param[in]  colors     - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
 *  \param[in]  Y_white    - luminance of white point
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int xyz_to_lab(const float* X, const float* Y, const float* Z, float* L, float* a, float* b, int width, int height, int p, int colors, float Y_white)
{
	float X_white, Z_white, f_x, f_y, f_z;
	int x, y;

	/* check parameters: */
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (L == NULL || a == NULL || b == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* XYZ coordinates of the white point: */
	X_white = Y_white * xy_prim[colors][3].x / xy_prim[colors][3].y;
	Z_white = Y_white * (1.f - xy_prim[colors][3].x - xy_prim[colors][3].y) / xy_prim[colors][3].y;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* compand normalized tristimulus values: */
		f_x = cie_f(X[(p + y) * (width + 2 * p) + p + x] / X_white);
		f_y = cie_f(Y[(p + y) * (width + 2 * p) + p + x] / Y_white);
		f_z = cie_f(Z[(p + y) * (width + 2 * p) + p + x] / Z_white);
		/* store L*, a*, b*: */
		L[(p + y) * (width + 2 * p) + p + x] = 116.f * f_y - 16.f;
		a[(p + y) * (width + 2 * p) + p + x] = 500.f * (f_x - f_y);
		b[(p + y) * (width + 2 * p) + p + x] = 200.f * (f_y - f_z);
	}

	/* fill the padding of the output planes & exit: */
	return pad_image(L, a, b, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief CIE L*a*b* to XYZ conversion.
 *
 *   Inverse of xyz_to_lab() (ISO/CIE 11664-4:2019, Annex A):
 *
 *      f_y = (L* + 16) / 116,   f_x = f_y + a* / 500,   f_z = f_y - b* / 200,
 *      X = X_n f^-1(f_x),   Y = Y_n f^-1(f_y),   Z = Z_n f^-1(f_z).
 *
 *  \param[in]  L, a, b    - L*, a*, and b* channels in CIE L*a*b* image
 *  \param[out] X, Y, Z    - X,Y,Z channels in XYZ image
 *  \param[in]  width      - image width
 *  \param[in]  height     - image height
 *  \param[in]  p          - padding parameter
 *  \param[in]  colors     - CIE 1931 xy primaries + white point to use (see COLORS_* enum)
 *  \param[in]  Y_white    - luminance of white point
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int lab_to_xyz(const float* L, const float* a, const float* b, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white)
{
	float X_white, Z_white, f_x, f_y, f_z;
	int x, y;

	/* check parameters: */
	if (X == NULL || Y == NULL || Z == NULL) return DFX_INVARG;
	if (L == NULL || a == NULL || b == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;
	if (Y_white <= 0.f) return DFX_INVARG;

	/* XYZ coordinates of the white point: */
	X_white = Y_white * xy_prim[colors][3].x / xy_prim[colors][3].y;
	Z_white = Y_white * (1.f - xy_prim[colors][3].x - xy_prim[colors][3].y) / xy_prim[colors][3].y;

	/* perform the conversion: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* recover companded values: */
		f_y = (L[(p + y) * (width + 2 * p) + p + x] + 16.f) / 116.f;
		f_x = f_y + a[(p + y) * (width + 2 * p) + p + x] / 500.f;
		f_z = f_y - b[(p + y) * (width + 2 * p) + p + x] / 200.f;
		/* invert the compander and scale by the white point: */
		X[(p + y) * (width + 2 * p) + p + x] = X_white * cie_f_inv_t(f_x);
		Y[(p + y) * (width + 2 * p) + p + x] = Y_white * cie_f_inv_t(f_y);
		Z[(p + y) * (width + 2 * p) + p + x] = Z_white * cie_f_inv_t(f_z);
	}

	/* fill the padding of the output planes & exit: */
	return pad_image(X, Y, Z, width, height, p, PAD_REPLICATE);
}

/***************
 * 
 *  Color adjustment operations:  
 *    
 *    gain_sat()           - gain + saturation function
 *    luv_sat_scale()      - computes saturation and gap to gamut boundary in L*u*v* space
 *	  adjust_brightness()  - brightness adjustment function
 *    adjust_saturation()  - saturation adjustment function
 *    adjust_hue()         - hue adjustment function
 */

 /*!
  *  \brief Gain + saturation function.
  *
  *  Apply gain factor alpha, and ensure that the result stays in the [0,max] range. 
  *  General formula is:
  * 
  *     f(x) = (max^(-k) + (alpha * x)^(-k) )^(-1/k)
  * 
  *  When k=0, we replace it by hard clipping:
  * 
  *     f(x) = min(alpha * x, max) 
  *
  *  \param[in]  x       - input value (must be the range [0,1])
  *  \param[in]  alpha   - gain factor (must be > 0)
  *  \param[in]  max     - maximum allowed value on output
  *  \param[in]  k       - smoothness parameter (0 = hard clipping, 1+ - smooth variants)
  *
  *  \returns    f(x)    - gain adjusted and saturated value
  */
static float gain_sat(float x, float alpha, float max, int k)
{
	float y;
	
	/* constrain parameters: */
	if (alpha < 0.f) alpha = 0.f;
	if (x < 0.f) x = 0.f; 
	
	/* apply gain & limit it to max:  */
	y = alpha * x;
	if (k <= 0 || x < 1e-6)
		y = (y < max) ? y : max;  /* use hard clipping */
	else 
		y = powf(powf(max, -(float)k) + powf(y, -(float)k), -1.f / (float)k);
	
	return y;
}

/*!
 *  \brief Compute L*u*v* color saturation and a scale that is needed to map this color to the RGB gamut boundary.
 *
 *  The input point may be inside or outside of the RGB gamut.
 *  Scale < 1 will indicate that it will need to be desaturated to be brought back to gamut.
 *  Scale > 1 will indicate that there is still some gap between this point and the gamut boundary.
 *  For near black or white points this function returns saturation 0 and scale 1.
 *
 *  The function relies on uv_prim[] values and assumes that this table is already generated.
 *
 *  \param[in]  L_star, u_star, v_star - input color (x)
 *  \param[out] p_sat                  - variable to receive saturation of this color
 *  \param[out] p_scale                - variable to received the scale to boundary
 *  \param[in]  colors                 - color space index (see COLORS_X enum)
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
static int luv_sat_scale(float L_star, float u_star, float v_star, float* p_sat, float* p_scale, int colors)
{
	float delta_xw_u, delta_xw_v, delta_ji_u, delta_ji_v, delta_iw_u, delta_iw_v, D;
	float s_uv, scale, min_scale;
	int i, j;

	/* check parameters: */
	if (p_sat == NULL || p_scale == NULL) return DFX_INVARG;
	if (colors < 0 || colors >= N_COLORS) return DFX_INVARG;

	/* are we near black? */
	if (L_star < 1e-4) {
		/* assume 0 saturation, unit scale & return: */
		*p_sat = 0.0;
		*p_scale = 1.0f;
		return DFX_SUCCESS;
	}

	/* compute delta_u,v and saturation: */
	delta_xw_u = u_star / (13.f * L_star);
	delta_xw_v = v_star / (13.f * L_star);
	s_uv = (float)sqrt(delta_xw_u * delta_xw_u + delta_xw_v * delta_xw_v);

	/* check if s_uv ~ 0 */
	if (s_uv < 1e-4) {
		*p_sat = s_uv;
		*p_scale = 1.0f;
		return DFX_SUCCESS;
	}

	/* 
	 * Scan boundaries of the RGB triangle, and find the smallest positive 
	 * scale that needs to be applied to a vector from white to our input point (x) 
	 * such that it touches one of those boundaries. 
	 */
	min_scale = FLT_MAX;
	for (i = 0; i < 3; i++)
	{
		/* pick the 2nd edge point: */
		j = i + 1; if (j == 3) j = 0;
		/* compute determinant: */
		delta_ji_u = uv_prim[colors][j].u - uv_prim[colors][i].u;
		delta_ji_v = uv_prim[colors][j].v - uv_prim[colors][i].v;
		D = delta_xw_u * delta_ji_v - delta_xw_v * delta_ji_u;
		/* check if the w->x is near parallel to line j->i: */
		if (fabs(D) < 1e-6) continue;
		/* compute the scale that is needed to reach the intersection point: */
		delta_iw_u = uv_prim[colors][i].u - uv_prim[colors][3].u;
		delta_iw_v = uv_prim[colors][i].v - uv_prim[colors][3].v;
		scale = (delta_iw_u * delta_ji_v - delta_iw_v * delta_ji_u) / D;
		/*******
		 * the coordinates of the intersection point are:
		 *   u_x_edge = uv_prim[colors][3].u + scale * delta_xw_u;
		 *   v_x_edge = uv_prim[colors][3].v + scale * delta_xw_v;
		 **/
		/* reject vectors positioned in the opposite directions: */
		if (scale <= 0.f) continue;
		/* pick the smallest positive scale needed to hit the triangle: */
		if (scale < min_scale)
			min_scale = scale;
	}

	/* set scale to 1 if the above search was not successful: */
	if (min_scale == FLT_MAX)
		min_scale = 1.0f;

	/* report the results and exit: */
	*p_sat = s_uv;
	*p_scale = min_scale;
	return DFX_SUCCESS;
}

/*!
 *  \brief Brightness adjustment function.
 * 
 *  Everything is done in CIE L*u*v* space for the sake of conceptual simplicity. 
 *  Brightness-scaled L* channel values are companded by the gain_sat() function. This is a non-linear operation. 
 *  We then compute the effective scale factor that became realized in the L* channel and apply it to u* and v*,
 *  which keeps the chromaticity (u',v') of each pixel unchanged (since u* = 13 L* (u' - u'_n)). 
 * 
 *  \param[in]  R_in, G_in, B_in     - input linear RGB image
 *  \param[out] R_out, G_out, B_out  - output RGB image
 *  \param[in]  width                - image width
 *  \param[in]  height               - image height
 *  \param[in]  p                    - padding parameter
 *  \param[in]  brightness           - brightness scale adjustment parameter (0..4)
 *  \param[in]  k                    - clipping shape parameter (0 - hard clipping, 1+ - use generalized harmonic mean compander)
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise (DFX_INVRANGE if some colors fall outside of the RGB gamut)
 */
int adjust_brightness(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float brightness, int k)
{
	/* local constants: */
	const int colors = COLORS_SRGB;			/* input/output color space */
	const float Y_white = 1.0f;				/* luminance of white, used in RGB<->XYZ conversions */

	/* temp buffers & local variables: */
	float *X = NULL, *Y = NULL, *Z = NULL;
	float *L = NULL, *u = NULL, *v = NULL;
	float L_in, u_in, v_in, L_out, u_out, v_out, L_scale;
	int x, y, err = DFX_SUCCESS;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (brightness < 0.f || brightness > 4.0f) return DFX_INVARG;

	/* allocate intermediate frames: */
	if (alloc_image(&X, &Y, &Z, width, height, p) != DFX_SUCCESS) return DFX_NOMEM;
	if (alloc_image(&L, &u, &v, width, height, p) != DFX_SUCCESS) { free_image(X, Y, Z); return DFX_NOMEM; }

	/* convert RGB image to XYZ and then Luv: */
	if ((err = linear_rgb_to_xyz(R_in, G_in, B_in, X, Y, Z, width, height, p, colors, Y_white)) == DFX_SUCCESS
	 && (err = xyz_to_luv(X, Y, Z, L, u, v, width, height, p, colors, Y_white)) == DFX_SUCCESS) {

		/* perform brightness adjustment: */
		for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
			/* retrieve L* u* v* values: */
			L_in = L[(p + y) * (width + 2 * p) + p + x];
			u_in = u[(p + y) * (width + 2 * p) + p + x];
			v_in = v[(p + y) * (width + 2 * p) + p + x];
			/* apply gain + saturation to L* channel: */
			L_out = gain_sat(L_in, brightness, 100.f, k);
			/* compute scale: */
			if (L_in > 1e-5) L_scale = L_out / L_in;
			else L_scale = brightness;
			/* propagate scale to u* and v*: */
			u_out = u_in * L_scale;
			v_out = v_in * L_scale;
			/* store the results: */
			L[(p + y) * (width + 2 * p) + p + x] = L_out;
			u[(p + y) * (width + 2 * p) + p + x] = u_out;
			v[(p + y) * (width + 2 * p) + p + x] = v_out;
		}

		/* convert Luv image back to XYZ and then RGB: */
		err = luv_to_xyz(L, u, v, X, Y, Z, width, height, p, colors, Y_white);
		if (err == DFX_SUCCESS)
			err = xyz_to_linear_rgb(X, Y, Z, R_out, G_out, B_out, width, height, p, colors, Y_white);
	}

	/* free intermediate buffers & exit: */
	free_image(X, Y, Z);
	free_image(L, u, v);
	return err;
}

/*!
 *  \brief Saturation adjustment function.
 *
 *  Everything is done in CIE L*u*v* space for the sake of conceptual simplicity. 
 *  For each color, we compute its current saturation and the distance to RGB gamut. 
 *  Then we apply saturation scale and gain_sat() compander, such that the resulting 
 *  value stays in the gamut. The scale is propagated to u* and v*. 
 * 
 *  \param[in]  R_in, G_in, B_in     - input linear RGB image
 *  \param[out] R_out, G_out, B_out  - output RGB image
 *  \param[in]  width                - image width
 *  \param[in]  height               - image height
 *  \param[in]  p                    - padding parameter
 *  \param[in]  saturation           - saturation scale adjustment parameter (0..4)
 *  \param[in]  k                    - clipping shape parameter (0 - hard clipping, 1+ - use generalized harmonic mean compander)
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int adjust_saturation(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float saturation, int k)
{ 
	/* local constants: */
	const int colors = COLORS_SRGB;			/* input/output color space */
	const float Y_white = 1.0f;				/* luminance of white, used in RGB<->XYZ conversions */

	/* temp buffers & local variables: */
	float* X = NULL, * Y = NULL, * Z = NULL;
	float* L = NULL, * u = NULL, * v = NULL;
	float L_in, u_in, v_in, u_out, v_out;
	float s_uv_in, s_uv_out, s_scale, s_max_scale;
	int x, y, err = DFX_SUCCESS;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (saturation < 0 || saturation > 4.0f) return DFX_INVARG;

	/* initialize uv color tables: */
	gen_uv_prim();

	/* allocate intermediate frames: */
	if (alloc_image(&X, &Y, &Z, width, height, p) != DFX_SUCCESS) return DFX_NOMEM;
	if (alloc_image(&L, &u, &v, width, height, p) != DFX_SUCCESS) { free_image(X, Y, Z); return DFX_NOMEM; }

	/* convert RGB image to XYZ and then Luv: */
	if ((err = linear_rgb_to_xyz(R_in, G_in, B_in, X, Y, Z, width, height, p, colors, Y_white)) == DFX_SUCCESS
		&& (err = xyz_to_luv(X, Y, Z, L, u, v, width, height, p, colors, Y_white)) == DFX_SUCCESS) {

		/* perform saturation adjustment: */
		for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
			/* retrieve L* u* v* values: */
			L_in = L[(p + y) * (width + 2 * p) + p + x];
			u_in = u[(p + y) * (width + 2 * p) + p + x];
			v_in = v[(p + y) * (width + 2 * p) + p + x];
			/* compute saturation & gap of this point to gamut: */
			luv_sat_scale(L_in, u_in, v_in, &s_uv_in, &s_max_scale, colors);
			/* apply saturation scale and compander: */
			s_uv_out = gain_sat(s_uv_in, saturation, s_uv_in * s_max_scale, k);
			/* compute the effective scale: */
			if (s_uv_in > 1e-5)  s_scale = s_uv_out / s_uv_in;
			else s_scale = saturation;
			/* apply it to u* and v*: */
			u_out = u_in * s_scale;
			v_out = v_in * s_scale;
			/* store the results: */
			L[(p + y) * (width + 2 * p) + p + x] = L_in;
			u[(p + y) * (width + 2 * p) + p + x] = u_out;
			v[(p + y) * (width + 2 * p) + p + x] = v_out;
		}

		/* convert Luv image back to XYZ and then RGB: */
		err = luv_to_xyz(L, u, v, X, Y, Z, width, height, p, colors, Y_white);
		if (err == DFX_SUCCESS)
			err = xyz_to_linear_rgb(X, Y, Z, R_out, G_out, B_out, width, height, p, colors, Y_white);
	}

	/* free intermediate buffers & exit: */
	free_image(X, Y, Z);
	free_image(L, u, v);
	return err;
}

/*!
 *  \brief Hue adjustment function.
 *
 *  This function performs hue rotation in CIE L*u*v* space.
 *  The colors pushed out of the gamut are forced in by saturation scaling. 
 *
 *  \param[in]  R_in, G_in, B_in     - input linear RGB image
 *  \param[out] R_out, G_out, B_out  - output RGB image
 *  \param[in]  width                - image width
 *  \param[in]  height               - image height
 *  \param[in]  p                    - padding parameter
 *  \param[in]  theta                - hue rotation angle in radians (-pi..pi)
 *
 *  \returns    DFX_SUCCESS if success, error code otherwise
 */
int adjust_hue(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float theta)
{
	/* local constants: */
	const int colors = COLORS_SRGB;			/* assumed input/output color space */
	const float Y_white = 1.0f;				/* assumed luminance of white in RGB<->XYZ conversions */

	/* temp buffers & local variables: */
	float* X = NULL, * Y = NULL, * Z = NULL;
	float* L = NULL, * u = NULL, * v = NULL;
	float L_in, u_in, v_in, u_out, v_out;
	float cos_theta, sin_theta;
	float s_uv_out, s_scale;
	int x, y, err = DFX_SUCCESS;

	/* check parameters: */
	if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
	if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
	if (height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (theta < -(float)M_PI || theta > (float)M_PI) return DFX_INVARG;	/* compare as floats: (float)M_PI > M_PI in double */

	/* initialize uv color tables: */
	gen_uv_prim();

	/* compute rotation factors: */
	cos_theta = cosf(theta);
	sin_theta = sinf(theta);

	/* allocate intermediate frames: */
	if (alloc_image(&X, &Y, &Z, width, height, p) != DFX_SUCCESS) return DFX_NOMEM;
	if (alloc_image(&L, &u, &v, width, height, p) != DFX_SUCCESS) { free_image(X, Y, Z); return DFX_NOMEM; }

	/* convert RGB image to XYZ and then L*u*v*: */
	if ((err = linear_rgb_to_xyz(R_in, G_in, B_in, X, Y, Z, width, height, p, colors, Y_white)) == DFX_SUCCESS
		&& (err = xyz_to_luv(X, Y, Z, L, u, v, width, height, p, colors, Y_white)) == DFX_SUCCESS) {

		/* scan image: */
		for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
			/* retrieve L* u* v* values: */
			L_in = L[(p + y) * (width + 2 * p) + p + x];
			u_in = u[(p + y) * (width + 2 * p) + p + x];
			v_in = v[(p + y) * (width + 2 * p) + p + x];
			/* rotate u*v* by theta: */
			u_out = u_in * cos_theta - v_in * sin_theta;
			v_out = u_in * sin_theta + v_in * cos_theta;
			/* compute saturation & gap of this point to gamut: */
			luv_sat_scale(L_in, u_out, v_out, &s_uv_out, &s_scale, colors);
			/* check if the rotated color is outside of the gamut: */
			if (s_scale < 1.0) {
				/* desaturate color:  */
				u_out = u_out * s_scale;
				v_out = v_out * s_scale;
			}
			/* store the results: */
			L[(p + y) * (width + 2 * p) + p + x] = L_in;
			u[(p + y) * (width + 2 * p) + p + x] = u_out;
			v[(p + y) * (width + 2 * p) + p + x] = v_out;
		}

		/* convert Luv image back to XYZ and then RGB: */
		err = luv_to_xyz(L, u, v, X, Y, Z, width, height, p, colors, Y_white);
		if (err == DFX_SUCCESS)
			err = xyz_to_linear_rgb(X, Y, Z, R_out, G_out, B_out, width, height, p, colors, Y_white);
	}

	/* free intermediate buffers & exit: */
	free_image(X, Y, Z);
	free_image(L, u, v);
	return err;
}

/* dfx_colors.c -- end of file */
