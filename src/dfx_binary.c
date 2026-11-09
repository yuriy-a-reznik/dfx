/*!
 *  \file   dfx_binary.c
 *  \brief  Binary images: thresholding and morphological operations.
 *
 *  A binary image is a set of pixel coordinates: the "ON" (foreground) pixels. We store it as a
 *  plane of unsigned chars holding the values 0 (OFF) and 1 (ON), with exactly the same geometry
 *  as the float planes of linear images: a plane for a width x height image with padding p holds
 *  (width + 2p) x (height + 2p) values, and the image pixel (x,y) is stored at index
 *  (p + y) * (width + 2p) + (p + x). The padding plays the same role as for linear filtering:
 *  it supplies the out-of-image neighbors needed by a structuring element of radius n <= p.
 *
 *  Binary images are deliberately kept apart from linear planes. The operations of mathematical
 *  morphology (dilation, erosion, hit-or-miss, thinning, ...) are set operations, that is, logical
 *  operations on memberships, and have nothing to do with light; the two worlds meet only in the
 *  thresholding functions (plane -> binary image) and in binary_to_plane() (binary image -> plane
 *  of 0.0 and 1.0 values, which can be displayed through luminance_to_grayscale_image()).
 *
 *  Structuring elements (SEs) are small binary arrays of (2n+1) x (2n+1) values with the origin
 *  at the center: B[(j + n) * (2n+1) + (i + n)] is the value at offset (i,j). Following
 *  S. Birchfield, Image Processing and Analysis (Cengage, 2017), chapter 4, we define
 *
 *     dilation:  (A (+) B)(z) = OR  over b in B of A(z - b),
 *     erosion:   (A (-) B)(z) = AND over b in B of A(z + b),
 *
 *  so that erosion with B undoes dilation with B as far as possible (erosion uses the reflected SE,
 *  which makes no difference for the symmetric SEs produced by structuring_element()). The
 *  hit-or-miss operator takes a ternary SE whose values are SE_ON (the pixel must be ON), SE_OFF
 *  (the pixel must be OFF), and SE_DONTCARE; it is the erosion of A with the ON pattern intersected
 *  with the erosion of the complement of A with the OFF pattern.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.2
 *  \date    November 8, 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dfx.h"

/**************
 *
 *  Allocation, initialization, padding, and copying of binary images:
 *
 *   binary_size()
 *   alloc_binary()
 *   free_binary()
 *   fill_binary()
 *   pad_binary()
 *   copy_binary()
 */

/*!
 *  \brief Compute the size of a padded binary image (in pixels).
 */
static size_t binary_size(int width, int height, int p)
{
	return (size_t)(2 * p + height) * (size_t)(2 * p + width);
}

/*!
 *  \brief Allocate a binary image.
 *
 *  \param[out]     pA     - pointer to a variable to receive the pointer to the allocated image
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *
 *  \returns        DFX_X error code
 */
int alloc_binary(unsigned char** pA, int width, int height, int p)
{
	size_t size;
	if (pA == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	if ((*pA = (unsigned char*)malloc(size > 0 ? size : 1)) == NULL) return DFX_NOMEM;	/* malloc(0) may return NULL */
	return DFX_SUCCESS;
}

/*!
 *  \brief Free a binary image.
 */
int free_binary(unsigned char* A)
{
	if (A != NULL) free(A);
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill a binary image (padding included) with a constant: 0 or 1.
 */
int fill_binary(unsigned char* A, int width, int height, int p, int value)
{
	if (A == NULL || height < 0 || width < 0 || p < 0 || (value != 0 && value != 1)) return DFX_INVARG;
	memset((void*)A, value, binary_size(width, height, p));
	return DFX_SUCCESS;
}

/*!
 *  \brief Fill the padding border of a binary image.
 *
 *  Same rules as pad_plane(): PAD_ZERO fills the border with OFF pixels (the outside of the image
 *  is background), PAD_REPLICATE repeats the boundary pixels (the image continues outside),
 *  PAD_REFLECT mirrors the image about its boundary (requires p <= width, height).
 *
 *  \param[in,out]  A            - binary image; its interior is left unchanged
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter
 *  \param[in]      padding_type - PAD_ZERO, PAD_REPLICATE, or PAD_REFLECT
 */
int pad_binary(unsigned char* A, int width, int height, int p, int padding_type)
{
	int width_p = width + 2 * p;	/* width of the padded image */
	int height_p = height + 2 * p;	/* height of the padded image */
	int x, y, xs, ys;

	/* check parameters: */
	if (A == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (padding_type < 0 || padding_type >= N_PADDINGS) return DFX_INVARG;
	if (padding_type == PAD_REFLECT && (p > width || p > height)) return DFX_INVARG;
	if (p == 0) return DFX_SUCCESS;

	if (width == 0 || height == 0) {
		/* nothing to replicate or reflect: */
		memset((void*)A, 0, binary_size(width, height, p));
		return DFX_SUCCESS;
	}

	/* fill the border, pixel by pixel: */
	for (y = 0; y < height_p; y++) for (x = 0; x < width_p; x++) {
		/* skip the interior: */
		if (y >= p && y < p + height && x >= p && x < p + width) continue;
		switch (padding_type) {
		case PAD_ZERO:
			A[y * width_p + x] = 0;
			break;
		case PAD_REPLICATE:
			/* coordinates of the nearest boundary pixel: */
			xs = x < p ? p : (x >= p + width ? p + width - 1 : x);
			ys = y < p ? p : (y >= p + height ? p + height - 1 : y);
			A[y * width_p + x] = A[ys * width_p + xs];
			break;
		default: /* PAD_REFLECT: the boundary pixel is repeated */
			xs = x < p ? 2 * p - 1 - x : (x >= p + width ? 2 * (p + width) - 1 - x : x);
			ys = y < p ? 2 * p - 1 - y : (y >= p + height ? 2 * (p + height) - 1 - y : y);
			A[y * width_p + x] = A[ys * width_p + xs];
			break;
		}
	}
	return DFX_SUCCESS;
}

/*!
 *  \brief Copy a binary image (padding included).
 */
int copy_binary(const unsigned char* A_in, unsigned char* A_out, int width, int height, int p)
{
	if (A_in == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (A_in != A_out)
		memmove((void*)A_out, (const void*)A_in, binary_size(width, height, p));
	return DFX_SUCCESS;
}

/**************
 *
 *  Set (logical) operations:
 *
 *   not_binary()
 *   and_binary()
 *   or_binary()
 *   xor_binary()
 *   subtract_binary()
 *   count_binary()
 *
 *  These operations process the padded area as well, so the padding of their results stays valid.
 *  The output may be the same as one of the inputs.
 */

/*!
 *  \brief Complement of a binary image: A_out = NOT A.
 */
int not_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A_out[i] = (unsigned char)(A[i] ? 0 : 1);
	return DFX_SUCCESS;
}

/*!
 *  \brief Intersection of two binary images: A_out = A AND B.
 */
int and_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || B == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A_out[i] = (unsigned char)((A[i] && B[i]) ? 1 : 0);
	return DFX_SUCCESS;
}

/*!
 *  \brief Union of two binary images: A_out = A OR B.
 */
int or_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || B == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A_out[i] = (unsigned char)((A[i] || B[i]) ? 1 : 0);
	return DFX_SUCCESS;
}

/*!
 *  \brief Symmetric difference of two binary images: A_out = A XOR B.
 */
int xor_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || B == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A_out[i] = (unsigned char)(((A[i] != 0) != (B[i] != 0)) ? 1 : 0);
	return DFX_SUCCESS;
}

/*!
 *  \brief Set difference of two binary images: A_out = A AND NOT B.
 */
int subtract_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || B == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A_out[i] = (unsigned char)((A[i] && !B[i]) ? 1 : 0);
	return DFX_SUCCESS;
}

/*!
 *  \brief Count the ON pixels of a binary image (the area of the foreground). Only the interior is counted.
 */
int count_binary(const unsigned char* A, int width, int height, int p, int* p_count)
{
	int width_p = width + 2 * p;
	int x, y, count = 0;
	if (A == NULL || p_count == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	for (y = 0; y < height; y++) for (x = 0; x < width; x++)
		if (A[(p + y) * width_p + p + x]) count++;
	*p_count = count;
	return DFX_SUCCESS;
}

/**************
 *
 *  Conversions between planes and binary images:
 *
 *   binary_to_plane()
 *   threshold_plane()
 *   histogram_plane()
 *   otsu_threshold()
 *   ridler_calvard_threshold()
 *   niblack_threshold_plane()
 *   double_threshold_plane()
 */

/*!
 *  \brief Convert a binary image to a plane of 0.0 (OFF) and 1.0 (ON) values, padding included.
 *
 *  Interpreted as luminance, such a plane renders the binary image as black and white
 *  (see luminance_to_grayscale_image() and linear_to_srgb()).
 */
int binary_to_plane(const unsigned char* A, float* X, int width, int height, int p)
{
	size_t i, size;
	if (A == NULL || X == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) X[i] = A[i] ? 1.0f : 0.0f;
	return DFX_SUCCESS;
}

/*!
 *  \brief Threshold a plane: A = (X > tau), padding included.
 *
 *  This is the only place (together with the adaptive and double thresholding functions below) where
 *  a quantity measured in linear space (luminance, a color component, a distance, a gradient magnitude...)
 *  turns into a set. The padding of X is thresholded as well, so the border of A follows the same
 *  convention as the border of X (e.g. replicate padding stays replicate padding).
 *
 *  \param[in]      X      - input plane
 *  \param[out]     A      - output binary image
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both
 *  \param[in]      tau    - threshold: pixels with values strictly greater than tau are set ON
 */
int threshold_plane(const float* X, unsigned char* A, int width, int height, int p, float tau)
{
	size_t i, size;
	if (X == NULL || A == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = binary_size(width, height, p);
	for (i = 0; i < size; i++) A[i] = (unsigned char)(X[i] > tau ? 1 : 0);
	return DFX_SUCCESS;
}

/*!
 *  \brief Compute the histogram of the values of a plane (interior only).
 *
 *  The range [lo, hi] is divided into nbins equal bins; bin k covers [lo + k*d, lo + (k+1)*d) with
 *  d = (hi - lo) / nbins. Values below lo (and NaNs) are counted in the first bin and values at or
 *  above hi in the last one, so that the counts always sum to width * height.
 *
 *  \param[in]      X      - input plane
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *  \param[in]      lo, hi - range of values covered by the histogram (lo < hi)
 *  \param[in]      nbins  - number of bins (>= 1)
 *  \param[out]     hist   - array of nbins counts
 */
int histogram_plane(const float* X, int width, int height, int p, float lo, float hi, int nbins, int* hist)
{
	int width_p = width + 2 * p;
	int x, y, k;
	double scale;

	/* check parameters: */
	if (X == NULL || hist == NULL || height < 0 || width < 0 || p < 0 || nbins < 1 || !(lo < hi)) return DFX_INVARG;

	/* clear and accumulate: */
	for (k = 0; k < nbins; k++) hist[k] = 0;
	scale = (double)nbins / ((double)hi - (double)lo);
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		double v = ((double)X[(p + y) * width_p + p + x] - (double)lo) * scale;
		if (!(v >= 0)) k = 0;							/* also catches NaN */
		else if (v >= (double)nbins) k = nbins - 1;
		else k = (int)v;
		hist[k]++;
	}
	return DFX_SUCCESS;
}

/*!
 *  \brief The threshold value separating levels 0..l from levels l+1..nbins-1.
 *
 *  Level l+1 starts at the bin edge lo + (l+1) d; since threshold_plane() turns the pixels with X > tau ON,
 *  tau is set to the largest float below that edge, so that exactly the pixels of levels l+1 and above are
 *  selected, also when the plane holds integer values that fall on bin edges. l may be fractional (the
 *  middle of a plateau of equivalent thresholds).
 */
static float level_threshold(double l, float lo, float hi, int nbins)
{
	float edge = (float)((double)lo + (l + 1.0) * ((double)hi - (double)lo) / nbins);
	return nextafterf(edge, -HUGE_VALF);
}

/*!
 *  \brief Compute the cumulative histogram moments m0[k] = sum_{l <= k} h[l] and m1[k] = sum_{l <= k} l * h[l].
 */
static int cumulative_moments(const float* X, int width, int height, int p, float lo, float hi, int nbins, double** p_m0, double** p_m1)
{
	int* hist;
	double* m0, * m1;
	int k, err;

	if ((hist = (int*)malloc((size_t)nbins * sizeof(int))) == NULL) return DFX_NOMEM;
	if ((err = histogram_plane(X, width, height, p, lo, hi, nbins, hist)) != DFX_SUCCESS) { free(hist); return err; }
	if ((m0 = (double*)malloc(2 * (size_t)nbins * sizeof(double))) == NULL) { free(hist); return DFX_NOMEM; }
	m1 = m0 + nbins;
	m0[0] = hist[0]; m1[0] = 0;
	for (k = 1; k < nbins; k++) {
		m0[k] = m0[k - 1] + hist[k];
		m1[k] = m1[k - 1] + (double)k * hist[k];
	}
	free(hist);
	*p_m0 = m0; *p_m1 = m1;
	return DFX_SUCCESS;
}

/*!
 *  \brief Find a global threshold by Otsu's method.
 *
 *  The values of the plane are quantized into nbins levels over [lo, hi] (see histogram_plane()), and
 *  the level l that splits the histogram into two classes with the maximum between-class variance
 *  (equivalently, the minimum within-class variance) is found:
 *
 *      sigma_b^2(l) = (m1[l] - mu * m0[l])^2 / (m0[l] * (N - m0[l])),
 *
 *  where m0[l] and m1[l] are the cumulative count and the cumulative sum of levels up to l, N is the
 *  number of pixels, and mu = m1[nbins-1] / N is the mean level (Birchfield, Algorithm 10.2). The
 *  returned threshold tau separates the two classes: threshold_plane(X, A, ..., tau) turns exactly the
 *  pixels of levels above l ON (tau is the largest value below the lower edge of level l+1). When several
 *  consecutive levels give the same maximum (which happens when no pixel falls between the two classes),
 *  the middle of that range is used. If all pixels fall into one level, tau is the upper edge of that level.
 *
 *  Reference: N. Otsu, "A threshold selection method from gray-level histograms," IEEE Trans. Syst.
 *  Man Cybern., vol. 9, no. 1, pp. 62-66, 1979.
 *
 *  \param[in]      X      - input plane
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *  \param[in]      lo, hi - range of values covered by the histogram, e.g. 0 and 1 for luminance (lo < hi)
 *  \param[in]      nbins  - number of histogram bins (>= 2), e.g. 256
 *  \param[out]     p_tau  - receives the threshold
 */
int otsu_threshold(const float* X, int width, int height, int p, float lo, float hi, int nbins, float* p_tau)
{
	double* m0, * m1, N, mu, s, best = -1.0;
	int l, l_first = -1, l_last = -1, err;

	/* check parameters: */
	if (X == NULL || p_tau == NULL || height <= 0 || width <= 0 || p < 0 || nbins < 2 || !(lo < hi)) return DFX_INVARG;

	/* cumulative moments of the histogram: */
	if ((err = cumulative_moments(X, width, height, p, lo, hi, nbins, &m0, &m1)) != DFX_SUCCESS) return err;
	N = m0[nbins - 1];
	mu = m1[nbins - 1] / N;

	/* maximize the between-class variance: */
	for (l = 0; l < nbins - 1; l++) {
		double w0 = m0[l], w1 = N - m0[l];
		if (w0 <= 0 || w1 <= 0) continue;				/* one of the classes is empty */
		s = m1[l] - mu * w0;
		s = s * s / (w0 * w1);
		if (s > best) { best = s; l_first = l_last = l; }
		else if (s == best) l_last = l;					/* a plateau of equivalent thresholds (empty levels between the classes) */
	}
	if (l_first < 0) {
		/* all pixels in one level: find it */
		for (l = 0; l < nbins && m0[l] <= 0; l++);
		l_first = l_last = l < nbins ? l : nbins - 1;
	}
	free(m0);

	/* the threshold is the class boundary; on a plateau, the one in the middle: */
	*p_tau = level_threshold(0.5 * (l_first + l_last), lo, hi, nbins);
	return DFX_SUCCESS;
}

/*!
 *  \brief Find a global threshold by the Ridler-Calvard iterative selection method.
 *
 *  The values of the plane are quantized into nbins levels over [lo, hi] (see histogram_plane()).
 *  Starting from the middle level, the threshold is repeatedly replaced by the average of the mean
 *  level of the pixels at or below it and the mean level of the pixels above it, until it stops
 *  changing (Birchfield, Algorithm 10.1). The average is rounded down: the levels at or below it form
 *  the lower class, exactly as in the continuous version of the algorithm (Algorithm 10.1 rounds to the
 *  nearest level, which shifts the result up by half a level on average). The returned tau separates the
 *  two classes, in the same sense as for otsu_threshold().
 *
 *  Reference: T. W. Ridler and S. Calvard, "Picture thresholding using an iterative selection method,"
 *  IEEE Trans. Syst. Man Cybern., vol. SMC-8, no. 8, pp. 630-632, 1978.
 *
 *  Parameters are the same as for otsu_threshold().
 */
int ridler_calvard_threshold(const float* X, int width, int height, int p, float lo, float hi, int nbins, float* p_tau)
{
	double* m0, * m1, N, mu_lo, mu_hi;
	int l, l_new, l_min, l_max, iter, err;

	/* check parameters: */
	if (X == NULL || p_tau == NULL || height <= 0 || width <= 0 || p < 0 || nbins < 2 || !(lo < hi)) return DFX_INVARG;

	/* cumulative moments of the histogram: */
	if ((err = cumulative_moments(X, width, height, p, lo, hi, nbins, &m0, &m1)) != DFX_SUCCESS) return err;
	N = m0[nbins - 1];

	/* occupied range of levels: thresholds in [l_min, l_max - 1] leave both classes nonempty */
	for (l_min = 0; l_min < nbins - 1 && m0[l_min] <= 0; l_min++);
	for (l_max = nbins - 1; l_max > l_min && m0[l_max - 1] >= N; l_max--);

	/* iterate (if all pixels fall into one level, tau becomes the upper edge of that level, as in otsu_threshold()): */
	l = nbins / 2;
	if (l < l_min) l = l_min;
	if (l > l_max - 1) l = l_max - 1;
	if (l_min >= l_max) l = l_min;
	for (iter = 0; l_min < l_max && iter < 1000; iter++) {
		mu_lo = m1[l] / m0[l];										/* mean level of the pixels at or below l */
		mu_hi = (m1[nbins - 1] - m1[l]) / (N - m0[l]);				/* mean level of the pixels above l */
		l_new = (int)floor(0.5 * (mu_lo + mu_hi));					/* levels at or below the mean of the means stay in the lower class */
		if (l_new < l_min) l_new = l_min;
		if (l_new > l_max - 1) l_new = l_max - 1;
		if (l_new == l) break;
		l = l_new;
	}
	free(m0);

	*p_tau = level_threshold((double)l, lo, hi, nbins);
	return DFX_SUCCESS;
}

/*!
 *  \brief Adaptive (local) thresholding by Niblack's method.
 *
 *  Each pixel is compared with a threshold computed from the statistics of its (2n+1) x (2n+1)
 *  neighborhood:
 *
 *      tau(x,y) = mu(x,y) - k * sigma(x,y),     A(x,y) = (X(x,y) > tau(x,y)),
 *
 *  where mu and sigma are the local mean and standard deviation. With k > 0 the threshold is lowered
 *  below the local mean; k = 0 compares each pixel with its local mean; k < 0 raises the threshold.
 *  (Niblack's original formulation, tau = mu + k' * sigma with k' = -0.2, is the same formula with
 *  k = -k'.) The input must be padded with p >= n; the output receives PAD_REPLICATE padding.
 *
 *  Reference: W. Niblack, An Introduction to Digital Image Processing, Prentice-Hall, 1986, pp. 115-116.
 *
 *  \param[in]      X      - input plane, padded
 *  \param[out]     A      - output binary image
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both
 *  \param[in]      n      - neighborhood radius: 1 <= n <= p
 *  \param[in]      k      - weight of the local standard deviation (about 0.2)
 */
int niblack_threshold_plane(const float* X, unsigned char* A, int width, int height, int p, int n, float k)
{
	int width_p = width + 2 * p;
	int x, y, i, j;
	double count, sum, sum2, mu, var, tau;

	/* check parameters: */
	if (X == NULL || A == NULL || height < 0 || width < 0 || p < 0 || n < 1 || n > p) return DFX_INVARG;

	count = (double)(2 * n + 1) * (2 * n + 1);
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		/* local mean and variance: */
		sum = sum2 = 0;
		for (j = -n; j <= n; j++) for (i = -n; i <= n; i++) {
			double v = X[(p + y + j) * width_p + p + x + i];
			sum += v; sum2 += v * v;
		}
		mu = sum / count;
		var = sum2 / count - mu * mu;
		if (var < 0) var = 0;
		tau = mu - (double)k * sqrt(var);
		A[(p + y) * width_p + p + x] = (unsigned char)((double)X[(p + y) * width_p + p + x] > tau ? 1 : 0);
	}
	return pad_binary(A, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Double thresholding (thresholding with hysteresis).
 *
 *  Pixels above the high threshold are certainly foreground; pixels above the low threshold are
 *  foreground only if they are connected (through pixels above the low threshold) to a pixel above the
 *  high threshold. The result is the set of connected components of (X > tau_lo) that contain at least
 *  one pixel of (X > tau_hi), computed with reconstruct_binary(). The output receives PAD_REPLICATE padding.
 *
 *  \param[in]      X            - input plane
 *  \param[out]     A            - output binary image
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter of both
 *  \param[in]      tau_lo       - low threshold
 *  \param[in]      tau_hi       - high threshold (>= tau_lo)
 *  \param[in]      connectivity - 4 or 8
 */
int double_threshold_plane(const float* X, unsigned char* A, int width, int height, int p, float tau_lo, float tau_hi, int connectivity)
{
	unsigned char* mask = NULL, * seeds = NULL;
	int err;

	/* check parameters: */
	if (X == NULL || A == NULL || height < 0 || width < 0 || p < 0 || !(tau_lo <= tau_hi)) return DFX_INVARG;
	if (connectivity != 4 && connectivity != 8) return DFX_INVARG;

	/* threshold twice and keep the components of the low-threshold image seeded by the high-threshold one: */
	if ((err = alloc_binary(&mask, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = alloc_binary(&seeds, width, height, p)) != DFX_SUCCESS) { free_binary(mask); return err; }
	threshold_plane(X, mask, width, height, p, tau_lo);
	threshold_plane(X, seeds, width, height, p, tau_hi);
	err = reconstruct_binary(mask, seeds, A, width, height, p, connectivity);
	free_binary(mask);
	free_binary(seeds);
	return err;
}

/**************
 *
 *  Morphological operations:
 *
 *   structuring_element()
 *   dilate_binary()
 *   erode_binary()
 *   open_binary()
 *   close_binary()
 *   majority_binary()
 *   hitmiss_binary()
 *   thin_binary()
 *   thicken_binary()
 *   boundary_binary()
 *
 *  The input image must be padded with p >= n, where n is the radius of the structuring element.
 *  The values in the padding decide what lies outside the image: pad with PAD_ZERO to treat the
 *  outside as background, or with PAD_REPLICATE to let objects touching the image boundary continue
 *  beyond it (for the symmetric structuring elements this makes dilation ignore the outside and
 *  erosion leave the boundary side of such objects alone). The output receives PAD_REPLICATE padding.
 *  The output image may be the same as the input.
 */

/*!
 *  \brief Generate a structuring element of radius n.
 *
 *  \param[out]     B      - array of (2n+1) x (2n+1) values to receive the structuring element
 *  \param[in]      n      - radius (>= 0); B[(j + n) * (2n+1) + (i + n)] is the value at offset (i,j)
 *  \param[in]      shape  - SE_CROSS: |i| + |j| <= n (the "city-block ball"; for n = 1 the 4-neighborhood B4),
 *                           SE_BOX:   max(|i|,|j|) <= n (the "chessboard ball"; for n = 1 the 8-neighborhood B8),
 *                           SE_DISK:  i^2 + j^2 <= n^2 (the discrete Euclidean ball).
 */
int structuring_element(unsigned char* B, int n, int shape)
{
	int i, j, m = 2 * n + 1, in;
	if (B == NULL || n < 0 || shape < 0 || shape >= N_SES) return DFX_INVARG;
	for (j = -n; j <= n; j++) for (i = -n; i <= n; i++) {
		switch (shape) {
		case SE_CROSS: in = abs(i) + abs(j) <= n; break;
		case SE_BOX:   in = 1; break;
		default:       in = i * i + j * j <= n * n; break;
		}
		B[(j + n) * m + (i + n)] = (unsigned char)in;
	}
	return DFX_SUCCESS;
}

/*!
 *  \brief Check a structuring element: values must be 0/1 (ternary = 0 allows SE_DONTCARE as well).
 */
static int check_se(const unsigned char* B, int n, int ternary)
{
	int i, m = (2 * n + 1) * (2 * n + 1);
	for (i = 0; i < m; i++)
		if (B[i] > (ternary ? SE_DONTCARE : SE_ON)) return DFX_INVARG;
	return DFX_SUCCESS;
}

/*!
 *  \brief Make a private copy of the input when an operation is requested in place.
 *
 *  Returns the pointer to use as the input, and sets *p_tmp to the buffer to free afterwards (or NULL).
 */
static const unsigned char* input_copy(const unsigned char* A, unsigned char* A_out, int width, int height, int p, unsigned char** p_tmp)
{
	*p_tmp = NULL;
	if (A != A_out) return A;
	if (alloc_binary(p_tmp, width, height, p) != DFX_SUCCESS) return NULL;
	copy_binary(A, *p_tmp, width, height, p);
	return *p_tmp;
}

/*!
 *  \brief Dilate a binary image: A_out(z) = OR over b in B of A(z - b).
 *
 *  \param[in]      A      - input binary image, padded with p >= n
 *  \param[out]     A_out  - output binary image (may be the same as A)
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both
 *  \param[in]      B      - structuring element of (2n+1) x (2n+1) values 0/1
 *  \param[in]      n      - radius of the structuring element: 0 <= n <= p
 */
int dilate_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n)
{
	int width_p = width + 2 * p, m = 2 * n + 1;
	int x, y, i, j, v;
	unsigned char* tmp;
	const unsigned char* in;

	/* check parameters: */
	if (A == NULL || A_out == NULL || B == NULL || height < 0 || width < 0 || p < 0 || n < 0 || n > p) return DFX_INVARG;
	if (check_se(B, n, 0) != DFX_SUCCESS) return DFX_INVARG;
	if ((in = input_copy(A, A_out, width, height, p, &tmp)) == NULL) return DFX_NOMEM;

	/* dilate: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		v = 0;
		for (j = -n; j <= n && !v; j++) for (i = -n; i <= n; i++)
			if (B[(j + n) * m + (i + n)] && in[(p + y - j) * width_p + p + x - i]) { v = 1; break; }
		A_out[(p + y) * width_p + p + x] = (unsigned char)v;
	}
	free_binary(tmp);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Erode a binary image: A_out(z) = AND over b in B of A(z + b). Same parameters as dilate_binary().
 */
int erode_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n)
{
	int width_p = width + 2 * p, m = 2 * n + 1;
	int x, y, i, j, v;
	unsigned char* tmp;
	const unsigned char* in;

	/* check parameters: */
	if (A == NULL || A_out == NULL || B == NULL || height < 0 || width < 0 || p < 0 || n < 0 || n > p) return DFX_INVARG;
	if (check_se(B, n, 0) != DFX_SUCCESS) return DFX_INVARG;
	if ((in = input_copy(A, A_out, width, height, p, &tmp)) == NULL) return DFX_NOMEM;

	/* erode: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		v = 1;
		for (j = -n; j <= n && v; j++) for (i = -n; i <= n; i++)
			if (B[(j + n) * m + (i + n)] && !in[(p + y + j) * width_p + p + x + i]) { v = 0; break; }
		A_out[(p + y) * width_p + p + x] = (unsigned char)v;
	}
	free_binary(tmp);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Morphological opening: erosion followed by dilation with the same structuring element.
 *
 *  Removes foreground details smaller than the structuring element (noise specks, thin protrusions).
 *  Same parameters as dilate_binary().
 */
int open_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n)
{
	int err;
	if ((err = erode_binary(A, A_out, width, height, p, B, n)) != DFX_SUCCESS) return err;
	return dilate_binary(A_out, A_out, width, height, p, B, n);
}

/*!
 *  \brief Morphological closing: dilation followed by erosion with the same structuring element.
 *
 *  Fills background details smaller than the structuring element (holes, gaps, notches).
 *  Same parameters as dilate_binary().
 */
int close_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n)
{
	int err;
	if ((err = dilate_binary(A, A_out, width, height, p, B, n)) != DFX_SUCCESS) return err;
	return erode_binary(A_out, A_out, width, height, p, B, n);
}

/*!
 *  \brief Majority filter: A_out(x,y) = 1 if more than half of the pixels in the (2n+1) x (2n+1)
 *  neighborhood of (x,y) are ON.
 *
 *  This is the median filter of a binary image: it removes isolated pixels and small specks of both
 *  colors while preserving edges. The input must be padded with p >= n; the output receives
 *  PAD_REPLICATE padding and may be the same as the input.
 */
int majority_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n)
{
	int width_p = width + 2 * p, half = (2 * n + 1) * (2 * n + 1) / 2;
	int x, y, i, j, count;
	unsigned char* tmp;
	const unsigned char* in;

	/* check parameters: */
	if (A == NULL || A_out == NULL || height < 0 || width < 0 || p < 0 || n < 0 || n > p) return DFX_INVARG;
	if ((in = input_copy(A, A_out, width, height, p, &tmp)) == NULL) return DFX_NOMEM;

	/* count ON neighbors and vote: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		count = 0;
		for (j = -n; j <= n; j++) for (i = -n; i <= n; i++)
			count += in[(p + y + j) * width_p + p + x + i] ? 1 : 0;
		A_out[(p + y) * width_p + p + x] = (unsigned char)(count > half ? 1 : 0);
	}
	free_binary(tmp);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Hit-or-miss transform: find all the places where a ternary pattern matches the image.
 *
 *  A_out(z) = 1 iff A(z + b) = 1 for every b with B(b) = SE_ON and A(z + b) = 0 for every b with
 *  B(b) = SE_OFF; positions with B(b) = SE_DONTCARE are ignored. This equals the intersection of the
 *  erosion of A with the ON pattern and the erosion of NOT A with the OFF pattern (Birchfield, Fig. 4.12).
 *  The input must be padded with p >= n; the output receives PAD_REPLICATE padding and may be the same
 *  as the input.
 *
 *  \param[in]      A      - input binary image, padded
 *  \param[out]     A_out  - output binary image
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both
 *  \param[in]      B      - ternary structuring element of (2n+1) x (2n+1) values SE_OFF, SE_ON, SE_DONTCARE
 *  \param[in]      n      - radius of the structuring element: 0 <= n <= p
 */
int hitmiss_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n)
{
	int width_p = width + 2 * p, m = 2 * n + 1;
	int x, y, i, j, v, b;
	unsigned char* tmp;
	const unsigned char* in;

	/* check parameters: */
	if (A == NULL || A_out == NULL || B == NULL || height < 0 || width < 0 || p < 0 || n < 0 || n > p) return DFX_INVARG;
	if (check_se(B, n, 1) != DFX_SUCCESS) return DFX_INVARG;
	if ((in = input_copy(A, A_out, width, height, p, &tmp)) == NULL) return DFX_NOMEM;

	/* match the pattern at every pixel: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		v = 1;
		for (j = -n; j <= n && v; j++) for (i = -n; i <= n; i++) {
			b = B[(j + n) * m + (i + n)];
			if (b == SE_DONTCARE) continue;
			if ((in[(p + y + j) * width_p + p + x + i] != 0) != (b == SE_ON)) { v = 0; break; }
		}
		A_out[(p + y) * width_p + p + x] = (unsigned char)v;
	}
	free_binary(tmp);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/*
 *  Structuring elements for thinning (Birchfield, Fig. 4.14): four "edge" SEs followed by four "corner"
 *  SEs. Each matches a foreground pixel that can be removed without breaking the connectivity of the
 *  region. Notation below: 0 = SE_OFF, 1 = SE_ON, x = SE_DONTCARE (= 2), rows from top to bottom.
 */
#define x_ SE_DONTCARE
static const unsigned char thinning_se[8][9] = {
	{ 0,  0,  0,    x_, 1, x_,    1,  1,  1  },		/* edge, removes the top side */
	{ 1,  x_, 0,    1,  1, 0,     1,  x_, 0  },		/* edge, removes the right side */
	{ 1,  1,  1,    x_, 1, x_,    0,  0,  0  },		/* edge, removes the bottom side */
	{ 0,  x_, 1,    0,  1, 1,     0,  x_, 1  },		/* edge, removes the left side */
	{ x_, 0,  0,    1,  1, 0,     x_, 1,  x_ },		/* corner, top-right */
	{ x_, 1,  x_,   1,  1, 0,     x_, 0,  0  },		/* corner, bottom-right */
	{ x_, 1,  x_,   0,  1, 1,     0,  0,  x_ },		/* corner, bottom-left */
	{ 0,  0,  x_,   0,  1, 1,     x_, 1,  x_ },		/* corner, top-left */
};

/*
 *  Structuring elements for thickening (Birchfield, Fig. 4.17). Each matches a background pixel
 *  that lies in a concavity of the foreground and can be added to it; applied until convergence,
 *  thickening approximates the convex hull of the regions.
 */
static const unsigned char thickening_se[8][9] = {
	{ 1,  1,  x_,   1,  0,  x_,   1,  x_, 0  },
	{ 1,  1,  1,    x_, 0,  1,    0,  x_, x_ },
	{ 0,  x_, 1,    x_, 0,  1,    x_, 1,  1  },
	{ x_, x_, 0,    1,  0,  x_,   1,  1,  1  },
	{ x_, 1,  1,    x_, 0,  1,    0,  x_, 1  },
	{ 0,  x_, x_,   x_, 0,  1,    1,  1,  1  },
	{ 1,  x_, 0,    1,  0,  x_,   1,  1,  x_ },
	{ 1,  1,  1,    1,  0,  x_,   x_, x_, 0  },
};
#undef x_

/*!
 *  \brief Thinning or thickening by a sequence of hit-or-miss operations (shared implementation).
 *
 *  \param[in]      thicken - 0: A = A - (A (*) B_k) for each SE; 1: A = A + (A (*) B_k) for each SE
 */
static int thin_or_thicken(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n_iter, int* p_iter, int thicken)
{
	const unsigned char(*se)[9] = thicken ? thickening_se : thinning_se;
	unsigned char* H = NULL;
	int iter, k, changed, count, err;

	/* check parameters: */
	if (A == NULL || A_out == NULL || height < 0 || width < 0 || p < 1) return DFX_INVARG;

	/* start from a copy of the input: */
	if ((err = copy_binary(A, A_out, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = alloc_binary(&H, width, height, p)) != DFX_SUCCESS) return err;

	/* iterate the sequence of 8 structuring elements until nothing changes (or n_iter passes are done): */
	for (iter = 0; n_iter <= 0 || iter < n_iter; ) {
		changed = 0;
		for (k = 0; k < 8; k++) {
			if ((err = hitmiss_binary(A_out, H, width, height, p, se[k], 1)) != DFX_SUCCESS) break;
			if ((err = count_binary(H, width, height, p, &count)) != DFX_SUCCESS) break;
			if (count == 0) continue;
			changed = 1;
			pad_binary(H, width, height, p, PAD_ZERO);				/* change the interior only: the padding of the input stays in force */
			if (thicken) err = or_binary(A_out, H, A_out, width, height, p);
			else err = subtract_binary(A_out, H, A_out, width, height, p);
			if (err != DFX_SUCCESS) break;
		}
		if (err != DFX_SUCCESS) break;
		iter++;
		if (!changed) break;
	}
	free_binary(H);
	if (p_iter != NULL) *p_iter = iter;
	if (err == DFX_SUCCESS) err = pad_binary(A_out, width, height, p, PAD_REPLICATE);
	return err;
}

/*!
 *  \brief Morphological thinning.
 *
 *  Repeatedly removes the pixels matched by the eight thinning structuring elements of Birchfield's
 *  Fig. 4.14, applied one after another: A <- A - (A (*) B_k), k = 1..8. One pass over the eight SEs
 *  peels one layer of pixels off every side of every region while preserving connectivity; iterated
 *  to convergence, thinning produces an approximate skeleton of the foreground. The input must be
 *  padded with p >= 1; its border decides what lies outside the image for the whole iteration (see the
 *  module description). The output receives PAD_REPLICATE padding and may be the same as the input.
 *
 *  \param[in]      A      - input binary image, padded
 *  \param[out]     A_out  - output binary image
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both (>= 1)
 *  \param[in]      n_iter - number of passes over the eight structuring elements; <= 0 to iterate until convergence
 *  \param[out]     p_iter - if not NULL, receives the number of passes performed
 */
int thin_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n_iter, int* p_iter)
{
	return thin_or_thicken(A, A_out, width, height, p, n_iter, p_iter, 0);
}

/*!
 *  \brief Morphological thickening.
 *
 *  Repeatedly adds the background pixels matched by the eight thickening structuring elements of
 *  Birchfield's Fig. 4.17: A <- A + (A (*) B_k), k = 1..8. Iterated to convergence, thickening
 *  approximates the convex hull of each region. Same parameters as thin_binary().
 */
int thicken_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n_iter, int* p_iter)
{
	return thin_or_thicken(A, A_out, width, height, p, n_iter, p_iter, 1);
}

/*!
 *  \brief Extract the (inner) boundary of the foreground: A_out = A - (A eroded by a 3 x 3 structuring element).
 *
 *  With connectivity = 8, the erosion uses the 4-neighborhood B4 and the result is the 8-connected
 *  boundary: the ON pixels having an OFF 4-neighbor. With connectivity = 4, the erosion uses the
 *  8-neighborhood B8 and the result is the thicker 4-connected boundary: the ON pixels having an OFF
 *  8-neighbor (Birchfield, Fig. 4.50). The input must be padded with p >= 1; the output receives
 *  PAD_REPLICATE padding and may be the same as the input.
 */
int boundary_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int connectivity)
{
	unsigned char B[9], * E = NULL;
	int err;

	/* check parameters: */
	if (A == NULL || A_out == NULL || height < 0 || width < 0 || p < 1) return DFX_INVARG;
	if (connectivity != 4 && connectivity != 8) return DFX_INVARG;

	/* erode and subtract: */
	structuring_element(B, 1, connectivity == 8 ? SE_CROSS : SE_BOX);
	if ((err = alloc_binary(&E, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = erode_binary(A, E, width, height, p, B, 1)) == DFX_SUCCESS)
		err = subtract_binary(A, E, A_out, width, height, p);
	free_binary(E);
	if (err == DFX_SUCCESS) err = pad_binary(A_out, width, height, p, PAD_REPLICATE);
	return err;
}
