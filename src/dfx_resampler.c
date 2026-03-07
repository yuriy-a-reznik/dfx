/*!
 *  \file   dfx_resampler.c
 *  \brief  DFX library: image resampling functions
 * 
 *  This module implements separable, polyphase-filter-based resampler. 
 *  The implementation is very basic and intends only to explain the principles. 
 * 
 *  A good read on polyphase filter design is:
 *   R. E. Crochiere, and L. R. Rabiner, "Interpolation and decimation of digital 
 *   signals — A tutorial review," Proc. IEEE 69, no. 3, 2005, pp. 300-331.
 *
 *  I mostly follow this paper, adding a specific design of the low-pass filter, 
 *  and some other small details. Everything is done in linear space, floats, 
 *  and padding reserved in input and output images.
 * 
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "dfx.h"

 /****************
  *
  *  Functions supporting polyphase filter construction:
  *
  *   gen_coeffs()             - generate lp filter coefficients
  *   free_coeffs()            - free memory occupied by coefficients
  *   gcd()                    - compute greatest common divisor
  *   apply_gain()             - apply gain to filter coefficients
  *   gen_polyphase_coeffs()   - generate polyphase filter coefficients
  */

/*!
 *  \brief Generate an array of LP filter coefficients.
 *
 *  This function uses the classic 2n+1=point sinc filter smoothed by Lanczos window.
 *  This array of coefficients is allocated and returned by a pointer *p_w.
 *  The calling application is responsible for freeing it when it is no longer needed.
 * 
 *  \param[out]  p_w  - pointer to an array of filter coefficients to generate
 *  \param[in]   n    - filter width (max lag in each direction)
 *  \param[in]   fc   - cutoff frequency
 *
 *  \returns     DFX_X error codes
 */
static int gen_coeffs(float **p_w, int n, float fc)
{
    double s, l;
    float* w = NULL;
    int x;

    /* check args: */
    if (p_w == NULL || n < 1 ||  fc <= 0 || fc > 0.5)
        return DFX_INVARG;

    /* allocate an array for LP filter coefficients: */
    *p_w = NULL;
    if ((w = (float*)malloc((2 * n + 2) * sizeof(float))) == NULL) 
        return DFX_NOMEM;

    /* generate filter coefficients: */
    for (x = -n; x <= n; x++) {
        if (x != 0) {
            s = sin(2 * M_PI * fc * x) / (M_PI * x);                   /* sinc filter */
            l = sin(M_PI * x / (double)n) / (M_PI * x / (double) n);   /* Lanczos window */
            w[n + x] = (float)(s * l);
        } else {
            w[n + x] = 2 * fc;                                         /* limit at x --> 0 */
        }
    }

    /* success: */
    *p_w = w;
    return DFX_SUCCESS;
}

/*!
 *  Free memory allocated for the coefficients.
 */
static int free_coeffs(float* w)
{
    if (w != NULL) free(w);
    return DFX_SUCCESS;
}

/*!
  *  \brief Find greatest common divisor.
  *
  *   Euclid's algorithm:
  *            gcd(0,b) = 0
  *            gcd(a,b) = gcd(b%a,a);
  */
static int gcd(int a, int b)
{
    while (a != 0) { int c = a; a = b % a;  b = c; }
    return b;
}

/*!
  *  \brief Apply gain to filter coefficients.
  *
  *  \param[in,out]  p    - filter coefficients
  *  \param[in]      n    - filter width (max lag in each direction)
  *  \param[in]      gain - gain to apply
  *
  *  \returns        DFX_X error codes
  */
static int apply_gain(float* w, int n, float gain)
{
    double s; int x;

    /* check args: */
    if (w == NULL || n < 1) return DFX_INVARG;

    /* compute current gain: */
    for (x = -n, s = 0.; x <= n; x++) s += w[n + x];
    /* force new gain: */
    for (x = -n; x <= n; x++) w[n + x] *= (float)(gain / s);

    /* success */
    return DFX_SUCCESS;
}

/*!
 *  \brief Compute a set of polyphase filter coefficients.
 *
 *  This function generates an array of filter coefficients, defining filters for each phase.
 *  This array of coefficients is allocated and returned by a pointer *p_W.
 *  The calling application is responsible for freeing it when it is no longer needed.
 *
 *  \param[in,out]  p_W       - pointer to an array of filter coefficients, defining filters for each phase 
 *  \param[in]      dim_in    - input dimension of an image (width or height) 
 *  \param[in]      dim_out   - output dimension of an image (width or height)
 *  \param[in]      n         - filter width (max lag in each direction)
 *
 *  \returns        DFX_X error codes
 */
static int gen_polyphase_coeffs(float **p_W, int dim_in, int dim_out, int n)
{
    float *w = NULL, *W = NULL, fc, gain;
    int i, j, k, phase, M, N;

    /* check args: */
    if (p_W == NULL || n < 1 || dim_in <= 1 || dim_out <= 1)
        return DFX_INVARG;

    /* compute min fraction M/N: */
    k = gcd(dim_out, dim_in);
    M = dim_out / k;
    N = dim_in / k;

    /* set cutoff and gain: */
    fc = (M < N) ? 0.5f / (float)N : 0.5f / (float)M; 
    gain = (float)M;

    /* generate filter coefficients: */
    if (gen_coeffs(&w, M * n, fc) != DFX_SUCCESS)
        return DFX_NOMEM;

    /* allocate an array for per-phase sorted LP filter coefficients: */
    *p_W = NULL;
    if ((W = (float*)malloc(M * (2 * n + 1) * sizeof(float))) == NULL) {
        free_coeffs(w);
        return DFX_NOMEM;
    }
    
    /* subsample for each phase: */
    for (phase = 0; phase < M; phase++)
    {
        /* set start index, first coeff: */
        if (!phase) { i = 0; j = 0; }
        else { W[phase * (2*n + 1) + 0] = 0; i = 1; j = M - phase; }

        /* move coeffs: */
        for (; j < M * 2 * n + 1; j += M, i++)
            W[phase*(2*n+1) + i] = w[j];

        /* apply gain: */
        apply_gain(&W[phase * (2*n + 1)], n, 1);
    }

    /* success: */
    *p_W = W;
    free_coeffs(w);
    return DFX_SUCCESS;
}

/*!
 *  \brief Apply polyphase filter to a row in an image.
 *
 *  \param[in]  x          - pointer to a row in source image 
 *  \param[out] y          - pointer to a row in output image
 *  \param[in]  width_in   - source image width
 *  \param[in]  width_out  - output image width
 *  \param[in]  n          - filter length
 *  \param[in]  W          - polyphase filter coefficients
 *
 *  \returns        DFX_X error codes
 */
static int polyphase_filter_row(float* x, float* y, int width_in, int width_out, int n, float *W)
{
    int M, N, j, j0, k, l, phase;
    float s;

    /* check parameters: */
    if (x == NULL || y == NULL || W == NULL) return DFX_INVARG;
    if (width_in < 1 || width_out < 1) return DFX_INVARG;

    /* compute min fraction M/N: */
    k = gcd(width_out, width_in);
    M = width_out / k;
    N = width_in / k;

    /* filter row: */
    for (j = 0, l = 0; j < width_out; j++, l += N) {
        /* map to input pixel position: */
        j0 = l / M;                     /* j0 = floor(j*N/M) */
        phase = l % M;                  /* phase = j*N - j0*M */
        /* apply filter: */
        for (k = -n, s = 0.; k <= n; k++) 
            s += x[j0 + k] * W[phase * (2 * n + 1) + n + k];
        /* save the result: */
        y[j] = s;
    }

    /* success: */
    return DFX_SUCCESS;
}

/*!
 *  \brief Apply vertical polyphase filter.
 *
 *  \param[in]  x          - pointer to a column in source image
 *  \param[out] y          - pointer to a column in output image
 *  \param[in]  width      - image width (the same for input & output, including padding)
 *  \param[in]  height_in  - source image width 
 *  \param[in]  height_out - output image width 
 *  \param[in]  n          - filter length
 *  \param[in]  W          - polyphase filter coefficients
 *
 *  \returns        DFX_X error codes
 */
static int polyphase_filter_col(float* x, float* y, int width, int height_in, int height_out, int n, float* W)
{
    int M, N, i, i0, k, l, phase;
    float s;

    /* check parameters: */
    if (x == NULL || y == NULL || W == NULL) return DFX_INVARG;
    if (width < 2*n+1 || height_in < 1 || height_out < 1) return DFX_INVARG;

    /* compute min fraction M/N: */
    k = gcd(height_out, height_in);
    M = height_out / k;
    N = height_in / k;

    /* apply filter across columns: */
    for (i = 0, l = 0; i < height_out; i++, l += N) {
        /* map to input pixel position: */
        i0 = l / M;                   /* i0 = floor(i*N/M) */
        phase = l % M;                /* phase = i*N - i0*M */
        /* apply filter: */
        for (k = -n, s = 0.; k <= n; k++)
            s += x[(i0 + k) * width] * W[phase * (2 * n + 1) + n + k];
        y[i * width] = s;
    }

    /* success: */
    return DFX_SUCCESS;
}

/*!
 *  \brief Image plane resampler.
 * 
 *  This function implements separable row/column resampler utilizing polyphase filters. 
 *  Both input and output images are assumed to be padded with boundary (p) at least as large as 
 *  filter length (n) employed. 
 *  
 *  \param[in]  X           - input image plane
 *  \param[out] Y           - output resampled image
 *  \param[in]  width_in    - input image width
 *  \param[in]  height_in   - input image height
 *  \param[in]  width_out   - input image width
 *  \param[in]  height_out  - input image height
 *  \param[in]  p           - padding parameter
 *  \param[in]  n           - filter width (max lag in each direction)
 *
 *  \returns    DFX_X error codes
 */
int resample_plane(float* X, float* Y, int width_in, int height_in, int width_out, int height_out, int p, int n)
{
    /* local buffers: */
    float *Z = NULL;		                /* temporary plane holding the results of row-filtering */
    float *W_row = NULL, *W_col = NULL;	    /* filter coeffs for scaling horisontally and vertically */
    int x, y;

    /* check parameters */
    if (X == NULL || Y == NULL) return DFX_INVARG;
    if (height_in < 0 || width_in < 0 || p < 1 || n < 1 || n > p) return DFX_INVARG;    /* n must be <= p! */

    /* allocate temp image of size [out_width x in_height] + padding: */
    if (alloc_plane(&Z, width_out, height_in, p) != DFX_SUCCESS) return DFX_NOMEM;

    /* generate polyphase filter coefficients: */
    if (gen_polyphase_coeffs(&W_row, width_in, width_out, n) != DFX_SUCCESS) { free_plane(Z); return DFX_NOMEM; }
    if (gen_polyphase_coeffs(&W_col, height_in, height_out, n) != DFX_SUCCESS) { free_plane(Z); free_coeffs(W_row); return DFX_NOMEM; }

    /* resample rows, including extra +-n rows in the padded region to support subsequent vertical processing: */
    for (y = -n; y < height_in + n; y++) {
        /* call 1D filter (convolution) function on each row: */
        polyphase_filter_row(
            &X[(y + p) * (width_in + 2 * p) + p],
            &Z[(y + p) * (width_out + 2 * p) + p],
            width_in, width_out, n, W_row);
    }

    /* filter columns: */
    for (x = 0; x < width_out; x++) {
        /* call 1D filter (convolution) function on each column: */
        polyphase_filter_col(
            &Z[p * (width_out + 2 * p) + x + p],
            &Y[p * (width_out + 2 * p) + x + p],
            width_out + 2 * p, height_in, height_out, n, W_col);
    }

    /* free intermediate buffers & exit: */
    free_plane(Z);
    free_coeffs(W_row);
    free_coeffs(W_col);

    return DFX_SUCCESS;
}

/*!
 *  \brief Resample linear RGB image.
 */
int resample_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, 
    int width_in, int height_in, int width_out, int height_out, int p, int n)
{
    /* check parameters */
    if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
    if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
    if (width_in < 0 || height_in < 0 || width_out < 0 || height_out < 0) return DFX_INVARG;
    if (p < 1 || n < 1 || n > p) return DFX_INVARG;

    /* resample planes: */
    resample_plane(R_in, R_out, width_in, height_in, width_out, height_out, p, n);
    resample_plane(G_in, G_out, width_in, height_in, width_out, height_out, p, n);
    resample_plane(B_in, B_out, width_in, height_in, width_out, height_out, p, n);

    /* success: */
    return DFX_SUCCESS;
}

/* dfx_resampler.c -- end of file */