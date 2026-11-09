/*!
 *  \file   dfx_filter.c
 *  \brief  Image filtering operations.
 *
 *  This module shows how to implement separable linear filters with symmetric kernels.
 *  Filtering of images is performed in linear RGB space, with identical kernels applied 
 *  separately for each channel.
 *
 *  Three low-pass kernels of 2n+1 taps are provided, all normalized to unit DC gain:
 *    FILT_SINC    - the ideal low-pass impulse response sin(2 pi fc x) / (pi x), truncated to |x| <= n
 *                   (a rectangular window, which produces Gibbs ringing near edges);
 *    FILT_LANCZOS - the same sinc multiplied by the Lanczos window sinc(x/n) = sin(pi x/n) / (pi x/n);
 *    FILT_GAUSS   - the Gaussian exp(-x^2 / (2 sigma^2)) with sigma chosen so that its frequency
 *                   response is -3 dB (1/sqrt(2)) at fc: since the Fourier transform of the Gaussian is
 *                   exp(-2 pi^2 sigma^2 f^2), this gives sigma = sqrt(ln 2) / (2 pi fc).
 *  Here fc is the cutoff frequency in cycles per pixel (0.5 = the Nyquist frequency).
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    November 8, 2026
 */

#include <math.h>
#include <stdlib.h>
#include "dfx.h"

/*****************
 *
 * Operations with filter coefficients:
 * 
 *  gen_coeffs()   - generates an array of filter coefficients
 *  free_coeffs()  - frees an array of filter coefficients 
 */

/*!
 *  \brief Generate an array of filter coefficients.
 *
 *  This function generates an array of filter coefficients for a specified low-pass 
 *  filter type and its cutoff frequency and length.
 * 
 *  \param[out]     p_w  - pointer to a variable to receive the pointer to the allocated array of 2n+1 coefficients 
 *  \param[in]      n    - filter half-length (max lag in each direction)
 *  \param[in]      t    - filter kernel type (see FILT_X constants)
 *  \param[in]      fc   - cutoff frequency in cycles per pixel, 0 < fc <= 0.5
 *
 *  \returns        DFX_X error code
 */
static int gen_coeffs(float** p_w, int n, int t, float fc)
{
    float *w = NULL;
    double s, l;
    int x;

    /* check args: */
    if (p_w == NULL || n < 1 || t < 0 || t >= N_FILTERS || fc <= 0 || fc > 0.5)
        return DFX_INVARG;

    /* allocate an array for filter coefficients: */
    *p_w = NULL;
    if ((w = (float*)malloc((2*n + 1) * sizeof(float))) == NULL) 
        return DFX_NOMEM;

    /* check filter type: */
    if (t == FILT_GAUSS)
    {
        /* Gaussian filter: */
        double sigma = sqrt(M_LN2) / (2 * M_PI * fc);                       /* sigma for -3 dB @ cutoff */
        for (x = -n; x <= n; x++)
            w[n + x] = (float)(1. / (sqrt(2 * M_PI) * sigma) * exp(-x * x / (2 * sigma * sigma))); 
    }
    else if (t == FILT_SINC)
    {
        /* sinc filter: */
        for (x = -n; x <= n; x++) {
            if (x != 0)  
                w[n + x] = (float)(sin(2 * M_PI * fc * x) / (M_PI * x));
            else         
                w[n + x] = 2 * fc;                                          /* limit at x --> 0 */
        }
    }
    else /* (t == FILT_LANCZOS) */
    {
        /* Lanczos filter: */
        for (x = -n; x <= n; x++) {
            if (x != 0) {
                s = sin(2 * M_PI * fc * x) / (M_PI * x);                    /* sinc filter */
                l = sin(M_PI * x / (double) n) / (M_PI * x / (double) n);   /* Lanczos window */
                w[n + x] = (float)(s*l);
            } else 
                w[n + x] = 2 * fc;                                          /* limit at x --> 0 */
        }
    }

    /* force unit gain: */
    for (x = -n, s = 0.; x <= n; x++) s += w[n + x];
    for (x = -n; x <= n; x++) w[n + x] = (float)(w[n + x] / s);

    /* success: */
    *p_w = w;
    return DFX_SUCCESS;
}

/*!
 *  Frees transform coefficients
 */
static int free_coeffs(float* w)
{
    if (w != NULL) free(w);
    return DFX_SUCCESS;
}

/*****************
 *
 * Internal functions supporting the implementation of a 2D filter:
 *
 *  filter_row()   - 1D filter along a row in an image
 *  filter_col()   - 1D filter along a column in an image
 */

/*!
 *  \brief Apply 1D filter along a row in an image.
 *
 *  \param[in]  x      - row in an input image (must have n valid samples on each side)
 *  \param[out] y      - row in a filtered image
 *  \param[in]  width  - image width
 *  \param[in]  n      - filter half-length
 *  \param[in]  w      - 2n+1 filter coefficients
 *
 *  \returns    DFX_X error code
 */
static int filter_row(const float* x, float* y, int width, int n, const float* w)
{
    int i, j; double s;

    /* check parameters */
    if (x == NULL || y == NULL || width < 0 || n<1 || w == NULL) 
        return DFX_INVARG;

    /* apply filter: */
    for (i = 0; i < width; i++) {
        for (j = -n, s = 0.; j <= n; j++)
            s += x[i + j] * w[j + n];
        y[i] = (float)s;
    }

    return DFX_SUCCESS;
}

/*!
 *  \brief Apply 1D filter along a column in an image
 *
 *  \param[in]  x      - column in an input image (must have n valid samples above and below)
 *  \param[out] y      - column in a filtered image
 *  \param[in]  width  - image width (including padding) = distance between consecutive samples of a column
 *  \param[in]  height - image height
 *  \param[in]  n      - filter half-length
 *  \param[in]  w      - 2n+1 filter coefficients
 *
 *  \returns    DFX_X error code
 */
static int filter_col(const float* x, float* y, int width, int height, int n, const float* w)
{
    int i, j; double s;

    /* check parameters */
    if (x == NULL || y == NULL || height < 0 || width < 0 || n < 1 || w == NULL) 
        return DFX_INVARG;

    /* apply filter: */
    for (i = 0; i < height; i++) {
        for (j = -n, s = 0.; j <= n; j++)
            s += x[(i + j) * width] * w[j + n];
        y[i * width] = (float)s;
    }

    return DFX_SUCCESS;
}

/*!
 *  \brief Apply 2-dimensional filter over an image plane
 * 
 *  This is a very basic, row-column separable implementation of a 2D filter. 
 *  The input image is anticipated to be padded, and the padding parameter p 
 *  must be at least as large as the filter half-length n. The output plane 
 *  receives PAD_REPLICATE padding, so it can be fed into further operations.
 *  The output plane may be the same as the input plane.
 *
 *  \param[in]  X       - input image plane 
 *  \param[out] Y       - filtered image plane
 *  \param[in]  width   - image width
 *  \param[in]  height  - image height
 *  \param[in]  p       - padding parameter
 *  \param[in]  n       - filter half-length (max lag in each direction), 1 <= n <= p
 *  \param[in]  t       - filter kernel type (see FILT_X constants)
 *  \param[in]  fc      - cutoff frequency in cycles per pixel, 0 < fc <= 0.5
 *
 *  \returns    DFX_X error code
 */
int filter_plane(const float* X, float* Y, int width, int height, int p, int n, int t, float fc)
{
    /* local buffers & vars: */
    float* Z = NULL;	        /* temporary plane holding the results of row-filtering */
    float* w = NULL;            /* filter coefficients  */
    int x, y, err;

    /* check parameters */
    if (X == NULL || Y == NULL) return DFX_INVARG;
    if (height < 0 || width < 0 || p < 1 || n < 1 || n > p) return DFX_INVARG;  /* n must be <= p - padding! */
    if (t < 0 || t >= N_FILTERS || fc <= 0 || fc > 0.5) return DFX_INVARG;

    /* allocate temporary image planes & arrays of filter coeffs: */
    if ((err = alloc_plane(&Z, width, height, p)) != DFX_SUCCESS) return err;
    if ((err = gen_coeffs(&w, n, t, fc)) != DFX_SUCCESS) { free_plane(Z); return err; }

    /* filter rows, include extra +-n rows in the padded region: */
    for (y = -n; y < height+n; y++) {
        /* call 1D filter (convolution) function on each row: */
        filter_row(
            &X[(y + p) * (width + 2 * p) + p],
            &Z[(y + p) * (width + 2 * p) + p],
            width, n, w);
    }

    /* filter columns: */
    for (x = 0; x < width; x++) {
        /* call 1D filter (convolution) function on each column: */
        filter_col(
            &Z[p * (width + 2 * p) + x + p],
            &Y[p * (width + 2 * p) + x + p],
            width + 2 * p, height, n, w);
    }

    /* free intermediate buffers: */
    free_plane(Z);
    free_coeffs(w);

    /* fill the padding of the output plane & exit: */
    return pad_plane(Y, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Filter linear RGB image (see filter_plane()).
 */
int filter_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, int n, int t, float fc)
{
    int err;

    /* check parameters */
    if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
    if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
    if (height < 0 || width < 0 || p < 1 || n < 1 || n > p) return DFX_INVARG;         /* n must be <= p! */
    if (t < 0 || t >= N_FILTERS || fc <= 0 || fc > 0.5) return DFX_INVARG;

    /* filter planes: */
    if ((err = filter_plane(R_in, R_out, width, height, p, n, t, fc)) != DFX_SUCCESS) return err;
    if ((err = filter_plane(G_in, G_out, width, height, p, n, t, fc)) != DFX_SUCCESS) return err;
    return filter_plane(B_in, B_out, width, height, p, n, t, fc);
}

/* dfx_filter.c -- end of file */
