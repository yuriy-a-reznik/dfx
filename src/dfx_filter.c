/*!
 *  \file   dfx_filter.c
 *  \brief  DFX library: linear filtering operations
 *
 *  This module shows how to implement separable linear filters with symmetric kernels.  
 *  The implementation is very simple, but this is not a production code.  
 *  One would want to remove mallocs and SIMD optimize it to turn it into 
 *  something workable in production context. 

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
 *  \param[in,out]  p_w  - pointer to a pointer to a buffer for filter coefficients 
 *  \param[in]      n    - filter width (max lag in each direction)
 *  \param[in]      t    - filter kernel type (see FILT_X constants)
 *  \param[in]      fc   - cutoff frequency
 *
 *  \returns        DFX_X error codes
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
        double sigma = sqrt(M_LN2) / (2 * M_PI * fc);                       /* sigma for -3db @ cutoff */
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
 * Internal functions supporting the implementaion of a 2D filter:
 *
 *  filter_row()   - 1D filter along a row in an image
 *  filter_col()   - 1D filter along a column in an image
 */

/*!
 *  \brief Apply 1D filter along a row in an image.
 *
 *  \param[in]  x      - row in an input image
 *  \param[out] y      - row in a filtered image
 *  \param[in]  width  - image width
 *  \param[in]  n      - filter width
 *  \param[in]  w      - filter coefficients
 *
 *  \returns    DFX_X error codes
 */
static int filter_row(float* x, float* y, int width, int n, float* w)
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
 *  \param[in]  x      - row in an input image
 *  \param[out] y      - row in a filtered image
 *  \param[in]  width  - image width (including padding)
 *  \param[in]  height - image height
 *  \param[in]  n      - filter width
 *  \param[in]  w      - filter coefficients
 *
 *  \returns    DFX_X error codes
 */
static int filter_col(float* x, float* y, int width, int height, int n, float* w)
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
 *  This is a very basic, row-column separable implenmentation of a 2D filter. 
 *  The input image is anticipated to be padded, and the padding parameter p 
 *  must be at least as large as filter width n. 
 *
 *  \param[in]  X       - input image plane 
 *  \param[out] Y       - filtered image plane
 *  \param[in]  height  - image height
 *  \param[in]  width   - image width
 *  \param[in]  p       - padding parameter
 *  \param[in]  n       - filter width (max lag in each direction)
 *  \param[in]  t       - filter kernel type (see FILT_X constants)
 *  \param[in]  fc      - cutoff frequency
 *
 *  \returns    DFX_X error codes
 */
int filter_plane(float* X, float* Y, int width, int height, int p, int n, int t, float fc)
{
    /* local buffers & vars: */
    float* Z = NULL;	        /* temporary plane holding the results of row-filtering */
    float* w = NULL;            /* filter coefficients  */
    int x, y;

    /* check parameters */
    if (X == NULL || Y == NULL) return DFX_INVARG;
    if (height < 0 || width < 0 || p < 1 || n < 1 || n > p) return DFX_INVARG;  /* n must be <= p - padding! */
    if (t < 0 || t >= N_FILTERS || fc < 0 || fc > 0.5) return DFX_INVARG;

    /* allocate temporary image planes & arrays of filter coeffs: */
    if (alloc_plane(&Z, width, height, p) != DFX_SUCCESS) return DFX_NOMEM;
    if (gen_coeffs(&w, n, t, fc) != DFX_SUCCESS) { free_plane(Z); return DFX_NOMEM;}

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

    /* free intermediate buffers & exit: */
    free_plane(Z);
    free_coeffs(w);

    return DFX_SUCCESS;
}

/*!
 *  \brief Filter linear RGB image.
 */
int filter_image(float* R_in, float* G_in, float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, int n, int t, float fc)
{
    /* check parameters */
    if (R_in == NULL || G_in == NULL || B_in == NULL) return DFX_INVARG;
    if (R_out == NULL || G_out == NULL || B_out == NULL) return DFX_INVARG;
    if (height < 0 || width < 0 || p < 1 || n < 1 || n > p) return DFX_INVARG;         /* n must be <= p! */
    if (t < 0 || t >= N_FILTERS || fc < 0 || fc > 0.5) return DFX_INVARG;

    /* filter planes: */
    filter_plane(R_in, R_out, width, height, p, n, t, fc);
    filter_plane(G_in, G_out, width, height, p, n, t, fc);
    filter_plane(B_in, B_out, width, height, p, n, t, fc);

    /* success: */
    return DFX_SUCCESS;
}

/* dfx_filter.c -- end of file */