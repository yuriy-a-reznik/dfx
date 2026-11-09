/*!
 *  \file   dfx.h
 *  \brief  DFX: Digital image processing Functions and eXamples.
 *
 *  This library provides a collection of basic image-processing primitives
 *  implemented entirely in C, without relying on external libraries.
 *  It is written as a companion to an image-processing course: every
 *  operation follows the textbook definition rather than a production
 *  shortcut, and all color, filtering, and resampling operations are
 *  performed on linear-light (not gamma-companded) pixel values.
 *
 *  Image representations
 *  ---------------------
 *
 *  *  sRGB images. read_bitmap() and write_bitmap() load and store 24-bit
 *     packed sRGB images. In memory, such an image is an array of unsigned
 *     chars holding B,G,R byte triplets, row by row, with each row padded to
 *     a multiple of 4 bytes and rows stored bottom-up, exactly as in a .bmp
 *     file. The parameters "width" and "height" specify the image dimensions
 *     in pixels. Use alloc_srgb_image() to allocate such an image.
 *
 *  *  Linear images. srgb_to_linear() and linear_to_srgb() convert between
 *     sRGB and linear RGB. In linear form, an image is stored as three float
 *     planes R, G, and B (the linear red, green, and blue components, each
 *     in the range [0,1] for in-gamut colors), with rows stored top-down.
 *     Linear representation is fundamental: all color mixing, filtering, and
 *     resampling operations in this library are implemented in linear space.
 *
 *  *  Padding. Every linear plane is allocated with a border of p pixels on
 *     each side, so a plane for a width x height image occupies
 *     (width + 2p) x (height + 2p) floats, and the pixel (x,y) of the image
 *     is stored at index (p + y) * (width + 2p) + (p + x). The parameter "p"
 *     is passed to every function working with planes. Padding supplies the
 *     out-of-image samples needed by filtering and resampling; pad_plane()
 *     and pad_image() fill the border by one of the PAD_* rules. Functions
 *     that produce new planes (srgb_to_linear(), filter_plane(),
 *     resample_plane(), the color-space conversions and adjustments) leave
 *     their results replicate-padded, so results can be fed into further
 *     operations directly. The DFT functions write only the interior of
 *     their output planes.
 *
 *  *  Single planes. Functions whose names end in _plane operate on one
 *     float plane; functions ending in _image apply the same operation to
 *     the three planes of an RGB image. Planes may also hold other
 *     quantities (luminance, XYZ or L*u*v* components, DFT spectra,
 *     distances).
 *
 *  *  Binary images. A binary image is a set of pixels, stored as a plane
 *     of unsigned chars with the values 0 (OFF, background) and 1 (ON,
 *     foreground) and exactly the same padded geometry as a float plane.
 *     Functions whose names end in _binary operate on binary images. The
 *     operations of mathematical morphology are set operations and are kept
 *     strictly apart from the arithmetic of linear planes: a plane becomes
 *     a binary image only through thresholding (threshold_plane() and its
 *     relatives), and a binary image becomes a plane again only through
 *     binary_to_plane(). Connected components are numbered in label planes,
 *     arrays of ints of the same geometry holding 0 for the background.
 *
 *  Modules
 *  -------
 *
 *  *  dfx_bmp.c       - bitmap file I/O
 *  *  dfx_srgb.c      - sRGB <-> linear conversions, 8-bit quantization and dithering
 *  *  dfx_image.c     - allocation, padding, copy/scale/add/blend/clip of planes
 *  *  dfx_filter.c    - separable low-pass filtering (sinc, Lanczos, Gaussian kernels)
 *  *  dfx_resampler.c - rational-factor polyphase resampling
 *  *  dfx_dft.c       - 2D DFT by direct matrix multiplication, magnitude and phase
 *  *  dfx_fft.c       - fast 1D FFT algorithms (radix-2 and prime-factor)
 *  *  dfx_dct.c       - fast 1D DCT-II and DCT-IV algorithms
 *  *  dfx_colors.c    - RGB <-> XYZ <-> L*u*v* / L*a*b* conversions, color adjustments
 *  *  dfx_snr.c       - MSE and SNR metrics
 *  *  dfx_binary.c    - binary images: thresholding, structuring elements, dilation, erosion, opening,
 *                       closing, majority filter, hit-or-miss, thinning, thickening, boundaries
 *  *  dfx_regions.c   - regions: floodfill, connected components, region properties, boundary tracing
 *                       and chain codes, chamfer distance transform, region growing
 *
 *  Error handling
 *  --------------
 *
 *  Every function returns DFX_SUCCESS (0) or one of the DFX_* error codes
 *  below. Unless stated otherwise, output arrays are left unspecified when
 *  an error is returned.
 *
 *  Examples demonstrating how to use these operations can be found in
 *  the dfx/examples directory.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.2
 *  \date    November 8, 2026
 */

#ifndef DFX_H
#define DFX_H  1

#include <math.h>

/* M_PI and M_LN2 are not part of standard C. They are provided by most
 * <math.h> implementations as extensions (MSVC requires _USE_MATH_DEFINES,
 * glibc hides them in strict -std=c11 mode), so we define fallbacks here. */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2  0.69314718055994530942
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! Filter kernels (parameter t of filter_plane() and filter_image()): */
enum {
	FILT_SINC = 0,		/* truncated (rectangular-windowed) sinc filter */
	FILT_LANCZOS = 1,	/* Lanczos-windowed sinc filter */
	FILT_GAUSS = 2,		/* Gaussian filter */
	/* ... */
	N_FILTERS			/* the number of supported filters */
};

/*! Image padding rules (parameter padding_type of pad_plane() and pad_image()): */
enum {
	PAD_ZERO = 0,		/* padding by zeros */
	PAD_REPLICATE = 1,	/* padding by replicating boundary pixels */
	PAD_REFLECT = 2,	/* padding by reflecting the image about its boundary (the boundary pixel is repeated) */
	N_PADDINGS			/* the number of supported padding rules */
};

/*! RGB primary colors corresponding to various color systems used in practice: */
enum {
	COLORS_BT709 = 0,   /* ITU-R Rec. BT.709 primaries, D65 white */
	COLORS_SRGB = 0,    /* sRGB primaries = the same as BT.709 */
	COLORS_NTSC = 1,    /* NTSC 1953 primaries, illuminant C white - legacy */
	COLORS_SMPTEC = 2,  /* SMPTE RP 145 (aka SMPTE-C), used by NTSC since 1987 */
	COLORS_EBU3213 = 3, /* EBU Tech 3213, used by PAL and SECAM systems since 1975 */
	COLORS_PAL = 3,     /* PAL/SECAM primaries = the same as EBU Tech 3213 */
	COLORS_BT2020 = 4,  /* ITU-R Rec. BT.2020 primaries */
	COLORS_P3D65 = 5,	/* DCI-P3 primaries + D65 white */
	COLORS_P3 = 6,		/* DCI-P3 primaries + DCI white */
	N_COLORS			/* the number of supported color systems */
};

/*! Error codes returned by DFX library functions: */
enum {
	DFX_SUCCESS = 0,	/* success */
	DFX_INVARG = 1,		/* invalid arguments */
	DFX_CANTOPEN = 2,	/* file cannot be opened/created */
	DFX_IOERR = 3,		/* file reading/writing error */
	DFX_NOTSUP = 4,		/* not supported input format, transform size, or operation */
	DFX_NOMEM = 5,		/* out of memory */
	DFX_INVRANGE = 6	/* pixel values are outside of the valid range */
};

/*! Shapes of structuring elements (parameter shape of structuring_element()): */
enum {
	SE_CROSS = 0,		/* |i| + |j| <= n: the city-block ball; for n = 1 the 4-neighborhood (B4) */
	SE_BOX = 1,			/* max(|i|, |j|) <= n: the chessboard ball; for n = 1 the 8-neighborhood (B8) */
	SE_DISK = 2,		/* i^2 + j^2 <= n^2: the discrete Euclidean ball */
	N_SES				/* the number of supported shapes */
};

/*! Values of a ternary structuring element (parameter B of hitmiss_binary()): */
enum {
	SE_OFF = 0,			/* the pixel must be OFF */
	SE_ON = 1,			/* the pixel must be ON */
	SE_DONTCARE = 2		/* the pixel is not examined */
};

/*! Boundary length estimators (parameter method of chain_code_length()): */
enum {
	CHAIN_FREEMAN = 0,		/* N_o + sqrt(2) N_d */
	CHAIN_PYTHAGOREAN = 1,	/* sqrt(N_d^2 + (N_o + N_d)^2) */
	CHAIN_KIMURA = 2,		/* sqrt(N_d^2 + (N_d + N_o/2)^2) + N_o/2 */
	N_CHAIN_METHODS			/* the number of supported estimators */
};

/*! Properties of a region of a binary image (filled by region_properties()): */
typedef struct {
	int    area;				/* number of pixels (the zeroth-order moment m00) */
	int    x_min, x_max;		/* bounding box: horizontal extent, inclusive */
	int    y_min, y_max;		/* bounding box: vertical extent, inclusive */
	double x_c, y_c;			/* centroid (m10 / m00, m01 / m00) */
	double mu20, mu11, mu02;	/* central second-order moments (sums over pixels, not divided by the area) */
	double theta;				/* orientation of the major axis of the equivalent ellipse: the angle from the x axis towards the y axis, in (-pi/2, pi/2] */
	double eccentricity;		/* eccentricity of the equivalent ellipse: 0 for a disk, approaching 1 for a line */
	double a, b;				/* semi-axes of the equivalent ellipse, a >= b */
	int    perimeter;			/* number of boundary pixels: region pixels with a 4-neighbor outside the region */
	double perimeter_length;	/* length of the outer boundary estimated by Kimura's formula from its chain code */
	double compactness;			/* 4 pi area / perimeter_length^2: about 1 for a disk (slightly above 1 for digitized disks), less for other shapes */
} dfx_region;

/********************
 *
 *  sRGB images and bitmap files (dfx_bmp.c)
 */

/*!
 *  \brief Allocate an sRGB image of the given dimensions.
 *  \param[out] p_srgb  - receives a pointer to the allocated image; free it with free_srgb_image()
 *  \param[in]  width   - image width in pixels
 *  \param[in]  height  - image height in pixels
 */
extern int alloc_srgb_image(unsigned char** p_srgb, int width, int height);

/*! \brief Deallocate an sRGB image allocated by alloc_srgb_image() or read_bitmap(). */
extern int free_srgb_image(unsigned char* srgb);

/*!
 *  \brief Read a 24-bit uncompressed .bmp file into a newly allocated sRGB image.
 *  \param[in]  filename  - name of the bitmap file to read
 *  \param[out] p_srgb    - receives a pointer to the allocated sRGB image; free it with free_srgb_image()
 *  \param[out] p_width   - receives the image width in pixels
 *  \param[out] p_height  - receives the image height in pixels
 *  \param[out] p_ppm_x   - receives the horizontal pixel density (pixels per meter; 0 if the file does not specify it)
 *  \param[out] p_ppm_y   - receives the vertical pixel density (pixels per meter; 0 if the file does not specify it)
 *  \return DFX_CANTOPEN if the file cannot be opened, DFX_NOTSUP if it is not a 24-bit uncompressed bitmap,
 *          DFX_IOERR if it is truncated or cannot be read.
 */
extern int read_bitmap(const char* filename, unsigned char** p_srgb, int* p_width, int* p_height, int* p_ppm_x, int* p_ppm_y);

/*!
 *  \brief Write an sRGB image into a 24-bit uncompressed .bmp file.
 *  \param[in] filename  - name of the bitmap file to create (overwritten if it exists)
 *  \param[in] srgb      - the sRGB image
 *  \param[in] width     - image width in pixels
 *  \param[in] height    - image height in pixels
 *  \param[in] ppm_x     - horizontal pixel density in pixels per meter (0 = unspecified)
 *  \param[in] ppm_y     - vertical pixel density in pixels per meter (0 = unspecified)
 */
extern int write_bitmap(const char* filename, const unsigned char* srgb, int width, int height, int ppm_x, int ppm_y);

/********************
 *
 *  Conversions between sRGB and linear RGB (dfx_srgb.c)
 */

/*!
 *  \brief Convert an 8-bit sRGB image to linear RGB planes (inverse sRGB transfer function).
 *
 *  The output planes must be allocated for the given width, height, and padding p.
 *  The border of each output plane is filled with PAD_REPLICATE padding.
 *
 *  \param[in]  sRGB     - input sRGB image
 *  \param[out] R, G, B  - output linear planes
 *  \param[in]  width    - image width in pixels
 *  \param[in]  height   - image height in pixels
 *  \param[in]  p        - padding parameter of the output planes
 */
extern int srgb_to_linear(const unsigned char* sRGB, float* R, float* G, float* B, int width, int height, int p);

/*!
 *  \brief Convert linear RGB planes to an 8-bit sRGB image (sRGB transfer function and rounding).
 *
 *  Values outside of [0,1] are clipped. Padding is dropped.
 *
 *  \param[in]  R, G, B  - input linear planes
 *  \param[out] sRGB     - output sRGB image, allocated by alloc_srgb_image() for the same width and height
 *  \param[in]  width    - image width in pixels
 *  \param[in]  height   - image height in pixels
 *  \param[in]  p        - padding parameter of the input planes
 */
extern int linear_to_srgb(const float* R, const float* G, const float* B, unsigned char* sRGB, int width, int height, int p);

/*!
 *  \brief Convert linear RGB planes to an 8-bit sRGB image using Floyd-Steinberg error diffusion.
 *
 *  Same as linear_to_srgb(), but the quantization error of each pixel (measured in linear
 *  space) is diffused to its right and lower neighbors, which suppresses banding in smooth gradients.
 *  The input planes are not modified.
 */
extern int linear_to_srgb_dithered(const float* R, const float* G, const float* B, unsigned char* sRGB, int width, int height, int p);

/********************
 *
 *  Linear image planes: allocation, initialization, padding, arithmetic (dfx_image.c)
 *
 *  All functions take the image dimensions (width, height) and the padding
 *  parameter p of the planes they operate on. Arithmetic operations
 *  (copy, scale, add, subtract, blend, clip) process the padded area as
 *  well, so the padding of their results stays valid. The output plane of
 *  an arithmetic operation may be the same as one of its inputs.
 */

/*!
 *  \brief Allocate a plane for a width x height image with padding p.
 *  \param[out] pX - receives a pointer to the allocated plane; free it with free_plane()
 */
extern int alloc_plane(float** pX, int width, int height, int p);
/*! \brief Deallocate a plane. */
extern int free_plane(float* X);
/*! \brief Allocate the three planes of an RGB image. */
extern int alloc_image(float** pR, float** pG, float** pB, int width, int height, int p);
/*! \brief Deallocate the three planes of an RGB image. */
extern int free_image(float* R, float* G, float* B);

/*! \brief Fill a plane (padding included) with zeros. */
extern int zero_plane(float* X, int width, int height, int p);
/*! \brief Fill an image (padding included) with zeros. */
extern int zero_image(float* R, float* G, float* B, int width, int height, int p);
/*! \brief Fill a plane (padding included) with ones. */
extern int unit_plane(float* X, int width, int height, int p);
/*! \brief Fill an image (padding included) with ones. */
extern int unit_image(float* R, float* G, float* B, int width, int height, int p);

/*!
 *  \brief Fill the padding border of a plane.
 *  \param[in,out] X            - the plane; its interior is left unchanged
 *  \param[in]     padding_type - one of PAD_ZERO, PAD_REPLICATE, PAD_REFLECT
 *  \return DFX_INVARG if padding_type is unknown or if PAD_REFLECT is requested with p larger than width or height.
 */
extern int pad_plane(float* X, int width, int height, int p, int padding_type);
/*! \brief Fill the padding borders of the three planes of an image (see pad_plane()). */
extern int pad_image(float* R, float* G, float* B, int width, int height, int p, int padding_type);

/*! \brief Copy a plane: X_out = X_in. */
extern int copy_plane(const float* X_in, float* X_out, int width, int height, int p);
/*! \brief Copy an image. */
extern int copy_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p);
/*! \brief Multiply a plane by a constant: X_out = scale * X_in. */
extern int scale_plane(const float* X_in, float* X_out, int width, int height, int p, float scale);
/*! \brief Multiply an image by a constant. */
extern int scale_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float scale);
/*! \brief Add two planes: X_out = X_1 + X_2. */
extern int add_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p);
/*! \brief Add two images. */
extern int add_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p);
/*! \brief Subtract two planes: X_out = X_1 - X_2. */
extern int subtract_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p);
/*! \brief Subtract two images. */
extern int subtract_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p);
/*! \brief Blend two planes: X_out = alpha * X_1 + (1 - alpha) * X_2, with alpha in [0,1]. */
extern int blend_planes(const float* X_1, const float* X_2, float* X_out, int width, int height, int p, float alpha);
/*! \brief Blend two images. */
extern int blend_images(const float* R_1, const float* G_1, const float* B_1, const float* R_2, const float* G_2, const float* B_2, float* R_out, float* G_out, float* B_out, int width, int height, int p, float alpha);
/*! \brief Clip the values of a plane to the range [lo, hi]: X_out = min(max(X_in, lo), hi). */
extern int clip_plane(const float* X_in, float* X_out, int width, int height, int p, float lo, float hi);
/*! \brief Clip the values of an image to the range [lo, hi]. */
extern int clip_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float lo, float hi);

/********************
 *
 *  Linear filtering (dfx_filter.c)
 */

/*!
 *  \brief Apply a separable symmetric low-pass filter to a plane.
 *
 *  The 2D filter is the outer product of a 1D kernel of 2n+1 taps, applied along
 *  rows and then along columns. The kernel is normalized to unit DC gain. The input
 *  plane must be padded (see pad_plane()) with p >= n. The output plane receives
 *  PAD_REPLICATE padding.
 *
 *  \param[in]  X      - input plane
 *  \param[out] X_out  - output plane (may be the same as X)
 *  \param[in]  width  - image width in pixels
 *  \param[in]  height - image height in pixels
 *  \param[in]  p      - padding parameter of both planes
 *  \param[in]  n      - filter half-length: the kernel covers lags -n..n (2n+1 taps); 1 <= n <= p
 *  \param[in]  t      - kernel type: FILT_SINC, FILT_LANCZOS, or FILT_GAUSS
 *  \param[in]  fc     - cutoff frequency in cycles per pixel, 0 < fc <= 0.5 (0.5 = Nyquist);
 *                       for FILT_GAUSS it is the -3 dB frequency
 */
extern int filter_plane(const float* X, float* X_out, int width, int height, int p, int n, int t, float fc);
/*! \brief Apply the same low-pass filter to the three planes of an image (see filter_plane()). */
extern int filter_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, int n, int t, float fc);

/********************
 *
 *  Resampling (dfx_resampler.c)
 */

/*!
 *  \brief Resample a plane to new dimensions using separable polyphase filters.
 *
 *  Each dimension is resampled by the rational factor dim_out/dim_in, using a
 *  Lanczos-windowed sinc low-pass prototype with cutoff at the smaller of the input and
 *  output Nyquist frequencies. Output sample j is placed at input coordinate j*dim_in/dim_out.
 *  Both planes must be allocated with the same padding p >= n; the input must be padded,
 *  and the output receives PAD_REPLICATE padding.
 *
 *  \param[in]  X           - input plane
 *  \param[out] X_out       - output plane
 *  \param[in]  width_in    - input image width in pixels (> 1)
 *  \param[in]  height_in   - input image height in pixels (> 1)
 *  \param[in]  width_out   - output image width in pixels (> 1)
 *  \param[in]  height_out  - output image height in pixels (> 1)
 *  \param[in]  p           - padding parameter of both planes
 *  \param[in]  n           - filter half-length at the output sampling rate (each output sample uses 2n+1 input samples); 1 <= n <= p
 */
extern int resample_plane(const float* X, float* X_out, int width_in, int height_in, int width_out, int height_out, int p, int n);
/*! \brief Resample the three planes of an image (see resample_plane()). */
extern int resample_image(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width_in, int height_in, int width_out, int height_out, int p, int n);

/********************
 *
 *  Discrete Fourier Transform (dfx_dft.c)
 */

/*!
 *  \brief Compute the 2D DFT of a real plane by direct (matrix-multiply) evaluation.
 *
 *  F(u,v) = sum_{x,y} f(x,y) exp(-2 pi i (u x / width + v y / height)), computed separably
 *  in O(width*height*(width+height)) operations. Only the interior of the plane is transformed.
 *
 *  \param[in]  f         - input plane (real)
 *  \param[out] reF, imF  - real and imaginary parts of the spectrum, stored as planes of the same geometry
 *  \param[in]  width     - image width in pixels
 *  \param[in]  height    - image height in pixels
 *  \param[in]  p         - padding parameter of all planes
 *  \param[in]  center_dc - if nonzero, f is multiplied by (-1)^(x+y) before the transform, which moves
 *                          the DC term to (width/2, height/2); exact for even width and height
 */
extern int dft_plane(const float* f, float* reF, float* imF, int width, int height, int p, int center_dc);
/*! \brief Compute the magnitude |F| = sqrt(reF^2 + imF^2) of a spectrum, element by element. */
extern int dft_magnitude(const float* reF, const float* imF, float* magF, int width, int height, int p);
/*! \brief Compute the phase atan2(imF, reF) in (-pi, pi] of a spectrum, element by element. */
extern int dft_phase(const float* reF, const float* imF, float* phaseF, int width, int height, int p);

/********************
 *
 *  Fast 1D transforms (dfx_fft.c, dfx_dct.c)
 */

/*!
 *  \brief In-place complex FFT of length n: X[k] = sum_j x[j] exp(-2 pi i j k / n), unnormalized.
 *
 *  Uses fft_dit() when n is a power of two and fft_pfa() otherwise.
 *
 *  \param[in,out] re, im - real and imaginary parts of the input, replaced by those of the spectrum
 *  \param[in]     n      - transform length
 *  \return DFX_NOTSUP if the length is supported by neither algorithm.
 */
extern int fft(float* re, float* im, int n);
/*! \brief In-place radix-2 decimation-in-time FFT; n must be a power of two. Same interface as fft(). */
extern int fft_dit(float* re, float* im, int n);
/*! \brief In-place radix-2 decimation-in-frequency FFT; n must be a power of two. Same interface as fft(). */
extern int fft_dif(float* re, float* im, int n);
/*!
 *  \brief In-place prime-factor FFT built from Winograd short-length modules. Same interface as fft().
 *
 *  n must be a product of mutually prime factors taken from {2,3,4,5,7,8,9,16}, that is,
 *  n = 2^a * 3^b * 5^c * 7^d with a <= 4, b <= 2, c <= 1, d <= 1 (for example 720, 1008, 5040).
 */
extern int fft_pfa(float* re, float* im, int n);

/*!
 *  \brief In-place DCT-II of length n: X[k] = sum_j x[j] cos(pi (2j+1) k / (2n)), unnormalized.
 *
 *  Computed by recursive splitting into half-length DCT-II and DCT-IV transforms down to
 *  short-length modules. Supported lengths: n = 2^m * n1 with n1 in {2, 3, 4, 5, 7, 9, 15}.
 *
 *  \param[in,out] x - input vector, replaced by its transform
 *  \param[out]    y - scratch buffer of at least n floats
 *  \param[in]     n - transform length
 *  \return DFX_NOTSUP if the length is not supported.
 */
extern int dct2_rec(float* x, float* y, int n);
/*!
 *  \brief In-place DCT-IV of length n: X[k] = sum_j x[j] cos(pi (2j+1)(2k+1) / (4n)), unnormalized.
 *  Same supported lengths and parameters as dct2_rec().
 */
extern int dct4_rec(float* x, float* y, int n);

/********************
 *
 *  Color spaces (dfx_colors.c)
 *
 *  The parameter "colors" selects the RGB primaries and white point (see the
 *  COLORS_* enum). The parameter "Y_white" is the luminance assigned to RGB
 *  white (1,1,1) in XYZ space; it also serves as the reference white
 *  luminance Y_n of the L*u*v* and L*a*b* spaces. Use Y_white = 1 for a
 *  relative (normalized) representation.
 */

/*!
 *  \brief Convert linear RGB planes to CIE 1931 XYZ planes.
 *  \return DFX_INVRANGE (and no output) if any RGB value is outside of [0,1]; see clip_image().
 */
extern int linear_rgb_to_xyz(const float* R, const float* G, const float* B, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white);
/*!
 *  \brief Convert CIE 1931 XYZ planes to linear RGB planes.
 *  \return DFX_INVRANGE if any resulting RGB value is outside of [0,1] (the output is still produced).
 */
extern int xyz_to_linear_rgb(const float* X, const float* Y, const float* Z, float* R, float* G, float* B, int width, int height, int p, int colors, float Y_white);
/*! \brief Convert XYZ planes to CIE 1976 L*u*v* planes (ISO/CIE 11664-5). */
extern int xyz_to_luv(const float* X, const float* Y, const float* Z, float* L, float* u, float* v, int width, int height, int p, int colors, float Y_white);
/*! \brief Convert CIE 1976 L*u*v* planes to XYZ planes. Returns DFX_INVRANGE for colors outside of the spectral locus. */
extern int luv_to_xyz(const float* L, const float* u, const float* v, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white);
/*! \brief Convert XYZ planes to CIE 1976 L*a*b* planes (ISO/CIE 11664-4). */
extern int xyz_to_lab(const float* X, const float* Y, const float* Z, float* L, float* a, float* b, int width, int height, int p, int colors, float Y_white);
/*! \brief Convert CIE 1976 L*a*b* planes to XYZ planes. */
extern int lab_to_xyz(const float* L, const float* a, const float* b, float* X, float* Y, float* Z, int width, int height, int p, int colors, float Y_white);

/*!
 *  \brief Compute the luminance (the Y component of XYZ) of a linear RGB image.
 *  \return DFX_INVRANGE (and no output) if any RGB value is outside of [0,1].
 */
extern int linear_rgb_to_luminance(const float* R, const float* G, const float* B, float* Y, int width, int height, int p, int colors, float Y_white);
/*! \brief Same as linear_rgb_to_luminance() with sRGB primaries and Y_white = 1. */
extern int linear_to_luminance(const float* R, const float* G, const float* B, float* Y, int width, int height, int p);
/*! \brief Produce a gray linear RGB image (R = G = B = Y) from a luminance plane normalized to Y_white = 1. */
extern int luminance_to_grayscale_image(const float* Y, float* R, float* G, float* B, int width, int height, int p);

/********************
 *
 *  Color adjustments (dfx_colors.c)
 *
 *  These operate on sRGB-primaries linear images (COLORS_SRGB, Y_white = 1)
 *  by converting to CIE L*u*v*, modifying L* or (u*,v*), and converting back.
 *  Parameter k selects how values are kept inside their valid range:
 *  k = 0 clips hard, k >= 1 uses the smooth compander
 *  f(x) = (max^-k + x^-k)^(-1/k), which approaches hard clipping as k grows.
 */

/*!
 *  \brief Scale the lightness L* of every pixel by the factor brightness (0..4) at constant chromaticity.
 *  \return DFX_INVRANGE if some adjusted colors fall outside of the RGB gamut (the output is still produced, with out-of-range values).
 */
extern int adjust_brightness(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float brightness, int k);
/*!
 *  \brief Scale the CIE 1976 u,v saturation of every pixel by the factor saturation (0..4) at constant lightness,
 *  limiting each color to the chromaticity gamut of the RGB primaries.
 */
extern int adjust_saturation(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float saturation, int k);
/*!
 *  \brief Rotate the hue of every pixel by the angle theta (radians, -pi..pi) in the (u*,v*) plane at constant lightness,
 *  desaturating colors that leave the chromaticity gamut of the RGB primaries.
 */
extern int adjust_hue(const float* R_in, const float* G_in, const float* B_in, float* R_out, float* G_out, float* B_out, int width, int height, int p, float theta);

/********************
 *
 *  Signal-to-noise metrics (dfx_snr.c)
 */

/*!
 *  \brief Compute the mean squared error and the signal-to-noise ratio between two planes.
 *  \param[in]  X_ref   - reference plane
 *  \param[in]  X_test  - test plane
 *  \param[out] p_MSE   - receives MSE = mean of (X_ref - X_test)^2 over the image interior
 *  \param[out] p_SNR   - receives SNR = 10 log10(mean of X_ref^2 / MSE) in dB, capped at 100 dB
 */
extern int snr_plane(const float* X_ref, const float* X_test, int width, int height, int p, float* p_MSE, float* p_SNR);
/*! \brief Compute MSE and SNR between two images, pooling the squared errors and energies of the three planes. */
extern int snr_image(const float* R_ref, const float* G_ref, const float* B_ref, const float* R_test, const float* G_test, const float* B_test, int width, int height, int p, float* p_MSE, float* p_SNR);

/********************
 *
 *  Binary images (dfx_binary.c)
 *
 *  A binary image for a width x height image with padding p is an array of
 *  (width + 2p) x (height + 2p) unsigned chars holding 0 (OFF) or 1 (ON),
 *  the pixel (x,y) being stored at index (p + y) * (width + 2p) + (p + x),
 *  exactly like a float plane. Set operations and thresholding process the
 *  padded area as well. The morphological operations read the padding of
 *  their input as the outside of the image (pad_binary() with PAD_ZERO makes
 *  the outside background; PAD_REPLICATE, the convention of all outputs,
 *  continues the image beyond its boundary) and leave their output
 *  replicate-padded. Operations with a structuring element of radius n
 *  require p >= n, like filters.
 */

/*! \brief Allocate a binary image for a width x height image with padding p; free it with free_binary(). */
extern int alloc_binary(unsigned char** pA, int width, int height, int p);
/*! \brief Deallocate a binary image. */
extern int free_binary(unsigned char* A);
/*! \brief Fill a binary image (padding included) with 0 or 1. */
extern int fill_binary(unsigned char* A, int width, int height, int p, int value);
/*! \brief Fill the padding border of a binary image by one of the PAD_* rules (see pad_plane()). */
extern int pad_binary(unsigned char* A, int width, int height, int p, int padding_type);
/*! \brief Copy a binary image (padding included). */
extern int copy_binary(const unsigned char* A_in, unsigned char* A_out, int width, int height, int p);

/* Set operations process the padded area as well; their output may be the same as one of their inputs. */
/*! \brief Complement: A_out = NOT A. */
extern int not_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p);
/*! \brief Intersection: A_out = A AND B. */
extern int and_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p);
/*! \brief Union: A_out = A OR B. */
extern int or_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p);
/*! \brief Symmetric difference: A_out = A XOR B. */
extern int xor_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p);
/*! \brief Set difference: A_out = A AND NOT B. */
extern int subtract_binary(const unsigned char* A, const unsigned char* B, unsigned char* A_out, int width, int height, int p);
/*! \brief Count the ON pixels in the interior of a binary image (the area of the foreground). */
extern int count_binary(const unsigned char* A, int width, int height, int p, int* p_count);

/*!
 *  \brief Convert a binary image to a plane of 0.0 and 1.0 values (padding included).
 *  Interpreted as luminance (see luminance_to_grayscale_image()), the plane displays the image in black and white.
 */
extern int binary_to_plane(const unsigned char* A, float* X, int width, int height, int p);
/*! \brief Threshold a plane: A = (X > tau), padding included. */
extern int threshold_plane(const float* X, unsigned char* A, int width, int height, int p, float tau);
/*!
 *  \brief Compute the histogram of the interior of a plane: nbins equal bins over [lo, hi), with out-of-range
 *  values counted in the end bins. hist must hold nbins ints.
 */
extern int histogram_plane(const float* X, int width, int height, int p, float lo, float hi, int nbins, int* hist);
/*!
 *  \brief Find a global threshold by Otsu's method (maximum between-class variance of the histogram).
 *
 *  The plane is quantized into nbins levels over [lo, hi] (for luminance: 0, 1, and e.g. 256 levels). The threshold
 *  returned in *p_tau separates the two classes of levels: threshold_plane(X, A, ..., tau) turns exactly the pixels of
 *  the upper class ON (tau is the largest value below the first level of that class).
 */
extern int otsu_threshold(const float* X, int width, int height, int p, float lo, float hi, int nbins, float* p_tau);
/*!
 *  \brief Find a global threshold by the Ridler-Calvard iterative method: tau <- (mean below tau + mean above tau) / 2
 *  until convergence. Same parameters as otsu_threshold().
 */
extern int ridler_calvard_threshold(const float* X, int width, int height, int p, float lo, float hi, int nbins, float* p_tau);
/*!
 *  \brief Adaptive thresholding by Niblack's method: A = (X > mu - k sigma), where mu and sigma are the mean
 *  and standard deviation of the (2n+1) x (2n+1) neighborhood of each pixel (typically k = 0.2).
 *  The input must be padded with p >= n >= 1; the output is replicate-padded.
 */
extern int niblack_threshold_plane(const float* X, unsigned char* A, int width, int height, int p, int n, float k);
/*!
 *  \brief Double thresholding (hysteresis): the connected components of (X > tau_lo) that contain a pixel of (X > tau_hi).
 *  \param[in] connectivity - 4 or 8
 */
extern int double_threshold_plane(const float* X, unsigned char* A, int width, int height, int p, float tau_lo, float tau_hi, int connectivity);

/*!
 *  \brief Generate a structuring element of radius n and the given shape (SE_CROSS, SE_BOX, SE_DISK).
 *  \param[out] B - array of (2n+1) x (2n+1) values; B[(j + n) * (2n+1) + (i + n)] is the value at offset (i,j)
 */
extern int structuring_element(unsigned char* B, int n, int shape);
/*!
 *  \brief Dilation: A_out(z) = OR over b in B of A(z - b).
 *
 *  \param[in]  A      - input binary image, padded with p >= n
 *  \param[out] A_out  - output binary image (may be the same as A); receives PAD_REPLICATE padding
 *  \param[in]  B      - structuring element of (2n+1) x (2n+1) values 0/1 (see structuring_element())
 *  \param[in]  n      - radius of the structuring element, 0 <= n <= p
 */
extern int dilate_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n);
/*! \brief Erosion: A_out(z) = AND over b in B of A(z + b). Same parameters as dilate_binary(). */
extern int erode_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n);
/*! \brief Opening: erosion followed by dilation (removes foreground details smaller than B). Same parameters as dilate_binary(). */
extern int open_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n);
/*! \brief Closing: dilation followed by erosion (fills background details smaller than B). Same parameters as dilate_binary(). */
extern int close_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n);
/*!
 *  \brief Majority (median) filter over the (2n+1) x (2n+1) neighborhood: A_out = 1 where more than half of the neighbors are ON.
 *  Requires p >= n; the output is replicate-padded and may be the same as the input.
 */
extern int majority_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n);
/*!
 *  \brief Hit-or-miss transform: A_out(z) = 1 where the ternary pattern B matches the image around z.
 *  Requires p >= n; the output is replicate-padded and may be the same as the input.
 *  \param[in] B - (2n+1) x (2n+1) values SE_ON (pixel must be ON), SE_OFF (must be OFF), SE_DONTCARE
 */
extern int hitmiss_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, const unsigned char* B, int n);
/*!
 *  \brief Thinning: repeatedly remove the pixels matched by the eight thinning structuring elements (Birchfield, Fig. 4.14).
 *
 *  Iterated to convergence, thinning yields an approximate skeleton. Requires p >= 1; the padding of the input decides
 *  what lies outside the image for the whole iteration (PAD_ZERO: background). The output is replicate-padded.
 *  \param[in]  n_iter - number of passes over the eight structuring elements, or <= 0 to iterate until convergence
 *  \param[out] p_iter - if not NULL, receives the number of passes performed
 */
extern int thin_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n_iter, int* p_iter);
/*!
 *  \brief Thickening: repeatedly add the background pixels matched by the eight thickening structuring elements
 *  (Birchfield, Fig. 4.17). Iterated to convergence, thickening approximates the convex hull. Same parameters as thin_binary().
 */
extern int thicken_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int n_iter, int* p_iter);
/*!
 *  \brief Inner boundary of the foreground: A minus its erosion by a 3 x 3 structuring element. Requires p >= 1.
 *  \param[in] connectivity - 8: the 8-connected boundary (erosion by B4); 4: the thicker 4-connected boundary (erosion by B8)
 */
extern int boundary_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int connectivity);

/********************
 *
 *  Regions of binary images (dfx_regions.c)
 *
 *  Connected components are numbered 1, 2, 3, ... in order of first
 *  appearance in row-by-row scanning and stored in label planes: arrays of
 *  ints with the geometry of a binary image, holding 0 for the background
 *  (and in the padding). The parameter "connectivity" is 4 or 8 throughout.
 */

/*! \brief Allocate a label plane; free it with free_labels(). */
extern int alloc_labels(int** pL, int width, int height, int p);
/*! \brief Deallocate a label plane. */
extern int free_labels(int* L);
/*! \brief Extract one region of a label plane as a binary image: A = (L == label); A is replicate-padded. */
extern int select_label(const int* L, unsigned char* A, int width, int height, int p, int label);
/*! \brief Render a label plane as a pseudocolored linear RGB image (black background, one hue per label). */
extern int labels_to_image(const int* L, float* R, float* G, float* B, int width, int height, int p);

/*!
 *  \brief Floodfill: replace the value of the connected region of L containing (x0,y0) by new_label (Birchfield, Algorithm 4.4).
 *  \param[out] p_count - if not NULL, receives the number of pixels filled (0 if the seed already holds new_label)
 */
extern int floodfill_labels(int* L, int width, int height, int p, int x0, int y0, int new_label, int connectivity, int* p_count);
/*!
 *  \brief Label the connected components of a binary image by the two-pass union-find algorithm (Birchfield, Algorithm 4.8).
 *  \param[out] L          - label plane: 0 for background pixels, the component number otherwise
 *  \param[out] p_n_labels - if not NULL, receives the number of components
 */
extern int label_binary(const unsigned char* A, int* L, int width, int height, int p, int connectivity, int* p_n_labels);
/*! \brief Label the connected components by repeated floodfills; gives the same labels as label_binary(). */
extern int label_binary_floodfill(const unsigned char* A, int* L, int width, int height, int p, int connectivity, int* p_n_labels);
/*! \brief Morphological reconstruction: the connected components of mask that contain an ON pixel of seeds. */
extern int reconstruct_binary(const unsigned char* mask, const unsigned char* seeds, unsigned char* A_out, int width, int height, int p, int connectivity);
/*! \brief Fill the holes of the foreground (the background components not connected to the image boundary); connectivity refers to the background. */
extern int fill_holes_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int connectivity);

/*!
 *  \brief Trace the boundary of the region containing (x0,y0) clockwise and return it as a Freeman chain code.
 *
 *  Codes 0..7 denote steps east, north-east, north (up), north-west, west, south-west, south, south-east.
 *  connectivity = 8 traces the 8-connected boundary (Moore's neighbor tracing, Jacob's stopping criterion),
 *  connectivity = 4 the 4-connected boundary (wall following). Pixels outside of the image count as OFF.
 *  The trace follows the boundary between the region and the background containing the OFF pixel above the start:
 *  from the first pixel of a region in scanning order, that is the outer boundary; from a pixel below a hole, the hole's.
 *
 *  \param[in]  x0, y0     - starting pixel: ON, with an OFF pixel or the image border above it
 *  \param[out] chain      - array for the codes, or NULL to only count them
 *  \param[in]  max_length - capacity of chain
 *  \param[out] p_length   - receives the number of codes (0 for an isolated pixel)
 *  \return DFX_INVARG for an invalid starting pixel, DFX_NOMEM if the chain does not fit into max_length codes.
 */
extern int trace_boundary_binary(const unsigned char* A, int width, int height, int p, int x0, int y0, int connectivity, unsigned char* chain, int max_length, int* p_length);
/*! \brief Estimate the length of a boundary from its chain code by CHAIN_FREEMAN, CHAIN_PYTHAGOREAN, or CHAIN_KIMURA. */
extern int chain_code_length(const unsigned char* chain, int length, int method, double* p_length);

/*!
 *  \brief Compute the properties of region "label" of a label plane: area, bounding box, centroid, central moments,
 *  orientation, eccentricity and semi-axes of the equivalent ellipse, boundary pixel count, boundary length, compactness.
 *  \return DFX_INVRANGE if no pixel carries the label.
 */
extern int region_properties(const int* L, int width, int height, int p, int label, dfx_region* r);
/*! \brief Compute the properties of the whole foreground of a binary image, treated as one region. */
extern int region_properties_binary(const unsigned char* A, int width, int height, int p, dfx_region* r);

/*!
 *  \brief Chamfer distance transform (Birchfield, Algorithm 4.10): D = distance from each pixel to the nearest ON pixel,
 *  with step costs a (horizontal/vertical) and b (diagonal): (1,1) chessboard, b >= 2a city-block, (1,1.35) quasi-Euclidean.
 *  \return DFX_INVRANGE if the image has no ON pixel.
 */
extern int chamfer_binary(const unsigned char* A, float* D, int width, int height, int p, float a, float b);
/*!
 *  \brief Seeded region growing: starting from (x0,y0), add neighboring pixels whose values stay within tolerance
 *  of the running mean of the region. A receives the region grown.
 */
extern int grow_region_plane(const float* X, unsigned char* A, int width, int height, int p, int x0, int y0, float tolerance, int connectivity);

#ifdef __cplusplus
}
#endif

#endif /* DFX_H */
