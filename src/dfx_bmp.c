/*!
 *  \file   dfx_bmp.c
 *  \brief  I/O operations with bitmap files and sRGB images.
 *
 *  This module implements reading and writing images from/to bitmap (.bmp) files.
 *  The implementation is limited to 3-channel, 24-bit packed sRGB images.
 *  Attempts to open bitmap files in any other formats will result in DFX_NOTSUP error.
 *
 *  In memory, an sRGB image follows the layout of the pixel array in a bitmap file:
 *  B,G,R byte triplets, rows padded to a multiple of 4 bytes, rows stored bottom-up.
 *
 *  Bitmap headers are read and written directly as packed C structures. The fields of a
 *  bitmap file are little-endian, so this only works on little-endian hosts, which today
 *  means every x86 and ARM desktop, laptop, and phone. Big-endian hosts (IBM z, AIX on
 *  POWER, SPARC, some embedded MIPS/ARM) are rejected at compile time.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    November 8, 2026
 */

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dfx.h"

/* refuse to build on big-endian hosts (these macros are only defined when the target is big-endian): */
#if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) || defined(__BIG_ENDIAN__) || defined(__ARMEB__) || defined(__MIPSEB__)
#error "dfx_bmp.c reads bitmap headers as little-endian structures and does not support big-endian hosts"
#endif

#pragma pack(push, 2)

/*! BITMAPFILEHEADER structure (adapted from Windows / WinGDI.h) */
typedef struct {
	unsigned short bfType;			/* File type; must be "BM" (0x4D42 in hexadecimal) */
	unsigned int   bfSize;			/* Size of the BMP file in bytes */
	unsigned short bfReserved1;		/* Reserved; must be zero */
	unsigned short bfReserved2;		/* Reserved; must be zero */
	unsigned int   bfOffBits;		/* Offset to the start of the pixel data */
} BITMAPFILEHEADER;

/*! BITMAPINFOHEADER structure (adapted from Windows / WinGDI.h) */
typedef struct
{
	unsigned int   biSize;			/* Specifies the number of bytes required by the structure. */
	int            biWidth;			/* Specifies the width of the bitmap, in pixels. */
	int            biHeight;		/* Specifies the height of the bitmap, in pixels. Negative for top-down bitmaps. */
	unsigned short biPlanes;		/* Specifies the number of planes for the target device. Must be set to 1. */
	unsigned short biBitCount;		/* Specifies the number of bits per pixel (bpp). */
	unsigned int   biCompression;	/* Image format. BI_RGB implies uncompressed RGB images. */
	unsigned int   biSizeImage;		/* Specifies the size, in bytes, of the image. */
	int            biXPelsPerMeter;	/* Specifies the horizontal resolution, in pixels per meter, of the target device for the bitmap. */
	int            biYPelsPerMeter; /* Specifies the vertical resolution, in pixels per meter, of the target device for the bitmap. */
	unsigned int   biClrUsed;		/* Specifies the number of color indices in the color palette/table. */
	unsigned int   biClrImportant;	/* Specifies the number of color indices that are considered important for displaying the bitmap. */
} BITMAPINFOHEADER;

#pragma pack(pop)

/* compile-time checks that the compiler packed the structures as the file format requires (14 and 40 bytes): */
typedef char dfx_check_bitmapfileheader_size[(sizeof(BITMAPFILEHEADER) == 14) ? 1 : -1];
typedef char dfx_check_bitmapinfoheader_size[(sizeof(BITMAPINFOHEADER) == 40) ? 1 : -1];

/* bitmap-related constants */
#define BF_TYPE            0x4D42	/* "BM" */
#define BI_RGB             0		/* biCompression value for uncompressed RGB images */
#define BMP_MAX_DIMENSION  32768	/* largest width or height we agree to work with */

/*!
 *  \brief Compute the number of bytes in one row of a 24-bit bitmap image (rows are padded to multiples of 4 bytes).
 */
static unsigned int srgb_row_size(int width)
{
	return ((unsigned int)width * 3 + 3) & ~3u;
}

/*!
 *  \brief Compute the size of a 24-bit bitmap image in bytes.
 */
static unsigned int srgb_image_size(int width, int height)
{
	return (unsigned int)height * srgb_row_size(width);
}

/*!
 *  \brief Allocate an sRGB image.
 *
 *  \param[out]  p_srgb   - pointer to a variable to receive the pointer to the allocated image
 *  \param[in]   width    - image width in pixels
 *  \param[in]   height   - image height in pixels
 *
 *  \return     DFX_X error code
 */
int alloc_srgb_image(unsigned char** p_srgb, int width, int height)
{
	unsigned char* srgb = NULL;
	unsigned int size;

	/* check parameters: */
	if (p_srgb == NULL || width < 0 || height < 0 || width > BMP_MAX_DIMENSION || height > BMP_MAX_DIMENSION) return DFX_INVARG;

	/* allocate new image: */
	size = srgb_image_size(width, height);
	if ((srgb = (unsigned char*)malloc(size > 0 ? size : 1)) == NULL) return DFX_NOMEM;	/* malloc(0) may return NULL */

	/* return pointer & exit: */
	*p_srgb = srgb;
	return DFX_SUCCESS;
}

/*!
 *  \brief Deallocate sRGB image.
 *
 *  \param[in]  srgb   - pointer to the allocated sRGB image
 *
 *  \return     DFX_X error code
 */
int free_srgb_image(unsigned char* srgb)
{
	if (srgb != NULL) free(srgb);
	return DFX_SUCCESS;
}

/*!
 *  \brief Read 24-bit packed sRGB image from a bitmap file.
 *
 *  This function opens the file, reads bitmap parameters, allocates memory to hold an sRGB image,
 *  reads this image, and returns a pointer to the allocated image. It also returns image width, height,
 *  and pixel density (pixels per meter) in horizontal and vertical directions.
 *  The caller is responsible for freeing the allocated image when it is no longer needed.
 *
 *  Files with BITMAPINFOHEADER or any of its later extensions (BITMAPV4HEADER, BITMAPV5HEADER) are
 *  accepted, as long as the image is 24-bit uncompressed RGB. Both bottom-up (positive height) and
 *  top-down (negative height) bitmaps are read; the image is always returned in bottom-up order.
 *
 *  \param[in]   filename - name of a bitmap file to read
 *  \param[out]  p_srgb   - pointer to a variable to receive the pointer to the allocated sRGB image
 *  \param[out]  p_width  - pointer to a variable to receive image width
 *  \param[out]  p_height - pointer to a variable to receive image height
 *  \param[out]  p_ppm_x  - pointer to a variable to receive pixel density (pixels per meter) in horizontal direction
 *  \param[out]  p_ppm_y  - pointer to a variable to receive pixel density (pixels per meter) in vertical direction
 *
 *  \return         DFX_X error code
 */
int read_bitmap(const char* filename, unsigned char** p_srgb, int* p_width, int* p_height, int* p_ppm_x, int* p_ppm_y)
{
	BITMAPFILEHEADER hdr;       /* bitmap file header */
	BITMAPINFOHEADER bihdr;     /* bitmap info header (first 40 bytes of it, if the file uses a longer variant) */
	int width, height, top_down, y;
	unsigned int row_size, size;
	unsigned char* srgb;
	FILE* fp;

	/* check parameters: */
	if (filename == NULL || p_srgb == NULL || p_width == NULL || p_height == NULL || p_ppm_x == NULL || p_ppm_y == NULL)
		return DFX_INVARG;

	/* open file: */
	if ((fp = fopen(filename, "rb")) == NULL)
		return DFX_CANTOPEN;

	/* read headers: */
	if (fread(&hdr, 1, sizeof(hdr), fp) != sizeof(hdr)
	 || fread(&bihdr, 1, sizeof(bihdr), fp) != sizeof(bihdr)) {
		fclose(fp);
		return DFX_IOERR;
	}
	width = bihdr.biWidth;
	height = bihdr.biHeight;

	/* check if we can work with this file (the height range is checked before it is negated, since -INT_MIN overflows): */
	if (hdr.bfType != BF_TYPE
	 || bihdr.biSize < sizeof(BITMAPINFOHEADER)
	 || hdr.bfOffBits < sizeof(BITMAPFILEHEADER) + bihdr.biSize
	 || bihdr.biPlanes != 1
	 || bihdr.biCompression != BI_RGB || bihdr.biBitCount != 24	/* we only support 24-bit RGB files now */
	 || width < 1 || width > BMP_MAX_DIMENSION
	 || height == 0 || height < -BMP_MAX_DIMENSION || height > BMP_MAX_DIMENSION) {
		fclose(fp);
		return DFX_NOTSUP;
	}
	top_down = (height < 0);
	if (top_down) height = -height;

	/* compute image size and check it against the header, if specified: */
	row_size = srgb_row_size(width);
	size = srgb_image_size(width, height);
	if (bihdr.biSizeImage != 0 && bihdr.biSizeImage < size) {
		fclose(fp);
		return DFX_NOTSUP;
	}

	/* allocate image: */
	if (alloc_srgb_image(&srgb, width, height) != DFX_SUCCESS) {
		fclose(fp);
		return DFX_NOMEM;
	}

	/* move to the pixel array (skipping the rest of a longer info header, color masks, or a color profile): */
	if (fseek(fp, (long)hdr.bfOffBits, SEEK_SET) != 0) {
		fclose(fp);
		free_srgb_image(srgb);
		return DFX_IOERR;
	}

	/* read the pixel array; a top-down file is stored row-reversed to keep the bottom-up layout in memory: */
	for (y = 0; y < height; y++) {
		unsigned char* row = srgb + (unsigned int)(top_down ? height - 1 - y : y) * row_size;
		if (fread(row, 1, row_size, fp) != row_size) {
			fclose(fp);
			free_srgb_image(srgb);
			return DFX_IOERR;
		}
	}

	/* return pointer to the allocated image: */
	*p_srgb = srgb;

	/* return image parameters: */
	*p_width  = width;
	*p_height = height;
	*p_ppm_x  = bihdr.biXPelsPerMeter;
	*p_ppm_y  = bihdr.biYPelsPerMeter;

	/* close file and exit: */
	fclose(fp);
	return DFX_SUCCESS;
}

/*!
 *  \brief Write 24-bit packed sRGB image into a bitmap file.
 *
 *  \param[in]  filename - name of a bitmap file to be created
 *  \param[in]  srgb     - pointer to the sRGB image
 *  \param[in]  width    - image width in pixels
 *  \param[in]  height   - image height in pixels
 *  \param[in]  ppm_x    - pixel density (pixels per meter) in horizontal direction (0 = unspecified)
 *  \param[in]  ppm_y    - pixel density (pixels per meter) in vertical direction (0 = unspecified)
 *
 *  \return     DFX_X error code
 */
int write_bitmap(const char* filename, const unsigned char* srgb, int width, int height, int ppm_x, int ppm_y)
{
	BITMAPFILEHEADER hdr;       /* bitmap file header */
	BITMAPINFOHEADER bihdr;     /* bitmap info header */
	static const unsigned char zeros[4] = {0, 0, 0, 0};	/* row padding bytes */
	unsigned int row_size, row_pad, size;
	int y;
	FILE* fp;

	/* check parameters: */
	if (filename == NULL || srgb == NULL || width < 1 || height < 1 || width > BMP_MAX_DIMENSION || height > BMP_MAX_DIMENSION || ppm_x < 0 || ppm_y < 0)
		return DFX_INVARG;

	/* compute image size: */
	row_size = srgb_row_size(width);
	row_pad = row_size - (unsigned int)width * 3;
	size = srgb_image_size(width, height);

	/* fill file header:  */
	memset(&hdr, 0, sizeof(hdr));
	hdr.bfType = BF_TYPE;
	hdr.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + size;
	hdr.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

	/* fill bitmap header:  */
	memset(&bihdr, 0, sizeof(bihdr));
	bihdr.biSize = sizeof(BITMAPINFOHEADER);
	bihdr.biWidth = width;
	bihdr.biHeight = height;
	bihdr.biPlanes = 1;
	bihdr.biBitCount = 24;
	bihdr.biCompression = BI_RGB;
	bihdr.biSizeImage = size;
	bihdr.biXPelsPerMeter = ppm_x;
	bihdr.biYPelsPerMeter = ppm_y;

	/* create file: */
	if ((fp = fopen(filename, "wb")) == NULL)
		return DFX_CANTOPEN;

	/* write headers: */
	if (fwrite(&hdr, 1, sizeof(hdr), fp) != sizeof(hdr)
	 || fwrite(&bihdr, 1, sizeof(bihdr), fp) != sizeof(bihdr)) {
		fclose(fp);
		return DFX_IOERR;
	}

	/* write the pixel array row by row, supplying the padding bytes explicitly: */
	for (y = 0; y < height; y++) {
		if (fwrite(srgb + (unsigned int)y * row_size, 1, (size_t)width * 3, fp) != (size_t)width * 3
		 || (row_pad > 0 && fwrite(zeros, 1, row_pad, fp) != row_pad)) {
			fclose(fp);
			return DFX_IOERR;
		}
	}

	/* close file and exit: */
	if (fclose(fp) != 0)
		return DFX_IOERR;
	return DFX_SUCCESS;
}

/* dfx_bmp.c -- end of file */
