/*!
 *  \file   dfx_bmp.c
 *  \brief  DFX library: I/O operations with bitmap files and sRGB images.
 * 
 *  This module implements reading and writing images from/to bitmap (.bmp) files. 
 *  The implementation is limited to 3-channel, 24-bit packed sRGB images. 
 *  Attempts to open bitmap files in any other formats will result in DFX_NOTSUP error.
 * 
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.01
 *  \date    March 7, 2026
 */

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dfx.h"

/*! BITMAPFILEHEADER structure (adapted from Windows / WinGDI.h) */
#pragma pack(2)
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
	int            biHeight;		/* Specifies the height of the bitmap, in pixels. */
	unsigned short biPlanes;		/* Specifies the number of planes for the target device. Must be set to 1. */
	unsigned short biBitCount;		/* Specifies the number of bits per pixel (bpp). */
	unsigned int   biCompression;	/* Image format. BI_RGB implies uncompressed RGB images. */
	unsigned int   biSizeImage;		/* Specifies the size, in bytes, of the image. */
	int            biXPelsPerMeter;	/* Specifies the horizontal resolution, in pixels per meter, of the target device for the bitmap. */
	int            biYPelsPerMeter; /* Specifies the vertical resolution, in pixels per meter, of the target device for the bitmap. */
	unsigned int   biClrUsed;		/* Specifies the number of color indices in the color palette/table. */
	unsigned int   biClrImportant;	/* Specifies the number of color indices that are considered important for displaying the bitmap. */
} BITMAPINFOHEADER;

/* bitmap-related constants */
#define BF_TYPE 0x4D42  
#define BI_RGB  0   

/*!
 * \brief Compute size of a 24-bit bitmap image:
 */
static unsigned int srgb_image_size(int width, int height)
{
	unsigned int w_srgb = (width * 3 + 3) & (~3);		/* w_srgb = offset to the next line in 24-bit packed bitmap image */
	unsigned int size = height * w_srgb;				/* bitmap image size */
	return size;
}

/*!
 *  \brief Allocate an sRGB image.
 *
 *  \param[in,out]  p_srgb   - pointer to pointer to the allocated sRGB image
 *
 *  \return     DFX_X error codes
 */
int alloc_srgb_image(unsigned char** p_srgb, int width, int height)
{
	int size = srgb_image_size(width, height);
	unsigned char* srgb = NULL;
	if (p_srgb == NULL || width < 0 || height < 0) return DFX_INVARG;
	if ((srgb = (unsigned char*)malloc(size)) == NULL) return DFX_NOMEM;
	*p_srgb = srgb;
	return DFX_SUCCESS;
}

/*!
 *  \brief Deallocate sRGB image.
 *
 *  \param[in]  srgb   - pointer to the allocated sRGB image
 *
 *  \return     DFX_X error codes
 */
int free_srgb_image(unsigned char* srgb)
{
	if (srgb != NULL) free(srgb);
	return DFX_SUCCESS;
}

/*!
 *  \brief Read 24-bit packed sRGB image from a bitmap file.
 *
 *  This function opens file, reads bitmap parameters, allocates memory to hold an sRGB image, 
 *  reads this image, and returns a pointer to the allocated image. It also returns image width, height, 
 *  and pixel density (pixels per meter) in horisontal and vertical directions. 
 *  The caller is responsible for freeing the allocated image when it is no longer needed. 
 *
 *  \param[in]      filename - name of a bitmap file to read
 *  \param[in,out]  p_srgb   - pointer to a pointer to the allocated sRGB image 
 *  \param[in,out]  p_width  - pointer to a variable storing image width
 *  \param[in,out]  p_height - pointer to a variable storing image height
 *  \param[in,out]  p_ppm_x  - pointer to a variable storing pixel density (pixels per meter) in horisontal direction
 *  \param[in,out]  p_ppm_y  - pointer to a variable storing pixel density (pixels per meter) in vertical direction
 *
 *  \return         DFX_X error codes
 */
int read_bitmap(char* filename, unsigned char** p_srgb, int* p_width, int* p_height, int* p_ppm_x, int* p_ppm_y)
{
	BITMAPFILEHEADER hdr;       /* bitmap file header */
	BITMAPINFOHEADER bihdr;     /* bitmap info header */
	unsigned int size;
	unsigned char* srgb;
	FILE* fp;

	/* check parameters: */
	if (filename == NULL || p_srgb == NULL || p_width == NULL || p_height == NULL || p_ppm_x == NULL || p_ppm_y == NULL)
		return DFX_INVARG;

	/* open file: */
	if ((fp = fopen(filename, "rb")) == NULL) 
		return DFX_CANTOPEN;

	/* read headers */
	if (fread((void*)&hdr, 1, sizeof(hdr), fp) != sizeof(hdr)
     || fread((void*)&bihdr, 1, sizeof(bihdr), fp) != sizeof(bihdr)) {
		fclose(fp);
		return DFX_IOERR;
	}

	/* compute image size: */
	size = srgb_image_size(bihdr.biWidth, bihdr.biHeight);

	/* check if we can work with this file:  */
	if (hdr.bfOffBits != sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)
		|| bihdr.biCompression != BI_RGB || bihdr.biBitCount != 24 /* we only support 24-bit RGB files now */
		|| bihdr.biWidth < 8 || bihdr.biHeight < 8 || bihdr.biWidth > 32768 || bihdr.biHeight > 32768 
		|| (bihdr.biSizeImage != 0 && bihdr.biSizeImage < size)) {
		fclose(fp);
		return DFX_NOTSUP;
	}

	/* allocate image */
	if (alloc_srgb_image(&srgb, bihdr.biWidth, bihdr.biHeight) != DFX_SUCCESS) { 
		fclose(fp); 
		return DFX_NOMEM; 
	}

	/* read the bitmap image */
	if (fread((void*)srgb, 1, size, fp) != size) {
		fclose(fp);
		free_srgb_image(srgb);
		return DFX_IOERR;
	}

	/* return pointer to the allocated image: */
	*p_srgb = srgb;

	/* return image parameters */
	*p_width  = bihdr.biWidth;
	*p_height = bihdr.biHeight;
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
 *  \param[in]  srgb     - pointer to the allocated sRGB image
 *  \param[in]  p_width  - image width
 *  \param[in]  p_height - image height
 *  \param[in]  p_ppm_x  - pixel density (pixels per meter) in horisontal direction
 *  \param[in]  p_ppm_y  - pixel density (pixels per meter) in vertical direction
 *
 *  \return     DFX_X error codes
 */
int write_bitmap(char* filename, unsigned char* srgb, int width, int height, int ppm_x, int ppm_y)
{
	static BITMAPFILEHEADER hdr;       /* bitmap file header */
	static BITMAPINFOHEADER bihdr;     /* bitmap info header */
	unsigned int size;
	FILE* fp;

	/* check parameters: */
	if (filename == NULL || srgb == NULL || width < 0 || height < 0 || ppm_x <= 0 || ppm_y <= 0)
		return DFX_INVARG;

	/* compute frame size */
	size = srgb_image_size(width, height);

	/* fill file header:  */
	memset(&hdr, 0, sizeof(hdr));
	hdr.bfType = 0x4d42;
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
	if ((fp = fopen(filename, "wb+")) == NULL)
		return DFX_CANTOPEN;

	/* write headers & image: */
	if (fwrite((void*)&hdr, 1, sizeof(hdr), fp) != sizeof(hdr)
     || fwrite((void*)&bihdr, 1, sizeof(bihdr), fp) != sizeof(bihdr)
	 || fwrite((void*)srgb, 1, size, fp) != size) {
		fclose(fp);
		return DFX_IOERR;
	}

	/* close file and exit: */
	fclose(fp);
	return DFX_SUCCESS;
}

/* dfx_bmp.c -- end of file */
