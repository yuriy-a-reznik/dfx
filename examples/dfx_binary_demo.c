/*!
 *  \file   dfx_binary_demo.c
 *  \brief  An example of binary image processing: thresholding, morphology, labeling, region properties, distances.
 *
 *  Without arguments, this example synthesizes a test scene of bright objects on a dark noisy
 *  background (disks, a ring, a crescent, an ellipse, a thin bar, and a disk with a strong brightness
 *  gradient), all in linear light, and then walks through the standard binary processing chain:
 *
 *    1. histogram of the luminance; global thresholds by Otsu's and the Ridler-Calvard methods;
 *    2. clean-up of the thresholded image by a majority filter, opening, and closing; hole filling;
 *    3. morphological operations with structuring elements of various sizes, boundary extraction;
 *    4. connected components (two algorithms, which must agree), pseudocolored labels, and a table
 *       of region properties (area, centroid, orientation, eccentricity, perimeter, compactness);
 *    5. chamfer distance transforms, inside and outside the objects;
 *    6. skeletons by thinning;
 *    7. double thresholding, adaptive (Niblack) thresholding, and region growing.
 *
 *  Since the true shapes are known, the results are checked against them (areas, perimeters,
 *  orientation and eccentricity of the ellipse, the Euclidean radius recovered by the chamfer
 *  distance, etc.); the program exits with a nonzero code if any check fails. With a bitmap file as
 *  argument, the same chain runs on the luminance of that image; the checks that compare the
 *  algorithms with each other (cascaded dilations, duality, the two labeling algorithms, thinning)
 *  still run, those that need the known shapes are skipped.
 *
 *  Output files: dfx-binary-input.bmp, -threshold.bmp, -clean.bmp, -filled.bmp, -boundary.bmp,
 *  -labels.bmp, -chamfer.bmp, -skeleton.bmp, -double.bmp, -niblack.bmp, -grown.bmp.
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

/* geometry of the synthetic scene: */
#define WIDTH   640
#define HEIGHT  480
#define PADDING 16

/* the objects of the synthetic scene: */
#define N_OBJECTS 7
static const struct {
	const char* name;
	double x, y;			/* center */
	double r1, r2;			/* radii / semi-axes */
	double angle;			/* orientation, degrees */
	double level;			/* luminance */
} objects[N_OBJECTS] = {
	{ "big disk",  120, 120, 60,  0, 0,  0.75 },
	{ "ring",      330, 110, 45, 14, 0,  0.70 },
	{ "small disk",520, 100, 32,  0, 0,  0.80 },
	{ "dim disk",  250, 200, 30,  0, 0,  0.90 },	/* its luminance falls linearly from 0.9 to 0.2 from left to right */
	{ "crescent",  140, 330, 75, 70, 0,  0.65 },	/* disk of radius 75 minus the disk of radius 70 centered 20 pixels to the right and 15 up */
	{ "ellipse",   420, 330, 90, 30, 30, 0.70 },
	{ "thin bar",  560, 430, 55,  3.5, -25, 0.80 },
};

/* a small deterministic pseudo-random generator, so that the test scene is the same on every platform: */
static unsigned int lcg_state = 12345u;
static float frand(void)
{
	lcg_state = lcg_state * 1664525u + 1013904223u;
	return (float)(lcg_state >> 8) / 16777216.0f;		/* uniform in [0,1) */
}

static int checks_failed = 0;
static void check(int condition, const char* what)
{
	printf("  [%s] %s\n", condition ? " ok " : "FAIL", what);
	if (!condition) checks_failed++;
}

/*!
 *  \brief Generate the synthetic scene in a luminance plane.
 */
static void generate_scene(float* Y, int width, int height, int p)
{
	int x, y, k, i;
	int width_p = width + 2 * p;

	/* dark background with a gentle gradient and noise: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		float noise = 0.0f;
		for (i = 0; i < 3; i++) noise += frand() - 0.5f;	/* roughly Gaussian, standard deviation 0.5 */
		Y[(p + y) * width_p + p + x] = 0.06f + 0.04f * (float)x / width + 0.03f * noise;
	}

	/* the objects, with a little shading: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		for (k = 0; k < N_OBJECTS; k++) {
			double dx = x - objects[k].x, dy = y - objects[k].y, a = objects[k].angle * M_PI / 180.0;
			double u = dx * cos(a) + dy * sin(a), v = -dx * sin(a) + dy * cos(a);
			double level = objects[k].level;
			int inside;
			switch (k) {
			case 1:  inside = dx * dx + dy * dy <= objects[k].r1 * objects[k].r1 && dx * dx + dy * dy > objects[k].r2 * objects[k].r2; break;
			case 3:  inside = dx * dx + dy * dy <= objects[k].r1 * objects[k].r1;
				level = 0.9 - 0.7 * (dx + objects[k].r1) / (2 * objects[k].r1); break;
			case 4:  inside = dx * dx + dy * dy <= objects[k].r1 * objects[k].r1
				&& (dx - 20) * (dx - 20) + (dy + 15) * (dy + 15) > objects[k].r2 * objects[k].r2; break;
			case 5: case 6: inside = u * u / (objects[k].r1 * objects[k].r1) + v * v / (objects[k].r2 * objects[k].r2) <= 1.0; break;
			default: inside = dx * dx + dy * dy <= objects[k].r1 * objects[k].r1; break;
			}
			if (inside) {
				/* shading: brighter towards the upper left */
				level *= 0.85 + 0.15 * (0.5 - 0.5 * (dx + dy) / (objects[k].r1 + 1.0));
				Y[(p + y) * width_p + p + x] = (float)level + 0.01f * (frand() - 0.5f);
			}
		}
	}

	/* salt specks on the background and pepper specks on the objects: */
	for (k = 0; k < 150; k++) {
		x = (int)(frand() * width); y = (int)(frand() * height);
		if (Y[(p + y) * width_p + p + x] < 0.2f) Y[(p + y) * width_p + p + x] = 0.9f;
	}
	for (k = 0; k < 600; k++) {
		x = (int)(frand() * width); y = (int)(frand() * height);
		if (Y[(p + y) * width_p + p + x] > 0.3f) Y[(p + y) * width_p + p + x] = 0.02f;
	}
	clip_plane(Y, Y, width, height, p, 0.0f, 1.0f);
	pad_plane(Y, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Write a luminance plane as a grayscale bitmap.
 */
static int write_plane(const char* filename, const float* Y, int width, int height, int p)
{
	float* R = NULL, * G = NULL, * B = NULL;
	unsigned char* srgb = NULL;
	int err;
	if ((err = alloc_image(&R, &G, &B, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = alloc_srgb_image(&srgb, width, height)) != DFX_SUCCESS) { free_image(R, G, B); return err; }
	luminance_to_grayscale_image(Y, R, G, B, width, height, p);
	linear_to_srgb(R, G, B, srgb, width, height, p);
	err = write_bitmap(filename, srgb, width, height, 0, 0);
	free_srgb_image(srgb);
	free_image(R, G, B);
	if (err != DFX_SUCCESS) printf("Error: cannot write %s\n", filename);
	else printf("  written %s\n", filename);
	return err;
}

/*!
 *  \brief Write a binary image as a black and white bitmap.
 */
static int write_binary(const char* filename, const unsigned char* A, int width, int height, int p)
{
	float* Y = NULL;
	int err;
	if ((err = alloc_plane(&Y, width, height, p)) != DFX_SUCCESS) return err;
	binary_to_plane(A, Y, width, height, p);
	err = write_plane(filename, Y, width, height, p);
	free_plane(Y);
	return err;
}

/*!
 *  \brief Write a linear RGB image as a bitmap.
 */
static int write_rgb(const char* filename, const float* R, const float* G, const float* B, int width, int height, int p)
{
	unsigned char* srgb = NULL;
	int err;
	if ((err = alloc_srgb_image(&srgb, width, height)) != DFX_SUCCESS) return err;
	linear_to_srgb(R, G, B, srgb, width, height, p);
	err = write_bitmap(filename, srgb, width, height, 0, 0);
	free_srgb_image(srgb);
	if (err != DFX_SUCCESS) printf("Error: cannot write %s\n", filename);
	else printf("  written %s\n", filename);
	return err;
}

/*!
 *  \brief Print a histogram as a bar chart.
 */
static void print_histogram(const int* hist, int nbins, float lo, float hi)
{
	int k, j, max = 1, bar;
	for (k = 0; k < nbins; k++) if (hist[k] > max) max = hist[k];
	for (k = 0; k < nbins; k++) {
		bar = (int)(50.0 * hist[k] / max + 0.5);
		printf("  %5.3f-%5.3f %7d |", lo + (hi - lo) * k / nbins, lo + (hi - lo) * (k + 1) / nbins, hist[k]);
		for (j = 0; j < bar; j++) putchar('#');
		putchar('\n');
	}
}

/*!
 *  \brief Find the label of the region whose centroid is nearest to (x,y); returns 0 if none is within 30 pixels.
 */
static int nearest_region(const dfx_region* r, int n, double x, double y)
{
	int k, best = 0;
	double d, d_best = 30.0 * 30.0;
	for (k = 0; k < n; k++) {
		d = (r[k].x_c - x) * (r[k].x_c - x) + (r[k].y_c - y) * (r[k].y_c - y);
		if (d < d_best) { d_best = d; best = k + 1; }
	}
	return best;
}

int main(int argc, char* argv[])
{
	int width = WIDTH, height = HEIGHT, p = PADDING, synthetic = (argc < 2);
	int x, y, k, n, err, count, count2, n_labels, n_labels2, iter, hist[16];
	int width_p;
	float* Y = NULL, * D = NULL, * R = NULL, * G = NULL, * B = NULL;
	unsigned char* A = NULL, * C = NULL, * T = NULL, * S = NULL, se[9], se_big[121];
	int* L = NULL, * L2 = NULL;
	dfx_region* regions = NULL;
	float tau_otsu, tau_rc, d_max;
	double length;

	/* get the input: a bitmap file, or the synthetic scene */
	if (!synthetic) {
		unsigned char* srgb = NULL;
		int ppm_x, ppm_y;
		if ((err = read_bitmap(argv[1], &srgb, &width, &height, &ppm_x, &ppm_y)) != DFX_SUCCESS) {
			printf("Error: cannot read bitmap %s (error %d)\n", argv[1], err);
			return 1;
		}
		if (alloc_image(&R, &G, &B, width, height, p) != DFX_SUCCESS || alloc_plane(&Y, width, height, p) != DFX_SUCCESS) {
			printf("Error: cannot allocate memory\n"); return 1;
		}
		srgb_to_linear(srgb, R, G, B, width, height, p);
		linear_to_luminance(R, G, B, Y, width, height, p);
		free_srgb_image(srgb);
		free_image(R, G, B); R = G = B = NULL;
		printf("Input: %s, %d x %d pixels\n", argv[1], width, height);
	} else {
		if (alloc_plane(&Y, width, height, p) != DFX_SUCCESS) { printf("Error: cannot allocate memory\n"); return 1; }
		generate_scene(Y, width, height, p);
		printf("Input: synthetic scene, %d x %d pixels, %d objects\n", width, height, N_OBJECTS);
	}
	width_p = width + 2 * p;

	/* allocate working images: */
	if (alloc_plane(&D, width, height, p) != DFX_SUCCESS || alloc_image(&R, &G, &B, width, height, p) != DFX_SUCCESS
		|| alloc_binary(&A, width, height, p) != DFX_SUCCESS || alloc_binary(&C, width, height, p) != DFX_SUCCESS
		|| alloc_binary(&T, width, height, p) != DFX_SUCCESS || alloc_binary(&S, width, height, p) != DFX_SUCCESS
		|| alloc_labels(&L, width, height, p) != DFX_SUCCESS || alloc_labels(&L2, width, height, p) != DFX_SUCCESS) {
		printf("Error: cannot allocate memory\n"); return 1;
	}
	write_plane("dfx-binary-input.bmp", Y, width, height, p);

	/*
	 *  1. Histogram and global thresholds
	 */
	printf("\n1. Luminance histogram and global thresholds\n");
	histogram_plane(Y, width, height, p, 0.0f, 1.0f, 16, hist);
	print_histogram(hist, 16, 0.0f, 1.0f);
	otsu_threshold(Y, width, height, p, 0.0f, 1.0f, 256, &tau_otsu);
	ridler_calvard_threshold(Y, width, height, p, 0.0f, 1.0f, 256, &tau_rc);
	printf("  Otsu's threshold:           %.4f\n", tau_otsu);
	printf("  Ridler-Calvard's threshold: %.4f\n", tau_rc);
	threshold_plane(Y, A, width, height, p, tau_otsu);
	count_binary(A, width, height, p, &count);
	printf("  thresholding with Otsu's value: %d foreground pixels (%.1f%%)\n", count, 100.0 * count / ((double)width * height));
	write_binary("dfx-binary-threshold.bmp", A, width, height, p);
	if (synthetic) {
		check(tau_otsu > 0.2f && tau_otsu < 0.5f, "Otsu's threshold lies between the background and the objects");
		check(fabs(tau_otsu - tau_rc) < 0.1, "Otsu and Ridler-Calvard agree within 0.1");
	}

	/*
	 *  2. Clean-up: majority filter, opening, closing, hole filling
	 */
	printf("\n2. Clean-up\n");
	majority_binary(A, C, width, height, p, 1);
	count_binary(C, width, height, p, &count2);
	printf("  majority filter (3 x 3):     %d pixels (%+d)\n", count2, count2 - count); count = count2;
	structuring_element(se_big, 2, SE_DISK);
	open_binary(C, C, width, height, p, se_big, 2);
	count_binary(C, width, height, p, &count2);
	printf("  opening (disk of radius 2):  %d pixels (%+d)\n", count2, count2 - count); count = count2;
	close_binary(C, C, width, height, p, se_big, 2);
	count_binary(C, width, height, p, &count2);
	printf("  closing (disk of radius 2):  %d pixels (%+d)\n", count2, count2 - count); count = count2;
	write_binary("dfx-binary-clean.bmp", C, width, height, p);
	fill_holes_binary(C, T, width, height, p, 4);
	count_binary(T, width, height, p, &count2);
	printf("  hole filling:                %d pixels (%+d)\n", count2, count2 - count);
	write_binary("dfx-binary-filled.bmp", T, width, height, p);
	if (synthetic) {
		double ring_hole = M_PI * objects[1].r2 * objects[1].r2;
		check(fabs((count2 - count) - ring_hole) < 0.1 * ring_hole + 20, "hole filling adds the area of the ring's hole");
	}

	/*
	 *  3. Morphological operations
	 */
	printf("\n3. Morphology\n");
	structuring_element(se, 1, SE_CROSS);
	copy_binary(C, T, width, height, p);
	for (k = 0; k < 5; k++) dilate_binary(T, T, width, height, p, se, 1);
	structuring_element(se_big, 5, SE_CROSS);
	dilate_binary(C, S, width, height, p, se_big, 5);
	xor_binary(T, S, T, width, height, p); count_binary(T, width, height, p, &count);
	printf("  5 cascaded dilations by B4 vs one dilation by the cross of radius 5: %d differing pixels\n", count);
	check(count == 0, "cascaded B4 dilations equal the dilation by the larger cross");
	structuring_element(se, 1, SE_BOX);
	copy_binary(C, T, width, height, p);
	for (k = 0; k < 5; k++) erode_binary(T, T, width, height, p, se, 1);
	structuring_element(se_big, 5, SE_BOX);
	erode_binary(C, S, width, height, p, se_big, 5);
	xor_binary(T, S, T, width, height, p); count_binary(T, width, height, p, &count);
	printf("  5 cascaded erosions by B8 vs one erosion by the box of radius 5:     %d differing pixels\n", count);
	check(count == 0, "cascaded B8 erosions equal the erosion by the larger box");
	structuring_element(se_big, 5, SE_DISK);
	dilate_binary(C, T, width, height, p, se_big, 5);
	not_binary(C, S, width, height, p); erode_binary(S, S, width, height, p, se_big, 5); not_binary(S, S, width, height, p);
	xor_binary(T, S, T, width, height, p); count_binary(T, width, height, p, &count);
	printf("  dilation vs complement of the erosion of the complement (disk of radius 5): %d differing pixels\n", count);
	check(count == 0, "duality of dilation and erosion");
	boundary_binary(C, T, width, height, p, 8); count_binary(T, width, height, p, &count);
	boundary_binary(C, S, width, height, p, 4); count_binary(S, width, height, p, &count2);
	printf("  boundary pixels: %d (8-connected boundary), %d (4-connected boundary)\n", count, count2);
	check(count2 >= count, "the 4-connected boundary is at least as thick as the 8-connected one");
	write_binary("dfx-binary-boundary.bmp", T, width, height, p);

	/*
	 *  4. Connected components and region properties
	 */
	printf("\n4. Connected components (8-connectivity)\n");
	label_binary(C, L, width, height, p, 8, &n_labels);
	label_binary_floodfill(C, L2, width, height, p, 8, &n_labels2);
	for (count = 0, y = 0; y < height; y++) for (x = 0; x < width; x++)
		if (L[(p + y) * width_p + p + x] != L2[(p + y) * width_p + p + x]) count++;
	printf("  union-find: %d regions; floodfill: %d regions; %d pixels labeled differently\n", n_labels, n_labels2, count);
	check(n_labels == n_labels2 && count == 0, "both labeling algorithms give the same labels");
	if (synthetic) check(n_labels == N_OBJECTS, "one region per object");
	labels_to_image(L, R, G, B, width, height, p);
	write_rgb("dfx-binary-labels.bmp", R, G, B, width, height, p);

	if ((regions = (dfx_region*)malloc(((size_t)n_labels + 1) * sizeof(dfx_region))) == NULL) { printf("Error: cannot allocate memory\n"); return 1; }
	printf("\n  label    area   centroid (x,y)  theta(deg) ecc.    a      b   bound.  length  compact.%s\n", synthetic ? "  object" : "");
	for (k = 1; k <= n_labels; k++) {
		if ((err = region_properties(L, width, height, p, k, &regions[k - 1])) != DFX_SUCCESS) { printf("Error: region_properties(%d) failed (%d)\n", k, err); return 1; }
		if (k > 40) continue;			/* print at most 40 regions */
		printf("  %5d %7d  (%6.1f,%6.1f)  %6.1f    %5.3f %6.1f %6.1f %6d %8.1f  %6.3f", k, regions[k - 1].area,
			regions[k - 1].x_c, regions[k - 1].y_c, regions[k - 1].theta * 180 / M_PI, regions[k - 1].eccentricity,
			regions[k - 1].a, regions[k - 1].b, regions[k - 1].perimeter, regions[k - 1].perimeter_length, regions[k - 1].compactness);
		if (synthetic) {
			/* name the object whose center is nearest to the centroid: */
			int best = -1; double d_best = 1e30;
			for (n = 0; n < N_OBJECTS; n++) {
				double dd = (objects[n].x - regions[k - 1].x_c) * (objects[n].x - regions[k - 1].x_c) + (objects[n].y - regions[k - 1].y_c) * (objects[n].y - regions[k - 1].y_c);
				if (dd < d_best) { d_best = dd; best = n; }
			}
			printf("  %s", objects[best].name);
		}
		printf("\n");
	}
	if (n_labels > 40) printf("  ... (%d regions)\n", n_labels);

	/* boundary lengths of the first region by the three estimators: */
	if (n_labels > 0) {
		unsigned char* chain;
		int x0 = -1, y0 = -1, len = 0, len4 = 0;
		select_label(L, T, width, height, p, 1);
		/* the first pixel of the region in scanning order has an OFF pixel (or the image border) above it: */
		for (y = 0; y < height && x0 < 0; y++) for (x = 0; x < width; x++) if (T[(p + y) * width_p + p + x]) { x0 = x; y0 = y; break; }
		if ((err = trace_boundary_binary(T, width, height, p, x0, y0, 8, NULL, 0, &len)) != DFX_SUCCESS
			|| (err = trace_boundary_binary(T, width, height, p, x0, y0, 4, NULL, 0, &len4)) != DFX_SUCCESS) {
			printf("Error: trace_boundary_binary() failed (%d)\n", err); return 1;
		}
		if ((chain = (unsigned char*)malloc((size_t)len + 1)) != NULL) {
			if ((err = trace_boundary_binary(T, width, height, p, x0, y0, 8, chain, len, &len)) != DFX_SUCCESS) { printf("Error: trace_boundary_binary() failed (%d)\n", err); return 1; }
			printf("\n  region 1: chain code of %d steps starting at (%d,%d):", len, x0, y0);
			for (k = 0; k < len && k < 24; k++) printf(" %d", chain[k]);
			printf("%s\n", len > 24 ? " ..." : "");
			chain_code_length(chain, len, CHAIN_FREEMAN, &length);     printf("  boundary length: Freeman %.1f,", length);
			chain_code_length(chain, len, CHAIN_PYTHAGOREAN, &length); printf(" Pythagorean %.1f,", length);
			chain_code_length(chain, len, CHAIN_KIMURA, &length);      printf(" Kimura %.1f\n", length);
			free(chain);
		}
		printf("  the 4-connected boundary of region 1 has %d steps\n", len4);
	}

	if (synthetic) {
		double true_ecc = sqrt(1.0 - (objects[5].r2 * objects[5].r2) / (objects[5].r1 * objects[5].r1));
		k = nearest_region(regions, n_labels, objects[0].x, objects[0].y);		/* the big disk */
		check(k > 0 && fabs(regions[k - 1].area - M_PI * objects[0].r1 * objects[0].r1) < 0.03 * M_PI * objects[0].r1 * objects[0].r1, "big disk: area = pi r^2 within 3%");
		check(k > 0 && regions[k - 1].eccentricity < 0.1, "big disk: eccentricity < 0.1");
		check(k > 0 && fabs(regions[k - 1].perimeter_length - 2 * M_PI * objects[0].r1) < 0.04 * 2 * M_PI * objects[0].r1, "big disk: Kimura length = 2 pi r within 4%");
		check(k > 0 && regions[k - 1].compactness > 0.9 && regions[k - 1].compactness < 1.1, "big disk: compactness = 1 within 10%");
		check(k > 0 && fabs(regions[k - 1].a - objects[0].r1) < 1.0 && fabs(regions[k - 1].b - objects[0].r1) < 1.0, "big disk: semi-axes = r within 1 pixel");
		k = nearest_region(regions, n_labels, objects[5].x, objects[5].y);		/* the ellipse */
		check(k > 0 && fabs(regions[k - 1].theta * 180 / M_PI - objects[5].angle) < 2.0, "ellipse: orientation = 30 degrees within 2 degrees");
		check(k > 0 && fabs(regions[k - 1].eccentricity - true_ecc) < 0.02, "ellipse: eccentricity within 0.02 of the truth");
		check(k > 0 && fabs(regions[k - 1].a - objects[5].r1) < 1.5 && fabs(regions[k - 1].b - objects[5].r2) < 1.5, "ellipse: semi-axes within 1.5 pixels");
		k = nearest_region(regions, n_labels, objects[6].x, objects[6].y);		/* the bar */
		check(k > 0 && regions[k - 1].eccentricity > 0.99 && fabs(regions[k - 1].theta * 180 / M_PI - objects[6].angle) < 2.0, "thin bar: eccentricity > 0.99, orientation -25 degrees");
		k = nearest_region(regions, n_labels, objects[1].x, objects[1].y);		/* the ring */
		check(k > 0 && fabs(regions[k - 1].area - M_PI * (objects[1].r1 * objects[1].r1 - objects[1].r2 * objects[1].r2)) < 0.05 * M_PI * objects[1].r1 * objects[1].r1, "ring: area = pi (r1^2 - r2^2) within 5% of the disk area");
	}

	/*
	 *  5. Chamfer distances
	 */
	printf("\n5. Chamfer distance transform (a = 1, b = 1.35)\n");
	if (chamfer_binary(C, D, width, height, p, 1.0f, 1.35f) != DFX_SUCCESS) {
		printf("  the image has no foreground pixels: no distances to compute\n");
		zero_plane(D, width, height, p);
	}
	for (d_max = 0, y = 0; y < height; y++) for (x = 0; x < width; x++) if (D[(p + y) * width_p + p + x] > d_max) d_max = D[(p + y) * width_p + p + x];
	printf("  distance to the nearest object: maximum %.1f pixels\n", d_max);
	if (d_max > 0) scale_plane(D, D, width, height, p, 1.0f / d_max);
	write_plane("dfx-binary-chamfer.bmp", D, width, height, p);
	/* inside the objects: distance to the background, and its maximum in every region (the center of the largest inscribed ball) */
	not_binary(C, T, width, height, p);
	if (chamfer_binary(T, D, width, height, p, 1.0f, 1.35f) != DFX_SUCCESS) {
		printf("  the image has no background pixels: no distances to compute\n");
		zero_plane(D, width, height, p);
	}
	printf("  inside the objects, distance to the background:\n");
	for (k = 1; k <= n_labels && k <= 40; k++) {
		int xm = 0, ym = 0;
		d_max = -1;
		for (y = regions[k - 1].y_min; y <= regions[k - 1].y_max; y++) for (x = regions[k - 1].x_min; x <= regions[k - 1].x_max; x++)
			if (L[(p + y) * width_p + p + x] == k && D[(p + y) * width_p + p + x] > d_max) { d_max = D[(p + y) * width_p + p + x]; xm = x; ym = y; }
		printf("  region %d: largest distance %.2f at (%d,%d); centroid (%.1f,%.1f) is %s the region\n", k, d_max, xm, ym,
			regions[k - 1].x_c, regions[k - 1].y_c,
			L[(p + (int)(regions[k - 1].y_c + 0.5)) * width_p + p + (int)(regions[k - 1].x_c + 0.5)] == k ? "inside" : "outside");
		if (synthetic && k == nearest_region(regions, n_labels, objects[0].x, objects[0].y))
			check(fabs(d_max - objects[0].r1) < 0.08 * objects[0].r1, "big disk: the largest inside distance equals the radius within 8%");
		if (synthetic && k == nearest_region(regions, n_labels, objects[4].x - 40, objects[4].y + 20))
			check(L[(p + (int)(regions[k - 1].y_c + 0.5)) * width_p + p + (int)(regions[k - 1].x_c + 0.5)] != k, "crescent: the centroid lies outside the region");
	}

	/*
	 *  6. Skeletons by thinning
	 */
	printf("\n6. Thinning\n");
	/* The padding of the input tells the morphological operations what lies outside the image. The clean
	 * image is replicate-padded (every operation leaves its output that way), which makes an object touching
	 * the image boundary continue beyond it, so that thinning would never remove its pixels along that
	 * boundary. For skeletons the outside is better treated as background: */
	copy_binary(C, T, width, height, p);
	pad_binary(T, width, height, p, PAD_ZERO);
	thin_binary(T, S, width, height, p, 0, &iter);
	count_binary(S, width, height, p, &count);
	printf("  thinning converged after %d passes: %d skeleton pixels\n", iter, count);
	label_binary(S, L2, width, height, p, 8, &n);
	printf("  the skeleton has %d connected components (the image has %d regions)\n", n, n_labels);
	check(n == n_labels, "thinning preserves the number of regions");
	/* the skeleton is thin: the eight structuring elements leave thick spots only where several branches meet */
	structuring_element(se, 1, SE_BOX);
	erode_binary(S, T, width, height, p, se, 1); count_binary(T, width, height, p, &count2);
	printf("  pixels of the skeleton with all 8 neighbors ON: %d\n", count2);
	if (synthetic) check(count2 == 0, "the skeleton contains no 3 x 3 block");
	/* show the skeleton in white over the objects in gray: */
	binary_to_plane(C, D, width, height, p); scale_plane(D, D, width, height, p, 0.15f);
	binary_to_plane(S, R, width, height, p); add_planes(D, R, D, width, height, p); clip_plane(D, D, width, height, p, 0.0f, 1.0f);
	write_plane("dfx-binary-skeleton.bmp", D, width, height, p);

	/*
	 *  7. Double thresholding, adaptive thresholding, region growing
	 */
	printf("\n7. Double thresholding, Niblack's adaptive thresholding, region growing\n");
	double_threshold_plane(Y, T, width, height, p, 0.2f, tau_otsu, 8);
	majority_binary(T, T, width, height, p, 1);
	count_binary(T, width, height, p, &count);
	printf("  double threshold (0.2, %.3f) + majority filter: %d pixels\n", tau_otsu, count);
	write_binary("dfx-binary-double.bmp", T, width, height, p);
	if (synthetic) {
		/* the dim disk: the single threshold cuts it where its luminance falls below tau; the double threshold recovers it */
		int area_single = 0, area_double = 0;
		for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
			double dx = x - objects[3].x, dy = y - objects[3].y;
			if (dx * dx + dy * dy > objects[3].r1 * objects[3].r1) continue;
			area_single += C[(p + y) * width_p + p + x];
			area_double += T[(p + y) * width_p + p + x];
		}
		printf("  the dim disk covers %.0f pixels: %d found by the single threshold, %d by the double threshold\n",
			M_PI * objects[3].r1 * objects[3].r1, area_single, area_double);
		check(area_single < 0.85 * M_PI * objects[3].r1 * objects[3].r1, "dim disk: the single threshold misses a part of it");
		check(area_double > 0.9 * M_PI * objects[3].r1 * objects[3].r1, "dim disk: the double threshold recovers it");
	}
	niblack_threshold_plane(Y, T, width, height, p, 15, 0.2f);
	count_binary(T, width, height, p, &count);
	printf("  Niblack's threshold (31 x 31 neighborhood, k = 0.2): %d pixels\n", count);
	printf("  (far from the objects the local threshold sits inside the background noise, so about half of those pixels turn ON)\n");
	write_binary("dfx-binary-niblack.bmp", T, width, height, p);
	if (synthetic) {
		/* seed: the brightest pixel within 3 pixels of the center of the big disk (the center itself may be a pepper speck) */
		int xc = (int)objects[0].x, yc = (int)objects[0].y, xs = xc, ys = yc;
		for (y = yc - 3; y <= yc + 3; y++) for (x = xc - 3; x <= xc + 3; x++)
			if (Y[(p + y) * width_p + p + x] > Y[(p + ys) * width_p + p + xs]) { xs = x; ys = y; }
		grow_region_plane(Y, T, width, height, p, xs, ys, 0.15f, 4);
		count_binary(T, width, height, p, &count);
		printf("  region grown from (%d,%d) near the center of the big disk (tolerance 0.15): %d pixels (pi r^2 = %.0f)\n", xs, ys, count, M_PI * objects[0].r1 * objects[0].r1);
		check(fabs(count - M_PI * objects[0].r1 * objects[0].r1) < 0.05 * M_PI * objects[0].r1 * objects[0].r1, "region growing recovers the big disk within 5%");
	} else {
		grow_region_plane(Y, T, width, height, p, width / 2, height / 2, 0.1f, 4);
		count_binary(T, width, height, p, &count);
		printf("  region grown from the image center (tolerance 0.1): %d pixels\n", count);
	}
	write_binary("dfx-binary-grown.bmp", T, width, height, p);

	/* summary */
	if (synthetic) printf("\n%s\n", checks_failed ? "SOME CHECKS FAILED" : "All checks passed.");

	/* free memory and exit: */
	free(regions);
	free_plane(Y); free_plane(D);
	free_image(R, G, B);
	free_binary(A); free_binary(C); free_binary(T); free_binary(S);
	free_labels(L); free_labels(L2);
	return checks_failed ? 1 : 0;
}
