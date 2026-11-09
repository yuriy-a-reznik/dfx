/*!
 *  \file   dfx_regions.c
 *  \brief  Regions in binary images: connectivity, labeling, region properties, boundaries, distances.
 *
 *  A region is a connected set of ON pixels of a binary image. Two ON pixels are 4-connected if one
 *  is the left, right, upper, or lower neighbor of the other, and 8-connected if they are 4-connected
 *  or diagonal neighbors; every function in this module taking the parameter "connectivity" accepts
 *  the values 4 and 8. The connected components of an image are numbered 1, 2, 3, ... in the order
 *  in which they are first met when scanning the image row by row, and stored in a "label plane":
 *  an array of ints with the same geometry as the binary image (padding included, the border being
 *  filled with zeros), holding 0 for background pixels and the component number for the others.
 *
 *  The algorithms follow S. Birchfield, Image Processing and Analysis (Cengage, 2017), chapter 4:
 *  floodfill with an explicit frontier stack (Algorithm 4.4), connected components by union-find
 *  (Algorithm 4.8), region moments and the equivalent ellipse, boundary tracing and Freeman chain
 *  codes (Fig. 4.27), and the two-pass chamfer distance transform (Algorithm 4.10).
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

/*
 *  Unit steps in the 8 directions, numbered as Freeman chain codes (Birchfield, Fig. 4.27):
 *  0 = east, 1 = north-east, 2 = north, 3 = north-west, 4 = west, 5 = south-west, 6 = south,
 *  7 = south-east, where "north" is up, i.e. towards decreasing y. Increasing codes turn
 *  counterclockwise as seen on the screen. The 4-neighbors are the even directions.
 */
static const int DX[8] = { 1,  1,  0, -1, -1, -1,  0,  1 };
static const int DY[8] = { 0, -1, -1, -1,  0,  1,  1,  1 };

/*!
 *  \brief Direction code of a unit step (dx,dy), or -1 if it is not one.
 */
static int direction_of(int dx, int dy)
{
	int k;
	for (k = 0; k < 8; k++) if (DX[k] == dx && DY[k] == dy) return k;
	return -1;
}

/*!
 *  \brief Fill the list of neighbor directions for the given connectivity; returns their number or 0 if invalid.
 */
static int neighbor_directions(int connectivity, int* dirs)
{
	int k, n = 0;
	if (connectivity != 4 && connectivity != 8) return 0;
	for (k = 0; k < 8; k++) if (connectivity == 8 || (k & 1) == 0) dirs[n++] = k;
	return n;
}

/**************
 *
 *  Label planes:
 *
 *   alloc_labels()
 *   free_labels()
 *   select_label()
 *   labels_to_image()
 */

/*!
 *  \brief Allocate a label plane for a width x height image with padding p.
 */
int alloc_labels(int** pL, int width, int height, int p)
{
	size_t size;
	if (pL == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	size = (size_t)(2 * p + height) * (size_t)(2 * p + width);
	if ((*pL = (int*)malloc((size > 0 ? size : 1) * sizeof(int))) == NULL) return DFX_NOMEM;
	return DFX_SUCCESS;
}

/*!
 *  \brief Free a label plane.
 */
int free_labels(int* L)
{
	if (L != NULL) free(L);
	return DFX_SUCCESS;
}

/*!
 *  \brief Extract one region from a label plane: A = (L == label). The output receives PAD_REPLICATE padding.
 */
int select_label(const int* L, unsigned char* A, int width, int height, int p, int label)
{
	int width_p = width + 2 * p;
	int x, y;
	if (L == NULL || A == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	for (y = 0; y < height; y++) for (x = 0; x < width; x++)
		A[(p + y) * width_p + p + x] = (unsigned char)(L[(p + y) * width_p + p + x] == label ? 1 : 0);
	return pad_binary(A, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Pseudocolor for a label: hues stepped by the golden ratio, so that consecutive labels differ as much as possible.
 */
static void label_color(int label, float* r, float* g, float* b)
{
	double h = fmod(0.618033988749895 * (double)label, 1.0) * 6.0;		/* hue sector in [0,6) */
	int i = (int)h;
	double f = h - i;
	double q = 1.0 - f, t = f;
	switch (i) {
	case 0:  *r = 1.0f;     *g = (float)t; *b = 0.0f;     break;
	case 1:  *r = (float)q; *g = 1.0f;     *b = 0.0f;     break;
	case 2:  *r = 0.0f;     *g = 1.0f;     *b = (float)t; break;
	case 3:  *r = 0.0f;     *g = (float)q; *b = 1.0f;     break;
	case 4:  *r = (float)t; *g = 0.0f;     *b = 1.0f;     break;
	default: *r = 1.0f;     *g = 0.0f;     *b = (float)q; break;
	}
}

/*!
 *  \brief Render a label plane as a pseudocolored linear RGB image: black background, a distinct hue per label.
 *
 *  The output planes receive PAD_REPLICATE padding.
 */
int labels_to_image(const int* L, float* R, float* G, float* B, int width, int height, int p)
{
	int width_p = width + 2 * p;
	int x, y, i;
	if (L == NULL || R == NULL || G == NULL || B == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		i = (p + y) * width_p + p + x;
		if (L[i] == 0) R[i] = G[i] = B[i] = 0.0f;
		else label_color(L[i], R + i, G + i, B + i);
	}
	return pad_image(R, G, B, width, height, p, PAD_REPLICATE);
}

/**************
 *
 *  Floodfill and connected components:
 *
 *   floodfill_labels()
 *   label_binary()
 *   label_binary_floodfill()
 *   reconstruct_binary()
 *   fill_holes_binary()
 */

/*!
 *  \brief Floodfill a label plane: replace the value of the connected region containing (x0,y0) by new_label.
 *
 *  Implements Birchfield's Algorithm 4.4 with an explicit frontier stack (the recursive version
 *  would overflow the C stack on large regions). Pixels are tested and relabeled when they are
 *  pushed, so that every pixel enters the stack at most once and the stack needs no more than
 *  width * height entries.
 *
 *  \param[in,out]  L            - label plane (or any plane of ints with the standard geometry)
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter
 *  \param[in]      x0, y0       - seed pixel
 *  \param[in]      new_label    - the value to fill with
 *  \param[in]      connectivity - 4 or 8
 *  \param[out]     p_count      - if not NULL, receives the number of pixels filled
 */
int floodfill_labels(int* L, int width, int height, int p, int x0, int y0, int new_label, int connectivity, int* p_count)
{
	int width_p = width + 2 * p;
	int dirs[8], n_dirs, k, old, count = 0;
	int* stack = NULL;
	size_t sp = 0;

	/* check parameters: */
	if (L == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (x0 < 0 || x0 >= width || y0 < 0 || y0 >= height) return DFX_INVARG;
	if ((n_dirs = neighbor_directions(connectivity, dirs)) == 0) return DFX_INVARG;

	old = L[(p + y0) * width_p + p + x0];
	if (old != new_label) {
		if ((stack = (int*)malloc((size_t)width * height * sizeof(int))) == NULL) return DFX_NOMEM;
		/* seed: */
		L[(p + y0) * width_p + p + x0] = new_label;
		stack[sp++] = y0 * width + x0;
		count = 1;
		/* grow: */
		while (sp > 0) {
			int q = stack[--sp], qx = q % width, qy = q / width;
			for (k = 0; k < n_dirs; k++) {
				int nx = qx + DX[dirs[k]], ny = qy + DY[dirs[k]];
				if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
				if (L[(p + ny) * width_p + p + nx] != old) continue;
				L[(p + ny) * width_p + p + nx] = new_label;
				stack[sp++] = ny * width + nx;
				count++;
			}
		}
		free(stack);
	}
	if (p_count != NULL) *p_count = count;
	return DFX_SUCCESS;
}

/*!
 *  \brief Union-find: find the representative (smallest equivalent label) of i, compressing the path on the way.
 */
static int uf_find(int* parent, int i)
{
	while (parent[i] != i) {
		parent[i] = parent[parent[i]];
		i = parent[i];
	}
	return i;
}

/*!
 *  \brief Union-find: declare labels i and j equivalent; the smaller representative becomes the parent.
 */
static void uf_union(int* parent, int i, int j)
{
	i = uf_find(parent, i);
	j = uf_find(parent, j);
	if (i < j) parent[j] = i;
	else if (j < i) parent[i] = j;
}

/*!
 *  \brief Label the connected components of a binary image by the two-pass union-find algorithm.
 *
 *  Implements Birchfield's Algorithm 4.8. In the first pass, every ON pixel inherits the smallest
 *  label among its already visited neighbors (left and upper for 4-connectivity, plus the two upper
 *  diagonal neighbors for 8-connectivity), or receives a new provisional label if there is none, and
 *  the labels of these neighbors are recorded as equivalent in a union-find structure. In the second
 *  pass, every provisional label is replaced by the representative of its equivalence class, and the
 *  classes are renumbered 1, 2, 3, ... in the order of their first appearance. The border of L is zero.
 *
 *  \param[in]      A            - input binary image
 *  \param[out]     L            - output label plane: 0 for background, the component number otherwise
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter of both
 *  \param[in]      connectivity - 4 or 8
 *  \param[out]     p_n_labels   - if not NULL, receives the number of components found
 */
int label_binary(const unsigned char* A, int* L, int width, int height, int p, int connectivity, int* p_n_labels)
{
	int width_p = width + 2 * p;
	int x, y, k, count, n_prov = 0, n_final = 0;
	int* parent = NULL, * final_label = NULL;

	/* check parameters: */
	if (A == NULL || L == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (connectivity != 4 && connectivity != 8) return DFX_INVARG;

	/* clear the label plane (padding included) and allocate the equivalence table: */
	memset((void*)L, 0, (size_t)(2 * p + height) * (size_t)(2 * p + width) * sizeof(int));
	count_binary(A, width, height, p, &count);
	if ((parent = (int*)malloc(((size_t)count + 1) * sizeof(int))) == NULL) return DFX_NOMEM;

	/* first pass: provisional labels and equivalences */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		int nb[4], nn = 0, label;
		if (!A[(p + y) * width_p + p + x]) continue;
		/* labels of the already visited ON neighbors: */
		if (x > 0 && L[(p + y) * width_p + p + x - 1]) nb[nn++] = L[(p + y) * width_p + p + x - 1];
		if (y > 0) {
			if (L[(p + y - 1) * width_p + p + x]) nb[nn++] = L[(p + y - 1) * width_p + p + x];
			if (connectivity == 8) {
				if (x > 0 && L[(p + y - 1) * width_p + p + x - 1]) nb[nn++] = L[(p + y - 1) * width_p + p + x - 1];
				if (x < width - 1 && L[(p + y - 1) * width_p + p + x + 1]) nb[nn++] = L[(p + y - 1) * width_p + p + x + 1];
			}
		}
		if (nn == 0) {
			/* a new region starts here: */
			label = ++n_prov;
			parent[label] = label;
		} else {
			/* take the smallest neighbor label and record the equivalences: */
			label = nb[0];
			for (k = 1; k < nn; k++) if (nb[k] < label) label = nb[k];
			for (k = 0; k < nn; k++) uf_union(parent, label, nb[k]);
		}
		L[(p + y) * width_p + p + x] = label;
	}

	/* second pass: replace provisional labels by final ones, numbered in order of first appearance */
	if ((final_label = (int*)calloc((size_t)n_prov + 1, sizeof(int))) == NULL) { free(parent); return DFX_NOMEM; }
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		int i = (p + y) * width_p + p + x, root;
		if (!L[i]) continue;
		root = uf_find(parent, L[i]);
		if (final_label[root] == 0) final_label[root] = ++n_final;
		L[i] = final_label[root];
	}
	free(final_label);
	free(parent);
	if (p_n_labels != NULL) *p_n_labels = n_final;
	return DFX_SUCCESS;
}

/*!
 *  \brief Label the connected components of a binary image by repeated floodfills.
 *
 *  The image is scanned row by row; every ON pixel not yet labeled starts a new component, which is
 *  then filled with the next label. This simpler algorithm produces exactly the same label plane as
 *  label_binary(). Same parameters as label_binary().
 */
int label_binary_floodfill(const unsigned char* A, int* L, int width, int height, int p, int connectivity, int* p_n_labels)
{
	int width_p = width + 2 * p;
	int x, y, k, dirs[8], n_dirs, n_labels = 0;
	int* stack = NULL;
	size_t sp;

	/* check parameters: */
	if (A == NULL || L == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if ((n_dirs = neighbor_directions(connectivity, dirs)) == 0) return DFX_INVARG;

	/* clear the label plane (padding included): */
	memset((void*)L, 0, (size_t)(2 * p + height) * (size_t)(2 * p + width) * sizeof(int));
	if ((stack = (int*)malloc(((size_t)width * height + 1) * sizeof(int))) == NULL) return DFX_NOMEM;

	/* scan the image: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		if (!A[(p + y) * width_p + p + x] || L[(p + y) * width_p + p + x]) continue;
		/* a new region: fill it */
		n_labels++;
		L[(p + y) * width_p + p + x] = n_labels;
		stack[0] = y * width + x;
		sp = 1;
		while (sp > 0) {
			int q = stack[--sp], qx = q % width, qy = q / width;
			for (k = 0; k < n_dirs; k++) {
				int nx = qx + DX[dirs[k]], ny = qy + DY[dirs[k]], i;
				if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
				i = (p + ny) * width_p + p + nx;
				if (!A[i] || L[i]) continue;
				L[i] = n_labels;
				stack[sp++] = ny * width + nx;
			}
		}
	}
	free(stack);
	if (p_n_labels != NULL) *p_n_labels = n_labels;
	return DFX_SUCCESS;
}

/*!
 *  \brief Morphological reconstruction: keep the connected components of mask that contain at least one ON pixel of seeds.
 *
 *  Every ON pixel of seeds that is also ON in mask starts a floodfill of mask; the union of the
 *  filled components is the result (in morphological terms, the geodesic dilation of seeds within
 *  mask, iterated to stability). The output receives PAD_REPLICATE padding and may be the same
 *  as mask or seeds.
 *
 *  \param[in]      mask         - binary image whose components are selected
 *  \param[in]      seeds        - binary image of markers
 *  \param[out]     A_out        - output binary image
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter of all
 *  \param[in]      connectivity - 4 or 8
 */
int reconstruct_binary(const unsigned char* mask, const unsigned char* seeds, unsigned char* A_out, int width, int height, int p, int connectivity)
{
	int width_p = width + 2 * p;
	int x, y, k, dirs[8], n_dirs, err;
	int* stack = NULL;
	unsigned char* out = NULL;
	size_t sp = 0;

	/* check parameters: */
	if (mask == NULL || seeds == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if ((n_dirs = neighbor_directions(connectivity, dirs)) == 0) return DFX_INVARG;

	/* work in a separate image, so that the output may alias an input: */
	if ((err = alloc_binary(&out, width, height, p)) != DFX_SUCCESS) return err;
	fill_binary(out, width, height, p, 0);
	if ((stack = (int*)malloc(((size_t)width * height + 1) * sizeof(int))) == NULL) { free_binary(out); return DFX_NOMEM; }

	/* push the seeds that lie within the mask: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		int i = (p + y) * width_p + p + x;
		if (seeds[i] && mask[i]) { out[i] = 1; stack[sp++] = y * width + x; }
	}
	/* grow them within the mask: */
	while (sp > 0) {
		int q = stack[--sp], qx = q % width, qy = q / width;
		for (k = 0; k < n_dirs; k++) {
			int nx = qx + DX[dirs[k]], ny = qy + DY[dirs[k]], i;
			if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
			i = (p + ny) * width_p + p + nx;
			if (!mask[i] || out[i]) continue;
			out[i] = 1;
			stack[sp++] = ny * width + nx;
		}
	}
	free(stack);
	copy_binary(out, A_out, width, height, p);
	free_binary(out);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/*!
 *  \brief Fill the holes of the foreground: the background components that do not touch the image boundary.
 *
 *  The background is flooded from the image boundary; whatever background remains unreached is a hole
 *  and is turned ON. The connectivity refers to the background (use 4 for 8-connected foreground
 *  regions and vice versa). The output receives PAD_REPLICATE padding and may be the same as the input.
 */
int fill_holes_binary(const unsigned char* A, unsigned char* A_out, int width, int height, int p, int connectivity)
{
	int width_p = width + 2 * p;
	int x, y, k, dirs[8], n_dirs, err;
	int* stack = NULL;
	unsigned char* bg = NULL;
	size_t sp = 0;

	/* check parameters: */
	if (A == NULL || A_out == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if ((n_dirs = neighbor_directions(connectivity, dirs)) == 0) return DFX_INVARG;

	/* bg marks the background reachable from the image boundary: */
	if ((err = alloc_binary(&bg, width, height, p)) != DFX_SUCCESS) return err;
	fill_binary(bg, width, height, p, 0);
	if ((stack = (int*)malloc(((size_t)width * height + 1) * sizeof(int))) == NULL) { free_binary(bg); return DFX_NOMEM; }

	/* seeds: OFF pixels on the image boundary */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		int i = (p + y) * width_p + p + x;
		if (x > 0 && x < width - 1 && y > 0 && y < height - 1) continue;
		if (!A[i]) { bg[i] = 1; stack[sp++] = y * width + x; }
	}
	/* flood the background: */
	while (sp > 0) {
		int q = stack[--sp], qx = q % width, qy = q / width;
		for (k = 0; k < n_dirs; k++) {
			int nx = qx + DX[dirs[k]], ny = qy + DY[dirs[k]], i;
			if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
			i = (p + ny) * width_p + p + nx;
			if (A[i] || bg[i]) continue;
			bg[i] = 1;
			stack[sp++] = ny * width + nx;
		}
	}
	free(stack);

	/* everything that is not reachable background is foreground: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		int i = (p + y) * width_p + p + x;
		A_out[i] = (unsigned char)(bg[i] ? 0 : 1);
	}
	free_binary(bg);
	return pad_binary(A_out, width, height, p, PAD_REPLICATE);
}

/**************
 *
 *  Boundaries:
 *
 *   trace_boundary_binary()
 *   chain_code_length()
 */

/*!
 *  \brief ON test with bounds checking: pixels outside of the image are OFF, whatever the padding holds.
 */
static int pixel_on(const unsigned char* A, int width, int height, int p, int x, int y)
{
	if (x < 0 || x >= width || y < 0 || y >= height) return 0;
	return A[(p + y) * (width + 2 * p) + p + x] != 0;
}

/*!
 *  \brief Trace the boundary of a region and encode it as a Freeman chain code.
 *
 *  Starting from an ON pixel (x0,y0) whose upper neighbor is OFF (for example the first pixel of
 *  a region in row-by-row scanning order), the boundary is followed clockwise (as seen on the
 *  screen) until it closes, and the direction of every step is recorded as a chain code 0..7
 *  (see the DX/DY tables above; even codes are horizontal and vertical steps, odd codes diagonal
 *  ones). The number of steps is the number of boundary pixels visited; a thin part of a region is
 *  traversed twice, once in each direction.
 *
 *  With connectivity = 8, the 8-connected boundary is traced by Moore's neighbor-tracing algorithm:
 *  at each pixel, the 8 neighbors are examined clockwise starting from the background pixel we
 *  came along, and the first ON pixel found is the next boundary pixel. The trace stops when the
 *  starting pixel is re-entered in the same way as initially (Jacob's stopping criterion), which
 *  makes the trace correct for regions one pixel wide. With connectivity = 4, the 4-connected
 *  boundary is traced by the wall-following algorithm: starting eastwards with the background on
 *  the left, at each step "if the pixel on your left is ON, turn left and move forward; otherwise,
 *  if the pixel in front of you is OFF, turn right; otherwise move forward", until the starting
 *  pixel is reached heading east again.
 *
 *  Pixels outside of the image are treated as OFF regardless of the padding. The trace follows the
 *  boundary between the region and the background component that contains the OFF pixel above the
 *  start: starting from the first pixel of a region in scanning order gives its outer boundary (holes
 *  are then ignored); starting from a pixel just below a hole traces the boundary of that hole.
 *
 *  \param[in]      A            - binary image
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter
 *  \param[in]      x0, y0       - starting pixel: must be ON, with an OFF pixel (or the image border) above it
 *  \param[in]      connectivity - 4 or 8
 *  \param[out]     chain        - array to receive the chain codes, or NULL to only count them
 *  \param[in]      max_length   - capacity of the chain array (ignored if chain is NULL)
 *  \param[out]     p_length     - receives the number of chain codes (0 for an isolated pixel)
 *  \return DFX_INVARG if (x0,y0) is not a valid starting pixel, DFX_NOMEM if the chain does not fit into max_length codes,
 *          DFX_INVRANGE if the trace does not close (which cannot happen for a valid start).
 */
int trace_boundary_binary(const unsigned char* A, int width, int height, int p, int x0, int y0, int connectivity, unsigned char* chain, int max_length, int* p_length)
{
	size_t iter, max_iter;
	int cx = x0, cy = y0, length = 0, d, k;

	/* check parameters: */
	if (A == NULL || p_length == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	if (connectivity != 4 && connectivity != 8) return DFX_INVARG;
	if (chain != NULL && max_length < 0) return DFX_INVARG;
	if (!pixel_on(A, width, height, p, x0, y0) || pixel_on(A, width, height, p, x0, y0 - 1)) return DFX_INVARG;

	max_iter = 8 * (size_t)width * (size_t)height + 8;	/* every boundary pixel is visited at most a few times: a safety net only */

	if (connectivity == 8) {
		/* Moore's neighbor tracing: b is the background pixel we arrived along ("backtrack") */
		int bx = x0, by = y0 - 1, p1x = 0, p1y = 0, nx = 0, ny = 0, lx, ly, found;
		for (iter = 0; ; iter++) {
			if (iter >= max_iter) return DFX_INVRANGE;
			d = direction_of(bx - cx, by - cy);
			lx = bx; ly = by;
			found = 0;
			for (k = 1; k <= 7; k++) {								/* examine the neighbors clockwise, starting after b */
				int dd = (d - k + 8) % 8;
				nx = cx + DX[dd]; ny = cy + DY[dd];
				if (pixel_on(A, width, height, p, nx, ny)) { d = dd; found = 1; break; }
				lx = nx; ly = ny;										/* the last background pixel examined */
			}
			if (!found) break;											/* an isolated pixel */
			if (length > 0 && cx == x0 && cy == y0 && nx == p1x && ny == p1y) break;	/* Jacob's stopping criterion */
			if (chain != NULL) {
				if (length >= max_length) return DFX_NOMEM;
				chain[length] = (unsigned char)d;
			}
			if (length == 0) { p1x = nx; p1y = ny; }
			length++;
			bx = lx; by = ly;
			cx = nx; cy = ny;
		}
	} else {
		/* wall following: heading d, the wall (background) on the left */
		d = 0;
		for (iter = 0; ; iter++) {
			int left = (d + 2) % 8, right = (d + 6) % 8, move = 0;
			if (iter >= 4 * max_iter) return DFX_INVRANGE;
			if (pixel_on(A, width, height, p, cx + DX[left], cy + DY[left])) { d = left; move = 1; }
			else if (!pixel_on(A, width, height, p, cx + DX[d], cy + DY[d])) d = right;
			else move = 1;
			if (move) {
				if (chain != NULL) {
					if (length >= max_length) return DFX_NOMEM;
					chain[length] = (unsigned char)d;
				}
				length++;
				cx += DX[d]; cy += DY[d];
			}
			if (cx == x0 && cy == y0 && d == 0) break;
		}
	}
	*p_length = length;
	return DFX_SUCCESS;
}

/*!
 *  \brief Estimate the length of a boundary from its chain code.
 *
 *  With N_o horizontal/vertical (isothetic) steps and N_d diagonal steps:
 *
 *     CHAIN_FREEMAN:      L = N_o + sqrt(2) N_d                          (the sum of the step lengths),
 *     CHAIN_PYTHAGOREAN:  L = sqrt(N_d^2 + (N_o + N_d)^2)                (a straight line with the same displacements),
 *     CHAIN_KIMURA:       L = sqrt(N_d^2 + (N_d + N_o / 2)^2) + N_o / 2  (Kimura's correction, close to the true length of smooth curves).
 *
 *  \param[in]      chain    - chain codes 0..7
 *  \param[in]      length   - number of codes
 *  \param[in]      method   - CHAIN_FREEMAN, CHAIN_PYTHAGOREAN, or CHAIN_KIMURA
 *  \param[out]     p_length - receives the estimated length
 */
int chain_code_length(const unsigned char* chain, int length, int method, double* p_length)
{
	int k, n_o = 0, n_d = 0;
	if ((chain == NULL && length > 0) || length < 0 || p_length == NULL) return DFX_INVARG;
	if (method < 0 || method >= N_CHAIN_METHODS) return DFX_INVARG;
	for (k = 0; k < length; k++) {
		if (chain[k] > 7) return DFX_INVARG;
		if (chain[k] & 1) n_d++; else n_o++;
	}
	switch (method) {
	case CHAIN_FREEMAN:     *p_length = n_o + sqrt(2.0) * n_d; break;
	case CHAIN_PYTHAGOREAN: *p_length = sqrt((double)n_d * n_d + (double)(n_o + n_d) * (n_o + n_d)); break;
	default:                *p_length = sqrt((double)n_d * n_d + (n_d + 0.5 * n_o) * (n_d + 0.5 * n_o)) + 0.5 * n_o; break;
	}
	return DFX_SUCCESS;
}

/**************
 *
 *  Region properties:
 *
 *   region_properties()
 *   region_properties_binary()
 */

/*!
 *  \brief Compute the properties of the set of ON pixels of a binary image (shared implementation).
 */
static int binary_properties(const unsigned char* A, int width, int height, int p, dfx_region* r)
{
	int width_p = width + 2 * p;
	int x, y, area = 0, perimeter = 0, x_first = 0, y_first = 0, length, err;
	double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0, S, T;
	unsigned char* chain;

	memset((void*)r, 0, sizeof(*r));

	/* first pass: area, bounding box, moments up to order 2, boundary pixel count */
	r->x_min = width; r->y_min = height; r->x_max = r->y_max = -1;
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		if (!A[(p + y) * width_p + p + x]) continue;
		if (area == 0) { x_first = x; y_first = y; }
		area++;
		if (x < r->x_min) r->x_min = x;
		if (x > r->x_max) r->x_max = x;
		if (y < r->y_min) r->y_min = y;
		if (y > r->y_max) r->y_max = y;
		sx += x; sy += y;
		sxx += (double)x * x; sxy += (double)x * y; syy += (double)y * y;
		if (!pixel_on(A, width, height, p, x - 1, y) || !pixel_on(A, width, height, p, x + 1, y)
			|| !pixel_on(A, width, height, p, x, y - 1) || !pixel_on(A, width, height, p, x, y + 1)) perimeter++;
	}
	if (area == 0) { r->x_min = r->y_min = 0; return DFX_INVRANGE; }
	r->area = area;
	r->perimeter = perimeter;

	/* centroid and central moments: */
	r->x_c = sx / area;
	r->y_c = sy / area;
	r->mu20 = sxx - r->x_c * sx;
	r->mu02 = syy - r->y_c * sy;
	r->mu11 = sxy - r->x_c * sy;

	/* equivalent ellipse: orientation, eccentricity, semi-axes */
	S = sqrt((r->mu20 - r->mu02) * (r->mu20 - r->mu02) + 4.0 * r->mu11 * r->mu11);
	T = r->mu20 + r->mu02;
	r->theta = 0.5 * atan2(2.0 * r->mu11, r->mu20 - r->mu02);
	r->eccentricity = (T + S > 0) ? sqrt(2.0 * S / (T + S)) : 0.0;
	r->a = sqrt(2.0 * (T + S) / area);
	r->b = (T - S > 0) ? sqrt(2.0 * (T - S) / area) : 0.0;

	/* outer boundary length (Kimura's estimate from the chain code) and compactness: */
	if ((err = trace_boundary_binary(A, width, height, p, x_first, y_first, 8, NULL, 0, &length)) != DFX_SUCCESS) return err;
	if (length > 0) {
		if ((chain = (unsigned char*)malloc((size_t)length)) == NULL) return DFX_NOMEM;
		err = trace_boundary_binary(A, width, height, p, x_first, y_first, 8, chain, length, &length);
		if (err == DFX_SUCCESS) err = chain_code_length(chain, length, CHAIN_KIMURA, &r->perimeter_length);
		free(chain);
		if (err != DFX_SUCCESS) return err;
		r->compactness = 4.0 * M_PI * area / (r->perimeter_length * r->perimeter_length);
	}
	return DFX_SUCCESS;
}

/*!
 *  \brief Compute the properties of one region of a label plane (see the dfx_region structure).
 *
 *  Area, bounding box, centroid and the central second-order moments are computed from the pixel
 *  coordinates (Birchfield, chapter 4):
 *
 *     m_pq = sum x^p y^q,  (x_c, y_c) = (m10 / m00, m01 / m00),
 *     mu20 = m20 - x_c m10,  mu02 = m02 - y_c m01,  mu11 = m11 - x_c m01.
 *
 *  The equivalent ellipse is the ellipse with the same second-order moments. Its orientation is
 *  theta = 1/2 atan2(2 mu11, mu20 - mu02), the angle from the x axis towards the y axis (clockwise
 *  on the screen); with S = sqrt((mu20 - mu02)^2 + 4 mu11^2) and T = mu20 + mu02, its semi-axes
 *  are a = sqrt(2 (T + S) / m00) and b = sqrt(2 (T - S) / m00), and its eccentricity
 *  sqrt(1 - b^2 / a^2) = sqrt(2 S / (T + S)). (For digitized disks and ellipses a and b come out within
 *  a few hundredths of a pixel of the true semi-axes.)
 *
 *  The perimeter is given in two ways: as the number of boundary pixels (region pixels having a
 *  4-neighbor outside the region or outside the image; unlike boundary_binary() on a replicate-padded
 *  image, this counts the pixels along the image border), and as the length of the outer boundary,
 *  estimated by Kimura's formula from the 8-connected chain code of the region traced from its
 *  first pixel (see trace_boundary_binary() and chain_code_length()). The compactness
 *  4 pi area / perimeter_length^2 is 1 for a disk and smaller for every other shape; for digitized
 *  disks it comes out slightly above 1 (1.00 to 1.07 for radii from 3 to 60), since Kimura's formula
 *  underestimates the length of a digital circle a little.
 *
 *  \param[in]      L      - label plane
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter
 *  \param[in]      label  - the region to measure (> 0)
 *  \param[out]     r      - receives the properties
 *  \return DFX_INVRANGE if no pixel carries the label.
 */
int region_properties(const int* L, int width, int height, int p, int label, dfx_region* r)
{
	unsigned char* A = NULL;
	int err;
	if (L == NULL || r == NULL || height < 0 || width < 0 || p < 0 || label <= 0) return DFX_INVARG;
	if ((err = alloc_binary(&A, width, height, p)) != DFX_SUCCESS) return err;
	if ((err = select_label(L, A, width, height, p, label)) == DFX_SUCCESS)
		err = binary_properties(A, width, height, p, r);
	free_binary(A);
	return err;
}

/*!
 *  \brief Compute the properties of the foreground of a binary image, treated as a single region.
 *
 *  Same as region_properties(); the outer boundary traced for perimeter_length is that of the
 *  component containing the first ON pixel. Returns DFX_INVRANGE if the image has no ON pixel.
 */
int region_properties_binary(const unsigned char* A, int width, int height, int p, dfx_region* r)
{
	if (A == NULL || r == NULL || height < 0 || width < 0 || p < 0) return DFX_INVARG;
	return binary_properties(A, width, height, p, r);
}

/**************
 *
 *  Distances and region growing:
 *
 *   chamfer_binary()
 *   grow_region_plane()
 */

/*!
 *  \brief Chamfer distance transform: the distance from every pixel to the nearest ON pixel.
 *
 *  Implements Birchfield's Algorithm 4.10. The distance between neighboring pixels is a for
 *  horizontal/vertical steps and b for diagonal ones, so the result is the length of the shortest
 *  path of such steps to the foreground: (a,b) = (1,1) gives the chessboard (D8) distance,
 *  b >= 2a the city-block (D4) distance, and (1, 1.35) or (3,4) a good approximation of the Euclidean
 *  distance. Two passes suffice: a forward pass propagating distances from the left and upper
 *  neighbors, and a backward pass from the right and lower ones. ON pixels have distance 0.
 *  The output plane receives PAD_REPLICATE padding.
 *
 *  \param[in]      A      - input binary image
 *  \param[out]     D      - output plane of distances
 *  \param[in]      width  - image width
 *  \param[in]      height - image height
 *  \param[in]      p      - padding parameter of both
 *  \param[in]      a      - length of a horizontal or vertical step (> 0)
 *  \param[in]      b      - length of a diagonal step (> 0)
 *  \return DFX_INVRANGE if the image has no ON pixel (every distance is then HUGE_VALF).
 */
int chamfer_binary(const unsigned char* A, float* D, int width, int height, int p, float a, float b)
{
	int width_p = width + 2 * p;
	int x, y, count;
	float d;

	/* check parameters: */
	if (A == NULL || D == NULL || height < 0 || width < 0 || p < 0 || !(a > 0) || !(b > 0)) return DFX_INVARG;

	/* first pass: */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		if (A[(p + y) * width_p + p + x]) { D[(p + y) * width_p + p + x] = 0; continue; }
		d = HUGE_VALF;
		if (x > 0 && a + D[(p + y) * width_p + p + x - 1] < d) d = a + D[(p + y) * width_p + p + x - 1];
		if (y > 0) {
			if (a + D[(p + y - 1) * width_p + p + x] < d) d = a + D[(p + y - 1) * width_p + p + x];
			if (x > 0 && b + D[(p + y - 1) * width_p + p + x - 1] < d) d = b + D[(p + y - 1) * width_p + p + x - 1];
			if (x < width - 1 && b + D[(p + y - 1) * width_p + p + x + 1] < d) d = b + D[(p + y - 1) * width_p + p + x + 1];
		}
		D[(p + y) * width_p + p + x] = d;
	}
	/* second pass: */
	for (y = height - 1; y >= 0; y--) for (x = width - 1; x >= 0; x--) {
		if (A[(p + y) * width_p + p + x]) continue;
		d = D[(p + y) * width_p + p + x];
		if (x < width - 1 && a + D[(p + y) * width_p + p + x + 1] < d) d = a + D[(p + y) * width_p + p + x + 1];
		if (y < height - 1) {
			if (a + D[(p + y + 1) * width_p + p + x] < d) d = a + D[(p + y + 1) * width_p + p + x];
			if (x < width - 1 && b + D[(p + y + 1) * width_p + p + x + 1] < d) d = b + D[(p + y + 1) * width_p + p + x + 1];
			if (x > 0 && b + D[(p + y + 1) * width_p + p + x - 1] < d) d = b + D[(p + y + 1) * width_p + p + x - 1];
		}
		D[(p + y) * width_p + p + x] = d;
	}
	pad_plane(D, width, height, p, PAD_REPLICATE);
	count_binary(A, width, height, p, &count);
	return count > 0 ? DFX_SUCCESS : DFX_INVRANGE;
}

/*!
 *  \brief Seeded region growing on a plane.
 *
 *  Starting from the seed pixel, neighboring pixels are added to the region as long as their values
 *  differ from the mean value of the region grown so far by at most the tolerance; the mean is
 *  updated after every addition. The region is grown breadth-first (a pixel rejected at one moment
 *  may still join later, through another neighbor, once the mean has moved). The output receives
 *  PAD_REPLICATE padding.
 *
 *  \param[in]      X            - input plane (e.g. luminance)
 *  \param[out]     A            - output binary image: the region grown
 *  \param[in]      width        - image width
 *  \param[in]      height       - image height
 *  \param[in]      p            - padding parameter of both
 *  \param[in]      x0, y0       - seed pixel
 *  \param[in]      tolerance    - maximum allowed deviation from the region mean (>= 0)
 *  \param[in]      connectivity - 4 or 8
 */
int grow_region_plane(const float* X, unsigned char* A, int width, int height, int p, int x0, int y0, float tolerance, int connectivity)
{
	int width_p = width + 2 * p;
	int k, dirs[8], n_dirs, count;
	int* queue = NULL;
	size_t head = 0, tail = 0, capacity;
	double sum;

	/* check parameters: */
	if (X == NULL || A == NULL || height < 0 || width < 0 || p < 0 || !(tolerance >= 0)) return DFX_INVARG;
	if (x0 < 0 || x0 >= width || y0 < 0 || y0 >= height) return DFX_INVARG;
	if ((n_dirs = neighbor_directions(connectivity, dirs)) == 0) return DFX_INVARG;

	/* FIFO queue; every pixel is enqueued at most once */
	capacity = (size_t)width * height + 1;
	if ((queue = (int*)malloc(capacity * sizeof(int))) == NULL) return DFX_NOMEM;
	fill_binary(A, width, height, p, 0);

	/* seed: */
	A[(p + y0) * width_p + p + x0] = 1;
	queue[tail++] = y0 * width + x0;
	sum = X[(p + y0) * width_p + p + x0];
	count = 1;

	/* grow: */
	while (head != tail) {
		int q = queue[head], qx = q % width, qy = q / width;
		head = (head + 1) % capacity;
		for (k = 0; k < n_dirs; k++) {
			int nx = qx + DX[dirs[k]], ny = qy + DY[dirs[k]], i;
			if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;
			i = (p + ny) * width_p + p + nx;
			if (A[i]) continue;
			if (fabs((double)X[i] - sum / count) > (double)tolerance) continue;
			A[i] = 1;
			sum += X[i];
			count++;
			queue[tail] = ny * width + nx;
			tail = (tail + 1) % capacity;
		}
	}
	free(queue);
	return pad_binary(A, width, height, p, PAD_REPLICATE);
}
