# Changelog

## 1.2 (2026-11-08)

Binary images and regions (new modules `dfx_binary.c` and `dfx_regions.c`)
- A binary image is a plane of unsigned chars holding 0 (OFF) or 1 (ON), with the same
  padded geometry as a float plane; connected-component labels are planes of ints. The
  set operations of morphology are kept apart from the arithmetic of linear planes: a
  plane becomes a binary image only through thresholding, and a binary image becomes a
  plane again only through `binary_to_plane()`.
- Thresholding: `threshold_plane()`, `histogram_plane()`, `otsu_threshold()`,
  `ridler_calvard_threshold()`, `niblack_threshold_plane()`, `double_threshold_plane()`.
- Morphology: `structuring_element()` (cross, box, disk of any radius), `dilate_binary()`,
  `erode_binary()`, `open_binary()`, `close_binary()`, `majority_binary()`,
  `hitmiss_binary()` (ternary structuring elements), `thin_binary()`, `thicken_binary()`,
  `boundary_binary()`; set operations `not/and/or/xor/subtract_binary()`, `count_binary()`.
- Regions: `floodfill_labels()`, `label_binary()` (union-find) and `label_binary_floodfill()`,
  `select_label()`, `labels_to_image()`, `reconstruct_binary()`, `fill_holes_binary()`,
  `region_properties()` (area, bounding box, centroid, central moments, orientation,
  eccentricity, equivalent ellipse, perimeter, compactness), `trace_boundary_binary()`
  (Moore tracing and wall following, Freeman chain codes), `chain_code_length()`
  (Freeman, Pythagorean, Kimura), `chamfer_binary()`, `grow_region_plane()`.
- New `dfx_binary_demo` (also a Visual Studio project): synthesizes a test scene, runs the
  whole chain, and checks the measurements against the known shapes; `make check` runs it.

## 1.1 (2026-11-08)

Portability
- Builds with strict `-std=c11` on gcc and clang (M_PI / M_LN2 fallbacks in `dfx.h`).
- All sources are plain ASCII.
- Bitmap reader: the "BM" signature is checked, `bfOffBits` is honored, V4/V5 info
  headers, top-down bitmaps and zero pixel density are accepted. The headers are still
  read as packed structures; big-endian hosts are rejected at compile time, and the
  structure sizes are checked at compile time.
- Makefile dependencies fixed; `make check` runs the self-tests.

Fixes
- `dft_magnitude()` and `dft_phase()` used the real part in place of the imaginary part.
- `linear_to_srgb_dithered()` diffused the red error into the blue channel, skipped
  the last row and column, and modified its inputs.
- `adjust_hue()` rejected theta = M_PI; `luv_sat_scale()` read past the color table
  for `colors == N_COLORS`; the RGB range check reported round-trip rounding errors.
- `dct2_rec()` / `dct4_rec()` support n = 1 and 2; the 7-point DCT-II constants are
  accurate to double precision.
- `snr_image()` pools energies over planes and accumulates in double.
- Error codes are propagated by all `_image` wrappers and by the filter and resampler.
- `dfx_color_adj_demo` no longer overflows its buffers for inputs larger than 1120 x 630.
- Citation years corrected (Crochiere-Rabiner 1981, Burrus-Parks 1985).

API
- All public functions documented in `dfx.h`; input pointers are `const`.
- `xyz_to_lyv()` / `lyv_to_xyz()` renamed to `xyz_to_luv()` / `luv_to_xyz()`.
- `xyz_to_lab()` and `lab_to_xyz()` implemented (ISO/CIE 11664-4).
- New `clip_plane()`, `clip_image()`, `fft()`, `fft_dit()`, `fft_dif()`, `fft_pfa()`;
  `dct2_rec()` and `dct4_rec()` declared and returning `DFX_*` codes.
- Filtering, resampling and color conversions replicate-pad their outputs.

Examples
- New `dfx_fft_demo` and `dfx_dct_demo` self-tests (also Visual Studio projects).
- Demos check return codes; the color demo reports out-of-gamut results.

## 1.02 (2026-03-15)

Initial public version.
