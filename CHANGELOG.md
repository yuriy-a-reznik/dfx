# Changelog

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
