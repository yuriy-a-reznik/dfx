# dfx

DFX: Digital image processing Functions and eXamples. Version 1.2.

This library provides a collection of basic image-processing primitives
implemented entirely in C, without relying on external libraries. It is
written as a companion to an image-processing course: every operation
follows the textbook definition rather than a production shortcut, and all
color, filtering, and resampling operations are performed on linear-light
(not gamma-companded) pixel values.

## What is in the library

* `read_bitmap()` and `write_bitmap()` load and store 24-bit packed sRGB
  images (`.bmp` files). In memory, these images are represented as arrays
  of unsigned chars. The parameters `width` and `height` specify the image
  dimensions.

* `srgb_to_linear()` and `linear_to_srgb()` convert between sRGB and linear
  RGB. In linear RGB form, images are stored as three float planes: R (red),
  G (green), and B (blue). Linear representation is fundamental: all color
  mixing and image-processing operations in this library are implemented in
  linear space. `linear_to_srgb_dithered()` adds Floyd-Steinberg error
  diffusion to the 8-bit quantization.

* When converting from sRGB to linear RGB, images are padded. The parameter
  `p` defines the number of pixels reserved around the boundary of the
  original image. Padding supplies the out-of-image samples needed by
  filtering and resampling operations; `pad_plane()` and `pad_image()` fill
  it by zeros, by replication, or by reflection.

* `linear_to_luminance()` and `luminance_to_grayscale_image()` extract image
  luminance (the Y channel in CIE 1931 XYZ space) and convert a luminance
  plane into a grayscale RGB image.

* `filter_plane()` and `filter_image()` implement separable low-pass
  filtering with truncated-sinc, Lanczos, and Gaussian kernels.

* `resample_plane()` and `resample_image()` implement rational-factor image
  resampling with polyphase filters.

* `dft_plane()`, `dft_magnitude()`, and `dft_phase()` compute the 2D Discrete
  Fourier Transform of an image plane by direct evaluation and extract
  magnitude and phase. `fft()` (with `fft_dit()`, `fft_dif()`, and
  `fft_pfa()`) and `dct2_rec()`, `dct4_rec()` implement fast 1D Fourier and
  cosine transforms.

* `linear_rgb_to_xyz()`, `xyz_to_luv()`, `xyz_to_lab()` and their inverses
  convert between linear RGB (for several sets of primaries), CIE XYZ,
  CIE L\*u\*v\*, and CIE L\*a\*b\*. `adjust_brightness()`,
  `adjust_saturation()`, and `adjust_hue()` perform color adjustments in
  L\*u\*v\* space.

* `snr_plane()` and `snr_image()` compute MSE and SNR between two images.

* Binary images (sets of pixels, stored as planes of 0/1 values with the
  same padded geometry as linear planes) are produced from luminance or any
  other plane by `threshold_plane()`, with global thresholds from
  `otsu_threshold()` and `ridler_calvard_threshold()`, locally adaptive
  thresholds from `niblack_threshold_plane()`, or hysteresis from
  `double_threshold_plane()`. `dilate_binary()`, `erode_binary()`,
  `open_binary()`, `close_binary()`, `majority_binary()`,
  `hitmiss_binary()`, `thin_binary()`, `thicken_binary()`, and
  `boundary_binary()` implement the morphological operations with
  structuring elements from `structuring_element()`.

* `label_binary()` (union-find) and `label_binary_floodfill()` find the
  connected components of a binary image; `floodfill_labels()`,
  `reconstruct_binary()`, and `fill_holes_binary()` are the related
  connectivity operations. `region_properties()` measures a region (area,
  centroid, moments, orientation, eccentricity, perimeter, compactness),
  `trace_boundary_binary()` and `chain_code_length()` trace boundaries as
  Freeman chain codes, `chamfer_binary()` computes distance transforms, and
  `grow_region_plane()` performs seeded region growing.

The public interface, with the conventions for image layout and padding, is
documented in `include/dfx.h`.

## Building

On Linux, macOS, or any system with a C compiler and `make`:

    make            # builds build/lib/libdfx.a and the examples in build/bin
    make check      # runs the self-tests (FFT, DCT, zone plate, binary images)

On Windows, open `build/dfx.sln` in Visual Studio.

The code is standard C (C99/C11) and compiles cleanly with
`-std=c11 -Wall -Wextra -pedantic`.

## Examples

Examples demonstrating how to use these operations can be found in the
`examples` directory:

| Program                | What it shows                                                        |
|------------------------|----------------------------------------------------------------------|
| `dfx_rgb_demo`         | splitting an image into R, G, B, and luminance components            |
| `dfx_gradient_demo`    | linear vs. gamma-companded gradients, banding, and dithering         |
| `dfx_zone_plate_demo`  | generating the zone-plate test pattern                               |
| `dfx_filter_demo`      | low-pass filtering with sinc, Lanczos, and Gaussian kernels          |
| `dfx_aliasing_demo`    | aliasing from subsampling, with and without a pre-filter             |
| `dfx_resampler_demo`   | resampling by rational factors                                       |
| `dfx_dft_demo`         | 2D DFT magnitude and phase images                                    |
| `dfx_fft_demo`         | verifying and timing the fast FFT algorithms                         |
| `dfx_dct_demo`         | verifying the fast DCT-II and DCT-IV algorithms                      |
| `dfx_color_adj_demo`   | brightness, saturation, and hue adjustments in CIE L\*u\*v\* space    |
| `dfx_binary_demo`      | thresholding, morphology, connected components, region properties, chamfer distances, skeletons |

The programs that take an image as input accept 24-bit `.bmp` files; test
images (the Kodak set at several resolutions, and the zone plate) are in the
`test` directory.

## Changes

See `CHANGELOG.md`.
