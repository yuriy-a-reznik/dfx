# dfx
This library provides a collection of basic image‑processing primitives implemented entirely in C, without relying on external libraries.

* read_bitmap() and write_bitmap() load and store 24‑bit packed sRGB images. In memory, these images are represented as arrays of unsigned chars. The parameters “width” and “height” specify the image dimensions.

* srgb_to_linear() and linear_to_srgb() convert between sRGB and linear RGB. In linear RGB form, images are stored as three float arrays: R (red), G (green), and B (blue). Linear representation is fundamental. All color mixing and image‑processing operations in this library are implemented in linear space.

* When converting from sRGB to linear RGB, images may be padded. The parameter “p” defines the number of pixels reserved around the boundary of the original image. Padding is useful for filtering and resampling operations.
 
* linear_to_luminance() and luminance_to_grayscale_image() extract image luminance (the Y channel in CIE 1931 XYZ space) and convert luminance channel into a grayscale RGB image.

* filter_plane() and filter_image() implement a variety of image‑filtering operations.

* resample_plane() and resample_image() implement image‑resampling operations.

*  dft_plane(), dft_magnitude(), and dft_phase() compute the Discrete Fourier Transform of an image plane and extract magnitude and phase.

Examples demonstrating how to use these operations can be found in the dfx/examples directory.
