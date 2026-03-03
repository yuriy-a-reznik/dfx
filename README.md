# dfx
DFX: Digital image processing Functions and eXamples.
 
This is a collection of basic image processing primitives implemented entirely in C, 
without reliance on any external libraries. 
 
   *  Functions: read_bitmap() / write_bitmap() load and store 24-bit-packed sRGB images. 
      In memory, such images are presented as unsigned character arrays. Parameters "width" 
      and "height" define image dimensions. 
      
   *  Functions: srgb_to_linear() and linear_to_srgb() perform conversions between sRGB and Linear RGB format. 
      In linear RGB form, images are stored in memory as 3 arrays of floats, called "R" (red), "G" (green), 
      and "B" (blue) planes. Linear representation is fundamental. Only in linear space we can mix colors. 
      All color and image processing operations are implemented in the linear RGB space.
      
   *  When images are converted from sRGB format to linear, they can be padded. 
      Parameter "p" defines the number of pixels to be reserved around the boundary of the original image. 
      Such padding becomes handy in implementing filtering and scaling operations. 
 
   *  Functions linear_to_luminance() and luminance_to_grayscale_image() allow the extraction of luminance 
      (Y channel in CIE 1931 XYZ space), and translation of luminance to a grayscale RGB image. 
 
   *  Functions filter_plane() and filter_image() implement various image filtering operations. 
 
   *  Functions dft_plane() and dft_magnitude() and dft_phase() compute Discrete Fourier Transform over an image plane.
 
 Examples showing how to use these operations can be found in the directory dfx/demos. 
