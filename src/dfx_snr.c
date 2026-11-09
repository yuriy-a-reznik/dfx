/*!
 *  \file   dfx_snr.c
 *  \brief  Signal to Noise Ratio (SNR) and related metrics.
 *
 *  For a reference plane X_ref and a test plane X_test with W x H pixels, this module computes
 *
 *      MSE = (1 / WH) sum (X_ref - X_test)^2,
 *      SNR = 10 log10( E_ref / MSE ),   where  E_ref = (1 / WH) sum X_ref^2,
 *
 *  i.e. the ratio of the mean energy of the reference signal to the mean energy of the error,
 *  in decibels. For an RGB image, the squared errors and energies of the three planes are pooled
 *  before taking the ratio. SNR is capped at 100 dB when the error vanishes.
 *
 *  Note that this is the SNR of linear-light values, which is not the PSNR commonly reported for
 *  8-bit images (PSNR uses the peak value 255 instead of E_ref, and gamma-companded samples).
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    November 8, 2026
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "dfx.h"

/* SNR reported when the error energy is zero (or negligible): */
#define DFX_SNR_MAX   100.0

/*!
 *  \brief Accumulate squared errors and reference energies over the interior of two planes.
 *
 *  \param[in]   X_ref   - reference image plane (padded)
 *  \param[in]   X_test  - test image plane (padded)
 *  \param[in]   width   - image width in pixels
 *  \param[in]   height  - image height in pixels
 *  \param[in]   p       - padding parameter
 *  \param[out]  p_se    - variable to receive the sum of squared errors
 *  \param[out]  p_e_ref - variable to receive the sum of squared reference values
 *
 *  \return      DFX_SUCCESS if success, or an error code otherwise.
 */
static int accumulate_errors(const float* X_ref, const float* X_test, int width, int height, int p, double* p_se, double* p_e_ref)
{
	double se = 0., e_ref = 0., ref, delta;
	int width_p = width + 2 * p;	/* width of padded plane */
	int x, y;

	/* check parameters: */
	if (X_ref == NULL || X_test == NULL || p_se == NULL || p_e_ref == NULL || height < 0 || width < 0 || p < 0)
		return DFX_INVARG;

	/* accumulate energies (in double precision, since the sums may run over millions of pixels): */
	for (y = 0; y < height; y++) for (x = 0; x < width; x++) {
		ref = X_ref[(p + y) * width_p + p + x];
		delta = ref - X_test[(p + y) * width_p + p + x];
		se += delta * delta;
		e_ref += ref * ref;
	}

	*p_se = se;
	*p_e_ref = e_ref;
	return DFX_SUCCESS;
}

/*!
 *  \brief Compute SNR (in dB) from pooled sums of squared errors and reference values.
 */
static float snr_from_sums(double se, double e_ref)
{
	if (se <= 0. || e_ref > 1e10 * se) return (float)DFX_SNR_MAX;	/* no error (or practically none) */
	if (e_ref <= 0.) return 0.f;									/* zero reference: SNR is undefined, report 0 dB */
	return (float)(10.0 * log10(e_ref / se));
}

 /*!
  *  \brief Compute MSE and SNR metrics between two image planes.
  *
  *  \param[in]      X_ref   - reference image plane (padded)
  *  \param[in]      X_test  - test image plane (padded)
  *  \param[in]      width   - image width in pixels
  *  \param[in]      height  - image height in pixels
  *  \param[in]      p       - padding parameter
  *  \param[out]     p_MSE   - pointer to where to store the computed MSE value
  *  \param[out]     p_SNR   - pointer to where to store the computed SNR value (dB)
  *
  *  \return         DFX_SUCCESS if success, or an error code otherwise.
  */
int snr_plane(const float* X_ref, const float* X_test, int width, int height, int p, float *p_MSE, float *p_SNR)
{
	double se, e_ref, n;
	int err;

	/* check parameters: */
	if (p_MSE == NULL || p_SNR == NULL) return DFX_INVARG;
	if ((err = accumulate_errors(X_ref, X_test, width, height, p, &se, &e_ref)) != DFX_SUCCESS) return err;

	/* normalize MSE by the number of pixels, compute SNR: */
	n = (double)width * (double)height;
	*p_MSE = (float)(n > 0. ? se / n : 0.);
	*p_SNR = snr_from_sums(se, e_ref);
	return DFX_SUCCESS;
}

/*!
 *  \brief Compute MSE and SNR metrics between two images.
 *
 *  Squared errors and reference energies are pooled over the three planes, so the result
 *  is the MSE/SNR of the image treated as a single signal with 3WH samples.
 */
int snr_image(const float* R_ref, const float* G_ref, const float* B_ref, const float* R_test, const float* G_test, const float* B_test, int width, int height, int p, float* p_MSE, float* p_SNR)
{
	double se_R, se_G, se_B, e_R, e_G, e_B, se, e_ref, n;
	int err;

	/* check parameters: */
	if (p_MSE == NULL || p_SNR == NULL) return DFX_INVARG;

	/* accumulate errors and energies for each plane: */
	if ((err = accumulate_errors(R_ref, R_test, width, height, p, &se_R, &e_R)) != DFX_SUCCESS) return err;
	if ((err = accumulate_errors(G_ref, G_test, width, height, p, &se_G, &e_G)) != DFX_SUCCESS) return err;
	if ((err = accumulate_errors(B_ref, B_test, width, height, p, &se_B, &e_B)) != DFX_SUCCESS) return err;

	/* pool across planes: */
	se = se_R + se_G + se_B;
	e_ref = e_R + e_G + e_B;
	n = 3.0 * (double)width * (double)height;

	/* store the results: */
	*p_MSE = (float)(n > 0. ? se / n : 0.);
	*p_SNR = snr_from_sums(se, e_ref);
	return DFX_SUCCESS;
}

/* dfx_snr.c -- end of file */
