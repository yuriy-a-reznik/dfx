/*!
 *  \file   dfx_fft.c
 *  \brief  DFX library: implementations of FFT transforms of various sizes.
 *
 *  This module includes examples of:
 * 
 *     Short length FFT/WFTA modules (n=2,3,4,5,7,8,9,16)
 *     PFA method (in-place construction)
 *     DIT and DIF methods
 *
 *  Reference:
 *    C.S.Burrus and T.W.Parks, "DFT/FFT and convolution algorithms: theory and implementation," John Wiley & Sons, Inc., 1991.
 * 
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.0
 *  \date    February 20, 2026
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "dfx.h"

/**********************
 * 
 * Short-length FFT/WFTA modules:
 * 
 *  wfta_2()
 *  wfta_3()
 *  wfta_4()
 *  wfta_5()
 *  wfta_7()
 *  wfta_8()
 *  wfta_16()
 */

/*!
 *  \brief In-place WFTA module for n = 2.
 *
 *  \param xr   Real part array of length N.
 *  \param xi   Imag part array of length N.
 *  \param idx  Index map for input positions.
 *  \param ip   Index map for output positions.
 */
static void wfta_2(float* xr, float* xi, int * idx, int * ip)
{
    float r1;

    r1 = xr[idx[0]];
    xr[idx[0]] = r1 + xr[idx[1]];
    xr[idx[1]] = r1 - xr[idx[1]];

    r1 = xi[idx[0]];
    xi[ip[0]] = r1 + xi[idx[1]];
    xi[ip[1]] = r1 - xi[idx[1]];
}

/*!
 *  \brief In-place WFTA module for n = 3.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_3(float* xr, float* xi, int * idx, int * ip)
{
    const float C31 = -0.86602540f;   /* cos(2*M_PI/3) */
    const float C32 = -1.5f;
    
    float r1, r2, s1, s2;

    r2 = (xr[idx[1]] - xr[idx[2]]) * C31;
    r1 = xr[idx[1]] + xr[idx[2]];
    xr[idx[0]] = xr[idx[0]] + r1;
    r1 = xr[idx[0]] + r1 * C32;

    s2 = (xi[idx[1]] - xi[idx[2]]) * C31;
    s1 = xi[idx[1]] + xi[idx[2]];
    xi[idx[0]] = xi[idx[0]] + s1;
    s1 = xi[idx[0]] + s1 * C32;

    xr[ip[1]] = r1 - s2;
    xr[ip[2]] = r1 + s2;
    xi[ip[1]] = s1 + r2;
    xi[ip[2]] = s1 - r2;
}

/*!
 *  \brief In-place WFTA module for n1 = 4.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_4(float* xr, float* xi, int * idx, int * ip)
{
    float r1, r2, t1, t2;

    r1 = xr[idx[0]] + xr[idx[2]];
    t1 = xr[idx[0]] - xr[idx[2]];
    r2 = xr[idx[1]] + xr[idx[3]];
    xr[ip[0]] = r1 + r2;
    xr[ip[2]] = r1 - r2;

    r1 = xi[idx[0]] + xi[idx[2]];
    t2 = xi[idx[0]] - xi[idx[2]];
    r2 = xi[idx[1]] + xi[idx[3]];
    xi[ip[0]] = r1 + r2;
    xi[ip[2]] = r1 - r2;

    r1 = xr[idx[1]] - xr[idx[3]];
    r2 = xi[idx[1]] - xi[idx[3]];

    xr[ip[1]] = t1 + r2;
    xr[ip[3]] = t1 - r2;
    xi[ip[1]] = t2 - r1;
    xi[ip[3]] = t2 + r1;
}

/*!
 *  \brief In-place WFTA module for n1 = 5.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_5(float* xr, float* xi, int * idx, int * ip)
{
    const float C51 = 0.95105652f;    /* sin(2*M_PI/5) */
    const float C52 = -1.53884180f;
    const float C53 = -0.36327126f;
    const float C54 = 0.55901699f;
    const float C55 = -1.25f;

    float r1, r2, r3, r4, s1, s2, s3, s4, t;

    r1 = xr[idx[1]] + xr[idx[4]];
    r4 = xr[idx[1]] - xr[idx[4]];
    r3 = xr[idx[2]] + xr[idx[3]];
    r2 = xr[idx[2]] - xr[idx[3]];

    t = (r1 - r3) * C54;
    r1 = r1 + r3;
    xr[idx[0]] = xr[idx[0]] + r1;
    r1 = xr[idx[0]] + r1 * C55;

    r3 = r1 - t;
    r1 = r1 + t;

    t = (r4 + r2) * C51;
    r4 = t + r4 * C52;
    r2 = t + r2 * C53;

    s1 = xi[idx[1]] + xi[idx[4]];
    s4 = xi[idx[1]] - xi[idx[4]];
    s3 = xi[idx[2]] + xi[idx[3]];
    s2 = xi[idx[2]] - xi[idx[3]];

    t = (s1 - s3) * C54;
    s1 = s1 + s3;
    xi[idx[0]] = xi[idx[0]] + s1;
    s1 = xi[idx[0]] + s1 * C55;

    s3 = s1 - t;
    s1 = s1 + t;

    t = (s4 + s2) * C51;
    s4 = t + s4 * C52;
    s2 = t + s2 * C53;

    xr[ip[1]] = r1 + s2;
    xr[ip[4]] = r1 - s2;
    xr[ip[2]] = r3 - s4;
    xr[ip[3]] = r3 + s4;

    xi[ip[1]] = s1 - r2;
    xi[ip[4]] = s1 + r2;
    xi[ip[2]] = s3 + r4;
    xi[ip[3]] = s3 - r4;
}

/*!
 *  \brief In-place WFTA module for n1 = 7.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map (I in Fortran).
 *  \param ip   Output index map (IP in Fortran).
 */
static void wfta_7(float* xr, float* xi, int * idx, int * ip)
{
    const float C71 = -1.16666667f;
    const float C72 = -0.79015647f;
    const float C73 = 0.055854267f;
    const float C74 = 0.7343022f;
    const float C75 = 0.44095855f;
    const float C76 = -0.34087293f;
    const float C77 = 0.53396936f;
    const float C78 = 0.87484229f;

    float r1, r2, r3, r4, r5, r6;
    float s1, s2, s3, s4, s5, s6;
    float t, t3;

    r1 = xr[idx[1]] + xr[idx[6]];
    r6 = xr[idx[1]] - xr[idx[6]];
    s1 = xi[idx[1]] + xi[idx[6]];
    s6 = xi[idx[1]] - xi[idx[6]];

    r2 = xr[idx[2]] + xr[idx[5]];
    r5 = xr[idx[2]] - xr[idx[5]];
    s2 = xi[idx[2]] + xi[idx[5]];
    s5 = xi[idx[2]] - xi[idx[5]];

    r3 = xr[idx[3]] + xr[idx[4]];
    r4 = xr[idx[3]] - xr[idx[4]];
    s3 = xi[idx[3]] + xi[idx[4]];
    s4 = xi[idx[3]] - xi[idx[4]];

    t3 = (r1 - r2) * C74;
    t = (r1 - r3) * C72;
    r1 = r1 + r2 + r3;
    xr[idx[0]] = xr[idx[0]] + r1;
    r1 = xr[idx[0]] + r1 * C71;
    r2 = (r3 - r2) * C73;
    r3 = r1 - t + r2;
    r2 = r1 - r2 - t3;
    r1 = r1 + t + t3;

    t = (r6 - r5) * C78;
    t3 = (r6 + r4) * C76;
    r6 = (r6 + r5 - r4) * C75;
    r5 = (r5 + r4) * C77;
    r4 = r6 - t3 + r5;
    r5 = r6 - r5 - t;
    r6 = r6 + t3 + t;

    t3 = (s1 - s2) * C74;
    t = (s1 - s3) * C72;
    s1 = s1 + s2 + s3;
    xi[idx[0]] = xi[idx[0]] + s1;
    s1 = xi[idx[0]] + s1 * C71;
    s2 = (s3 - s2) * C73;
    s3 = s1 - t + s2;
    s2 = s1 - s2 - t3;
    s1 = s1 + t + t3;

    t = (s6 - s5) * C78;
    t3 = (s6 + s4) * C76;
    s6 = (s6 + s5 - s4) * C75;
    s5 = (s5 + s4) * C77;
    s4 = s6 - t3 + s5;
    s5 = s6 - s5 - t;
    s6 = s6 + t3 + t;

    xr[ip[1]] = r3 + s4;  
    xr[ip[6]] = r3 - s4;  
    xr[ip[2]] = r1 + s6;  
    xr[ip[5]] = r1 - s6;  
    xr[ip[3]] = r2 - s5;  
    xr[ip[4]] = r2 + s5;  

    xi[ip[3]] = s2 + r5;  
    xi[ip[4]] = s2 - r5;  
    xi[ip[1]] = s3 - r4;  
    xi[ip[6]] = s3 + r4;  
    xi[ip[2]] = s1 - r6;  
    xi[ip[5]] = s1 + r6;  
}

/*!
 *  \brief In-place WFTA module for n1 = 8.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_8(float* xr, float* xi, int * idx, int * ip)
{
    const float C81 = 0.70710678f;    /* 1/sqrt(2) = cos(*M_PI/4) = sin(*M_PI/4) */

    float r1, r2, r3, r4, r5, r6, r7, r8;
    float s1, s2, s3;
    float t1, t2, t3, t4;

    r1 = xr[idx[0]] + xr[idx[4]];
    r2 = xr[idx[0]] - xr[idx[4]];
    r3 = xr[idx[1]] + xr[idx[7]];
    r4 = xr[idx[1]] - xr[idx[7]];
    r5 = xr[idx[2]] + xr[idx[6]];
    r6 = xr[idx[2]] - xr[idx[6]];
    r7 = xr[idx[3]] + xr[idx[5]];
    r8 = xr[idx[3]] - xr[idx[5]];

    t1 = r1 + r5;
    t2 = r1 - r5;
    t3 = r3 + r7;
    r3 = (r3 - r7) * C81;

    xr[ip[0]] = t1 + t3;
    xr[ip[4]] = t1 - t3;

    t1 = r2 + r3;
    t3 = r2 - r3;
    s1 = r4 - r8;
    r4 = (r4 + r8) * C81;
    s2 = r4 + r6;
    s3 = r4 - r6;

    r1 = xi[idx[0]] + xi[idx[4]];
    r2 = xi[idx[0]] - xi[idx[4]];
    r3 = xi[idx[1]] + xi[idx[7]];
    r4 = xi[idx[1]] - xi[idx[7]];
    r5 = xi[idx[2]] + xi[idx[6]];
    r6 = xi[idx[2]] - xi[idx[6]];
    r7 = xi[idx[3]] + xi[idx[5]];
    r8 = xi[idx[3]] - xi[idx[5]];

    t4 = r1 + r5;
    r1 = r1 - r5;
    r5 = r3 + r7;
    r3 = (r3 - r7) * C81;

    xi[ip[0]] = t4 + r5;
    xi[ip[4]] = t4 - r5;

    r5 = r2 + r3;
    r2 = r2 - r3;
    r3 = r4 - r8;
    r4 = (r4 + r8) * C81;
    r7 = r4 + r6;
    r4 = r4 - r6;

    xr[ip[1]] = t1 + r7;
    xr[ip[7]] = t1 - r7;
    xr[ip[2]] = t2 + r3;
    xr[ip[6]] = t2 - r3;
    xr[ip[3]] = t3 + r4;
    xr[ip[5]] = t3 - r4;

    xi[ip[1]] = r5 - s2;
    xi[ip[7]] = r5 + s2;
    xi[ip[2]] = r1 - s1;
    xi[ip[6]] = r1 + s1;
    xi[ip[3]] = r2 - s3;
    xi[ip[5]] = r2 + s3;
}

/*!
 *  \brief In-place WFTA module for n1 = 9.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_9(float* xr, float* xi, int * idx, int * ip)
{
    const float C31 = -0.86602540f;   /* cos(2*M_PI/3) */
    const float C95 = -0.50000000f;   /* -1/2 */
    const float C92 = 0.93969262f;    /* cos(M_PI/9) */
    const float C93 = -0.17364818f;   /* -sin(M_PI/18) */
    const float C94 = 0.76604444f;    /* cos(2*M_PI/9)  */
    const float C96 = -0.34202014f;   /* -sin(M_PI/9) */
    const float C97 = -0.98480775f;   /* -cos(M_PI/18) */
    const float C98 = -0.64278761f;   /* -sin(2*M_PI/9) */

    float r1, r2, r3, r4, r5, r6, r7, r8, r9;
    float t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;

    r1 = xr[idx[1]] + xr[idx[8]];
    r2 = xr[idx[1]] - xr[idx[8]];
    r3 = xr[idx[2]] + xr[idx[7]];
    r4 = xr[idx[2]] - xr[idx[7]];
    r5 = xr[idx[3]] + xr[idx[6]];
    t8 = (xr[idx[3]] - xr[idx[6]]) * C31;
    r7 = xr[idx[4]] + xr[idx[5]];
    r8 = xr[idx[4]] - xr[idx[5]];

    t0 = xr[idx[0]] + r5;
    t7 = xr[idx[0]] + r5 * C95;
    r5 = r1 + r3 + r7;
    xr[idx[0]] = t0 + r5;
    t5 = t0 + r5 * C95;

    t3 = (r3 - r7) * C92;
    r7 = (r1 - r7) * C93;
    r3 = (r1 - r3) * C94;

    t1 = t7 + t3 + r3;
    t3 = t7 - t3 - r7;
    t7 = t7 + r7 - r3;

    t6 = (r2 - r4 + r8) * C31;
    t4 = (r4 + r8) * C96;
    r8 = (r2 - r8) * C97;
    r2 = (r2 + r4) * C98;

    t2 = t8 + t4 + r2;
    t4 = t8 - t4 - r8;
    t8 = t8 + r8 - r2;

    r1 = xi[idx[1]] + xi[idx[8]];
    r2 = xi[idx[1]] - xi[idx[8]];
    r3 = xi[idx[2]] + xi[idx[7]];
    r4 = xi[idx[2]] - xi[idx[7]];
    r5 = xi[idx[3]] + xi[idx[6]];
    r6 = (xi[idx[3]] - xi[idx[6]]) * C31;
    r7 = xi[idx[4]] + xi[idx[5]];
    r8 = xi[idx[4]] - xi[idx[5]];

    t0 = xi[idx[0]] + r5;
    t9 = xi[idx[0]] + r5 * C95;
    r5 = r1 + r3 + r7;
    xi[idx[0]] = t0 + r5;
    r5 = t0 + r5 * C95;

    t0 = (r3 - r7) * C92;
    r7 = (r1 - r7) * C93;
    r3 = (r1 - r3) * C94;

    r1 = t9 + t0 + r3;
    t0 = t9 - t0 - r7;
    r7 = t9 + r7 - r3;

    r9 = (r2 - r4 + r8) * C31;
    r3 = (r4 + r8) * C96;
    r8 = (r2 - r8) * C97;
    r4 = (r2 + r4) * C98;

    r2 = r6 + r3 + r4;
    r3 = r6 - r8 - r3;
    r8 = r6 + r8 - r4;

    xr[ip[1]] = t1 - r2;  
    xr[ip[8]] = t1 + r2;  
    xi[ip[1]] = r1 + t2;
    xi[ip[8]] = r1 - t2;

    xr[ip[2]] = t3 + r3;  
    xr[ip[7]] = t3 - r3;  
    xi[ip[2]] = t0 - t4;
    xi[ip[7]] = t0 + t4;

    xr[ip[3]] = t5 - r9;  
    xr[ip[6]] = t5 + r9;  
    xi[ip[3]] = r5 + t6;
    xi[ip[6]] = r5 - t6;

    xr[ip[4]] = t7 - r8;  
    xr[ip[5]] = t7 + r8;  
    xi[ip[4]] = r7 + t8;
    xi[ip[5]] = r7 - t8;
}

/*!
 *  \brief In-place WFTA module for n = 16.
 *
 *  \param xr   Real part array.
 *  \param xi   Imag part array.
 *  \param idx  Input index map.
 *  \param ip   Output index map.
 */
static void wfta_16(float* xr, float* xi, int * idx, int * ip)
{
    const float C81 = 0.70710678f;    /* 1/sqrt(2) = cos(*M_PI/4) = sin(*M_PI/4) */
    const float C162 = 0.38268343f;   /* sin(M_PI/8) */
    const float C163 = 1.30656297f;   /* (1+sqrt(2)) * sin(M_PI/8) */
    const float C164 = 0.54119610f;   /* sqrt(2)-1 */
    const float C165 = 0.92387953f;   /* cos(M_PI/8) */

    float r1, r2, r3, r4, r5, r6, r7, r8;
    float r9, r10, r11, r12, r13, r14, r15, r16;
    float s2, s3, s4, s6, s7, s8, s9, s10;
    float s11, s12, s13, s14, s15, s16;
    float t1, t2, t3, t4, t5, t6, t7, t8;

    r1 = xr[idx[0]] + xr[idx[8]];
    r2 = xr[idx[0]] - xr[idx[8]];
    r3 = xr[idx[1]] + xr[idx[9]];
    r4 = xr[idx[1]] - xr[idx[9]];
    r5 = xr[idx[2]] + xr[idx[10]];
    r6 = xr[idx[2]] - xr[idx[10]];
    r7 = xr[idx[3]] + xr[idx[11]];
    r8 = xr[idx[3]] - xr[idx[11]];
    r9 = xr[idx[4]] + xr[idx[12]];
    r10 = xr[idx[4]] - xr[idx[12]];
    r11 = xr[idx[5]] + xr[idx[13]];
    r12 = xr[idx[5]] - xr[idx[13]];
    r13 = xr[idx[6]] + xr[idx[14]];
    r14 = xr[idx[6]] - xr[idx[14]];
    r15 = xr[idx[7]] + xr[idx[15]];
    r16 = xr[idx[7]] - xr[idx[15]];

    t1 = r1 + r9;
    t2 = r1 - r9;
    t3 = r3 + r11;
    t4 = r3 - r11;
    t5 = r5 + r13;
    t6 = r5 - r13;
    t7 = r7 + r15;
    t8 = r7 - r15;

    r1 = t1 + t5;
    r3 = t1 - t5;
    r5 = t3 + t7;
    r7 = t3 - t7;

    xr[ip[0]] = r1 + r5;
    xr[ip[8]] = r1 - r5;

    t1 = C81 * (t4 + t8);
    t5 = C81 * (t4 - t8);
    r9 = t2 + t5;
    r11 = t2 - t5;
    r13 = t6 + t1;
    r15 = t6 - t1;

    t1 = r4 + r16;
    t2 = r4 - r16;
    t3 = C81 * (r6 + r14);
    t4 = C81 * (r6 - r14);
    t5 = r8 + r12;
    t6 = r8 - r12;

    t7 = C162 * (t2 - t6);
    t2 = C163 * t2 - t7;
    t6 = C164 * t6 - t7;

    t7 = r2 + t4;
    t8 = r2 - t4;
    r2 = t7 + t2;
    r4 = t7 - t2;
    r6 = t8 + t6;
    r8 = t8 - t6;

    t7 = C165 * (t1 + t5);
    t2 = t7 - C164 * t1;
    t4 = t7 - C163 * t5;

    t6 = r10 + t3;
    t8 = r10 - t3;
    r10 = t6 + t2;
    r12 = t6 - t2;
    r14 = t8 + t4;
    r16 = t8 - t4;

    r1 = xi[idx[0]] + xi[idx[8]];
    s2 = xi[idx[0]] - xi[idx[8]];
    s3 = xi[idx[1]] + xi[idx[9]];
    s4 = xi[idx[1]] - xi[idx[9]];
    r5 = xi[idx[2]] + xi[idx[10]];
    s6 = xi[idx[2]] - xi[idx[10]];
    s7 = xi[idx[3]] + xi[idx[11]];
    s8 = xi[idx[3]] - xi[idx[11]];
    s9 = xi[idx[4]] + xi[idx[12]];
    s10 = xi[idx[4]] - xi[idx[12]];
    s11 = xi[idx[5]] + xi[idx[13]];
    s12 = xi[idx[5]] - xi[idx[13]];
    s13 = xi[idx[6]] + xi[idx[14]];
    s14 = xi[idx[6]] - xi[idx[14]];
    s15 = xi[idx[7]] + xi[idx[15]];
    s16 = xi[idx[7]] - xi[idx[15]];

    t1 = r1 + s9;
    t2 = r1 - s9;
    t3 = s3 + s11;
    t4 = s3 - s11;
    t5 = r5 + s13;
    t6 = r5 - s13;
    t7 = s7 + s15;
    t8 = s7 - s15;

    r1 = t1 + t5;
    s3 = t1 - t5;
    r5 = t3 + t7;
    s7 = t3 - t7;

    xi[ip[0]] = r1 + r5;
    xi[ip[8]] = r1 - r5;

    xr[ip[4]] = r3 + s7;
    xr[ip[12]] = r3 - s7;
    xi[ip[4]] = s3 - r7;
    xi[ip[12]] = s3 + r7;

    t1 = C81 * (t4 + t8);
    t5 = C81 * (t4 - t8);
    s9 = t2 + t5;
    s11 = t2 - t5;
    s13 = t6 + t1;
    s15 = t6 - t1;

    t1 = s4 + s16;
    t2 = s4 - s16;
    t3 = C81 * (s6 + s14);
    t4 = C81 * (s6 - s14);
    t5 = s8 + s12;
    t6 = s8 - s12;

    t7 = C162 * (t2 - t6);
    t2 = C163 * t2 - t7;
    t6 = C164 * t6 - t7;

    t7 = s2 + t4;
    t8 = s2 - t4;
    s2 = t7 + t2;
    s4 = t7 - t2;
    s6 = t8 + t6;
    s8 = t8 - t6;

    t7 = C165 * (t1 + t5);
    t2 = t7 - C164 * t1;
    t4 = t7 - C163 * t5;

    t6 = s10 + t3;
    t8 = s10 - t3;
    s10 = t6 + t2;
    s12 = t6 - t2;
    s14 = t8 + t4;
    s16 = t8 - t4;

    xr[ip[1]] = r2 + s10;
    xr[ip[15]] = r2 - s10;
    xi[ip[1]] = s2 - r10;
    xi[ip[15]] = s2 + r10;

    xr[ip[2]] = r9 + s13;
    xr[ip[14]] = r9 - s13;
    xi[ip[2]] = s9 - r13;
    xi[ip[14]] = s9 + r13;

    xr[ip[3]] = r8 - s16;
    xr[ip[13]] = r8 + s16;
    xi[ip[3]] = s8 + r16;
    xi[ip[13]] = s8 - r16;

    xr[ip[5]] = r6 + s14;
    xr[ip[11]] = r6 - s14;
    xi[ip[5]] = s6 - r14;
    xi[ip[11]] = s6 + r14;

    xr[ip[6]] = r11 - s15;
    xr[ip[10]] = r11 + s15;
    xi[ip[6]] = s11 + r15;
    xi[ip[10]] = s11 - r15;

    xr[ip[7]] = r4 - s12;
    xr[ip[9]] = r4 + s12;
    xi[ip[7]] = s4 + r12;
    xi[ip[9]] = s4 - r12;
}

/**************
 * 
 *  Prime Factor FFT 
 * 
 *   - compute_lp() - computes permutation table
 *   - compute_index_maps() - computes idex map
 *   - pfa_fft() - main PFA FFT program
 */ 

/*!
 *  \brief Compute the LP[] permutation table for one PFA stage.
 *
 *  \param lp   Output array of length n1.
 *  \param n1   Current factor size (e.g., 2,3,4,5,7,8,9,16).
 *  \param n2   N / n1, where N is the full FFT length.
 */
static void compute_lp(int* lp, int n1, int n2)
{
    int L = 0;
    int n3 = n2 - n1 * (n2 / n1);

    for (int j = 0; j < n1; ++j)
    {
        lp[j] = L;
        L += n3;
        if (L >= n1) {
            L -= n1;
        }
    }
}

/*!
 *  \brief Compute index maps idx[] and ip[] for one PFA block.
 *
 *  \param idx  Output array: idx[L] = index of X(I(L))
 *  \param ip   Output array: ip[LP(L)] = same index, permuted.
 *  \param lp   The LP[] permutation table from compute_lp().
 *  \param j0   Starting index of this block (0-based).
 *  \param n1   Current factor size.
 *  \param n2   N / n1.
 *  \param n    Total FFT length.
 */
static void compute_index_maps(int* idx, int* ip, int * lp, int j0, int n1, int n2, int n)
{
    int it = j0;

    for (int L = 0; L < n1; ++L)
    {
        idx[L] = it;
        ip[lp[L]] = it;
        it += n2;
        if (it >= n) {
            it -= n;
        }
    }
}

/*!
 *  \brief Prime Factor FFT (in-place, in-order).
 * 
 *  \param xr   Real part array of length N.
 *  \param xi   Imag part array of length N.
 *  \param n    Total FFT length.
 *  \param m    Number of prime factors.
 *  \param ni   Array of length m containing the factors (ni[0]*...*ni[m-1] = n).
 */
static void fft_pfa(float* xr, float* xi, int n, int m, int * ni)
{
    int idx[16];
    int ip[16];
    int lp[16];

    /* scan composite length factors: */
    for (int k = 0; k < m; ++k) 
    {
        int n1 = ni[k];
        int n2 = n / n1;

        compute_lp(lp, n1, n2);

        /* execute transforms of length n1: */
        for (int j = 0; j < n; j += n1) 
        {
            compute_index_maps(idx, ip, lp, j, n1, n2, n);
            switch (n1) {
                case 2:  wfta_2(xr, xi, idx, ip);  break;
                case 3:  wfta_3(xr, xi, idx, ip);  break;
                case 4:  wfta_4(xr, xi, idx, ip);  break;
                case 5:  wfta_5(xr, xi, idx, ip);  break;
                case 7:  wfta_7(xr, xi, idx, ip);  break;
                case 8:  wfta_8(xr, xi, idx, ip);  break;
                case 9:  wfta_9(xr, xi, idx, ip);  break;
                case 16: wfta_16(xr, xi, idx, ip);  break;
            }
        }
    }
}

/**************
 *
 *  Decimation in Time (DIT) and Decimation in Frequency (DIT) methods. 
 *
 *   - fft_dit() - decimation in time FFT
 *   - fft_dif() - decimation in frequency FFT
 */

/*!
 *  \brief Decimation-in-Frequency, Radix-2 FFT.
 *
 *  \param[in,out]  x - real part of the input/output signal 
 *  \param[in,out]  y - imaginary part of the input/output signal
 *  \param[in]      n - transform size
 *  \param[in]      m - logarithm of transform size
 */
static void fft_dif(float* x, float* y, int n, int m)
{
    int i, j, k, l;
    int n1, n2;
    float e, a, c, s;
    float xt, yt;

    /* main loop (radix-2, DIF process) */
    n2 = n;
    for (k = 0; k < m; ++k) 
    {
        n1 = n2;
        n2 = n2 / 2;
        e = (float)(2.0 * M_PI / (float)n1);
        a = 0.0;
        for (j = 0; j < n2; ++j) 
        {
            /* compute factors: */
            c = cosf(a);
            s = sinf(a);
            a = (j + 1) * e;
            /* butterflies: */
            for (i = j; i < n; i += n1) 
            {
                l = i + n2;
                xt = x[i] - x[l];
                x[i] = x[i] + x[l];
                yt = y[i] - y[l];
                y[i] = y[i] + y[l];
                x[l] = c * xt + s * yt;
                y[l] = c * yt - s * xt;
            }
        }
    }

    /* bit reversal process */
    j = 0;
    n1 = n - 1;
    for (i = 0; i < n1; ++i) 
    {
        if (i < j) {
            /* swap x[i],x[j] */
            xt = x[j];
            x[j] = x[i];
            x[i] = xt;
            /* swap y[i],y[j] */
            xt = y[j];
            y[j] = y[i];
            y[i] = xt;
        }
        k = n / 2;
        while (k <= j && k > 0) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

/*!
 *  \brief Decimation-in-Time, Radix-2 FFT.
 *
 *  \param[in,out]  x - real part of the input/output signal
 *  \param[in,out]  y - imaginary part of the input/output signal
 *  \param[in]      n - transform size
 *  \param[in]      m - logarithm of transform size
 */
static void fft_dit(float* x, float* y, int n, int m)
{
    int i, j, k, l;
    int n1, n2;
    float xt, yt;
    float e, a, c, s;

    /* bit-reversal operation */
    j = 0;  
    n1 = n - 1;
    for (i = 0; i < n1; ++i) 
    { 
        if (i < j) {
            /* swap X(i) <-> X(j) */
            xt = x[j];
            x[j] = x[i];
            x[i] = xt;
            /* swap Y(i) <-> Y(j) */
            xt = y[j];
            y[j] = y[i];
            y[i] = xt;
        }
        k = n >> 1;  /* N/2 */
        while (k <= j && k > 0) {
            j -= k;
            k >>= 1; /* K = K/2 */
        }
        j += k;
    }

    /* main FFT loops - bottom-up construction: */
    n2 = 1;
    for (k = 0; k < m; ++k) 
    {  
        e = (float)(M_PI / (float)n2);  /* 2*pi/(2*N2) = pi/N2 */
        a = 0.0;
        for (j = 0; j < n2; ++j) 
        {  
            /* compute factors: */
            c = cosf(a);
            s = sinf(a);
            a = (j + 1) * e;
            /* compute butteflies */
            for (i = j; i < n; i += 2 * n2) {
                l = i + n2;
                xt = c * x[l] + s * y[l];
                yt = c * y[l] - s * x[l];

                x[l] = x[i] - xt;
                x[i] = x[i] + xt;
                y[l] = y[i] - yt;
                y[i] = y[i] + yt;
            }
        }
        n2 <<= 1;  /* N2 = N2 + N2 */
    }
}

/* dfx_fft.c -- end of file */
