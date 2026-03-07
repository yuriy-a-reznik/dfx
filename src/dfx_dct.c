/*!
 *  \file   dfx_dct.c
 *  \brief  DFX library: implementations of DCT-II and DCT-IV transforms of various sizes.
 *
 *  This module includes examples of:
 * 
 *     Short length DCT-II modules (n1 = 2,3,4,5,7,9)
 *     PFA method (n1 = 3*5 =15)
 *     Short length DCT-IV models 
 *     Recursive construction of DCT-II of sizes n=2^m * n1, using split in DCT-II and DCT-IV of half-lengths.
 *
 *  Key references: 
 *    DCT via DFT
 *      N. Ahmed, T. Natarajan, K. R.; Rao, 1974, "Discrete Cosine Transform," IEEE Transactions on Computers. C-23 (1): 90-93.
 *      R. Haralick, "A Storage Efficient Way to Implement the Discrete Cosine Transform," IEEE Transactions on Computers, vol. C-25, no. 7, pp. 764-765, 1976
 *      M. T. Heideman, "Computation of an odd-length DCT from a real valued DFT of the same length," TSP, vol. 40,no. 1, pp. 54-61, 1992.
 *    Prime factor method
 *      C.S.Burrus and T.W.Parks, "DFT/FFT and convolution algorithms: theory and implementation," John Wiley & Sons, Inc., 1991.
 *      P. Yang and M. Narasimha, "Prime factor decomposition of the discrete cosine transform," ICASSP, 1985.   
 *      B. G. Lee, "Input and output index mapping for a prime-factor decomposed computation of discrete cosine transform," TASSP, vol. 37, no. 2, 1989.
 *    Connection between DCT-II and DCT-IV
 *      C.W. Kok, "Fast Algorithm for Computing Discrete Cosine Transform", IEEE Trans. Signal Proc., vol.45, no.3, pp.757-760, Mar 1997
 *    Dyadic DCT-II and DCT-VI factorizations
 *      W.-H. Chen, C. Smith and S. Fralick, "A Fast Computational Algorithm for the Discrete Cosine Transform," IEEE Trans. Communications, vol. 25, no. 9, pp. 1004-1009, 1977
 *      Z. Wang, On Computing the Discrete Fourier and Cosine Transforms, IEEE Trans. ASSP, Vol 33, No. 4, Oct. 1985.
 *      G.Plonka and M.Tashe, "Fast and numerically stable algorithms for discrete cosine transforms", Linear Algebra and Applications, 394 (1) 309-345 (2005). 
 *      V. Britanak, P. Yip, K. R. Rao, "Discrete Cosine and Sine Transforms: General Properties, Fast Algorithms and Integer Approximations", Academic Press, 2006.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dfx.h"

/*************
 * 
 *  Short-length, in-place, DCT-II modules for n1 = 2, 3, 4, 5, 7, 9.
 * 
 *  Odd-size transforms are derived from WFTA modules as described in:
 *    M. T. Heideman, "Computation of an odd-length DCT from a real valued DFT of the same length," TSP, vol. 40,no. 1, pp. 54-61, 1992.
 */

/*!
 *  \brief Length-2, DCT-II (1 multiplication, 2 adds)
 */
static void dct2_2(float* x)
{ 
    const float C21 = 0.707106781f;  /* cos(pi/4) = sqrt(2)/2 */
    float x0, x1;

    x0 = x[0];
    x1 = x[1];

    x[0] = x0 + x1; 
    x[1] = (x0 - x1) * C21; 
}

/*!
 *  \brief Length-3 DCT-II (1 multiplication, 4 adds, 1 shift)
 */
static void dct2_3(float* x)
{
    const float C31 = -0.866025404f;  /* -cos(pi/6) = -sqrt(3)/2 */
    const float C32 = 1.5f;

    float a1, a2, m1, m2;

    a1 = x[0] + x[2];
    a2 = x[2] - x[0];

    m1 = C31 * a2;
    m2 = C32 * a1;

    x[0] = x[1] + a1;
    x[1] = m1;
    x[2] = m2 - x[0];
}

/*!
 *  \brief Length-4 DCT-II (5 multiplications, 6 adds)
 */
static void dct2_4(float* x)
{
    const float C21 = 0.707106781f;  /* cos(pi/4) = sqrt(2)/2 */
    const float C41 = 0.923879532f;  /* cos(pi/8) */
    const float C42 = 0.382683432f;  /* sin(pi/8) */

    float a1, a2, a3, a4;

    a1 = x[0] + x[3];
    a2 = x[1] + x[2];
    a3 = x[0] - x[3];
    a4 = x[1] - x[2];

    x[0] = a1 + a2;
    x[2] = (a1 - a2) * C21;   
    x[1] = C41 * a3 + C42 * a4;
    x[3] = C42 * a3 - C41 * a4;
}

/*!
 *  \brief 5-point DCT-II (5 multiplications, 13 adds)
 */
static void dct2_5(float* x)
{
#if 0
    /* algorithm from Heideman's paper: */
    const float C51 = 0.95105651629515357f;    /* (-sqrt(10.0) + 2.0 * sqrt(5.0)) / 4.0 */
    const float C52 = 1.53884179671731191f;    /* (sqrt(5.0) + 2.0 * sqrt(5.0)) / 2.0   */
    const float C53 = -0.36327123629798958f;   /* (-sqrt(5.0) - 2.0 * sqrt(5.0)) / 2.0  */
    const float C54 = -0.55901699437494742f;   /* -sqrt(5.0) / 4.0 */
    const float C55 = -1.25;

    float a1, a2, a3, a4, a5, a6, a7, a8, m1, m2, m3, m4, m5, c0;

    a1 = x[1] + x[3];
    a2 = x[1] - x[3];
    a3 = x[4] + x[0];
    a4 = x[4] - x[0];

    a5 = a1 + a3;
    a6 = a1 - a3;
    a7 = a2 + a4;

    m1 = C54 * a6;
    m2 = C55 * a5;
    m3 = C51 * a7;
    m4 = C52 * a2;
    m5 = C53 * a4;

    c0 = a5 + x[2];
    a8 = c0 + m2;

    x[0] = c0;
    x[1] = m3 + m4;
    x[2] = m1 - a8;
    x[3] = m3 - m5;
    x[4] = m1 + a8;
#else
    /* slightly better factorization (Chivukula-Reznik, 2007) */
    float t0, t1, t3, t4, ta, tb;
    ta = x[0] + x[4];
    t4 = x[4] - x[0];
    t3 = x[3] - x[1];
    tb = x[3] + x[1];
    t0 = ta + tb;
    t1 = ta - tb;
    ta = -0.58778525229247313f * t4;		/* cos(3Pi/10) */
    tb = -0.58778525229247313f * t3;
    x[3] = ta + 0.95105651629515353f * t3;	/* cos(Pi/10) */
    x[1] = tb - 0.95105651629515353f * t4;
    ta = x[2] - 0.25f * t0;
    x[0] = x[2] + t0;
    tb = 0.55901699437494745f * t1;         /* cos(pi/5) - 1/4 */
    x[4] = ta + tb;
    x[2] = -ta + tb;
#endif
}

/*
 *  \brief 7-point DCT-II (8 multiplications, 30 adds)
 */
static void dct2_7(float* x)
{
    const float C71 = -1.1666666666666667f;  /* -7/6 */
    const float C72 = -0.7901560896121518f;  /* -(2cos u - cos2u - cos3u)/3, u=2pi/7 */
    const float C73 = -0.0558537874849908f;  /* -(cos u - 2cos2u + cos3u)/3, u=2pi/7 */
    const float C74 = 0.7343021342664287f;   /*  (cos u + cos2u - 2cos3u)/3, u=2pi/7 */
    const float C75 = 0.4409591763108903f;   /*  (sin u + sin2u - sin3u)/3, u=2pi/7  */
    const float C76 = -0.3408730336181886f;  /* -(2sin u - sin2u + sin3u)/3, u=2pi/7 */
    const float C77 = -0.5339694823505409f;  /*  (sin u - 2sin2u - sin3u)/3, u=2pi/7 */
    const float C78 = 0.8748419904388223f;   /*  (sin u + sin2u + 2sin3u)/3, u=2pi/7 */

    float a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16;
    float m1, m2, m3, m4, m5, m6, m7, m8, c0, a17, a18, a19, a20, a21, a22, a23;

    a1 = x[4] + x[2];
    a2 = x[4] - x[2];
    a3 = x[1] + x[5];
    a4 = x[1] - x[5];
    a5 = x[6] + x[0];
    a6 = x[6] - x[0];

    a7 = a1 - a3;
    a8 = a1 - a5;
    a9 = a1 + a3;
    a10 = a9 + a5;
    a11 = a5 - a3;

    a12 = a2 - a4;
    a13 = a2 + a6;
    a14 = a2 + a4;
    a15 = a14 - a6;
    a16 = a4 + a6;

    m1 = C74 * a7;
    m2 = C72 * a8;
    m3 = C71 * a10;
    m4 = C73 * a11;
    m5 = C78 * a12;
    m6 = C76 * a13;
    m7 = C75 * a15;
    m8 = C77 * a16;

    c0 = x[3] + a10;
    a17 = c0 + m3;
    a18 = a17 - m2;
    a19 = a17 + m4;
    a20 = a17 + m2;
    a21 = m7 - m6;
    a22 = m7 + m8;
    a23 = m7 + m6;

    x[0] = c0;
    x[1] = a22 - m5;
    x[2] = m4 - a18;
    x[3] = a23 + m5;
    x[4] = a20 + m1;
    x[5] = m8 - a21;
    x[6] = m1 - a19;
}

/*
 *  \brief 9-point DCT-II (10 multiplications, 34 adds)
 */
static void dct2_9(float* x)
{
    const float C31 = -0.8660254037844386f; /* -sqrt(3)/2        */
    const float C92 = 0.9396926207859084f;  /* cos(4u), u=2pi/9 */
    const float C93 = -0.1736481776746724f; /* -cos(2u), u=2pi/9 */
    const float C94 = -0.7660444431189780f; /* -cos(u), u=2pi/9  */
    const float C95 = 0.5f;
    const float C96 = -0.3420201433256687f; /* -sin(4u), u=2pi/9 */
    const float C97 = -0.9848077530122080f; /* -sin(2u), u=2pi/9 */
    const float C98 = -0.6427876096865393f; /* -sin(u), u=2pi/9  */

    float a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19;
    float m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, a20, a21, a22, a23, a24, a25, a26;

    a1 = x[3] + x[5];
    a2 = x[3] - x[5];
    a3 = x[6] + x[2];
    a4 = x[6] - x[2];
    a5 = x[1] + x[7];
    a6 = x[1] - x[7];
    a7 = x[8] + x[0];
    a8 = x[8] - x[0];
    a9 = x[4] + a5;

    a10 = a1 + a3;
    a11 = a10 + a7;
    a12 = a3 - a7;
    a13 = a1 - a7;
    a14 = a1 - a3;
    a15 = a2 - a4;
    a16 = a15 + a8;
    a17 = a4 + a8;
    a18 = a2 - a8;
    a19 = a2 + a4;

    m1 = C31 * a6;
    m2 = C95 * a5;
    m3 = C95 * a11;
    m4 = C92 * a12;
    m5 = C93 * a13;
    m6 = C94 * a14;
    m7 = C31 * a16;
    m8 = C96 * a17;
    m9 = C97 * a18;
    m10 = C98 * a19;

    a20 = x[4] - m2;
    a21 = a20 + m4;
    a22 = a20 - m4;
    a23 = a20 + m5;
    a24 = m1 + m8;
    a25 = m1 - m8;
    a26 = m1 + m9;

    x[0] = a9 + a11;
    x[1] = m10 - a26;
    x[2] = m6 - a21;
    x[3] = m7;
    x[4] = a22 - m5;
    x[5] = a25 - m9;
    x[6] = m3 - a9;
    x[7] = a24 + m10;
    x[8] = a23 + m6;
}

/*
 * Fast 15-point DCT-II (14 multiplications, 67 additions, 3 multiplications by dyadic rationals):
 *
 *  15-point (3x5) DCT-II factorization constructed by using PFA technique. 
 *  Full construction is described in:
 *    R. K. Chivukula, Fast Algorithms For MDCT And Low Delay Filterbanks Used In Audio Coding, MS Thesis, UT Arlington, 2008
 *    https://mavmatrix.uta.edu/context/electricaleng_theses/article/1147/type/native/viewcontent
 */
static void dct2_15(float* x)
{
    const float C1 = -1.25f;                /* 1/2*(cos(a) + cos(2a)), a=-2*Pi/5 => -1.25 */
    const float C2 = 0.559016994374947f;    /* 1/2*(cos(a) - cos(2a)), a=-2*Pi/5 */
    const float C3 = -1.53884176858763f;    /* sin(a) + sin(2a), a=-2*Pi/5 */
    const float C4 = 0.587785252292473f;    /* sin(2a), a=-2*Pi/5 */
    const float C5 = -0.36327126400268f;    /* sin(a) - sin(2a), a=-2*Pi/5 */
    const float C6 = -1.5f;                 /* cos(b)-1, b=-2*Pi/3 => -3/2 */
    const float C7 = 1.875f;                /* C1 * C6 => 15/8 */
    const float C8 = -0.838525491562421f;   /* C2 * C6 */
    const float C9 = 2.30826265288144f;     /* C3 * C6 */
    const float C10 = 0.88167787843871f;    /* C4 * C2 */
    const float C11 = 0.54490689600402f;    /* C5 * C6 */
    const float C12 = -0.866025403784439f;  /* sin(b), b=-2*Pi/3 */
    const float C13 = 1.08253175473055f;    /* C1 * C12 */
    const float C14 = -0.484122918275927f;  /* C2 * C12 */
    const float C15 = -1.33267606400146f;   /* -C3 * C12 */
    const float C16 = -0.509036960455127f;  /* C4 * C12  */
    const float C17 = -0.314602143091205f;  /* -C5 * C12 */

    float in[15], m[6][3], e[8], o[7];
    float s1, s2, s3, s4, s5;

    in[12] = x[12] + x[2];
    in[2] = x[12] - x[2];
    in[7] = x[7] + in[12];

    in[0] = x[0] + x[9];
    in[9] = x[0] - x[9];
    in[10] = x[10] + in[0];

    in[11] = x[11] + x[8];
    in[8] = x[11] - x[8];
    in[1] = x[1] + in[11];

    in[6] = x[6] + x[3];
    in[3] = x[6] - x[3];
    in[13] = x[13] + in[6];
    
    in[5] = x[5] + x[14];
    in[14] = x[5] - x[14];
    in[4] = x[4] + in[5];

    s1 = in[10] + in[4];
    s2 = in[10] - in[4];
    s3 = in[13] + in[1];
    s4 = in[13] - in[1];
    s5 = s1 + s3;
    m[0][0] = s5 + in[7];
    m[1][0] = s5 * C1;
    m[2][0] = (s1 - s3) * C2;
    m[3][0] = s2 * C3;;
    m[4][0] = -(s4 + s2) * C4;
    m[5][0] = s4 * C5;

    s1 = in[0] + in[5];
    s2 = in[0] - in[5];
    s3 = in[6] + in[11];
    s4 = in[6] - in[11];
    s5 = s1 + s3;
    m[0][1] = (s5 + in[12]) * C6;
    m[1][1] = s5 * C7;
    m[2][1] = (s1 - s3) * C8;
    m[3][1] = s2 * C9;
    m[4][1] = (s4 + s2) * C10;
    m[5][1] = s4 * C11;

    s1 = in[9] + in[14];
    s2 = in[9] - in[14];
    s3 = in[3] + in[8];
    s4 = in[3] - in[8];
    s5 = s1 + s3;
    m[0][2] = (s5 + in[2]) * C12;
    m[1][2] = s5 * C13;
    m[2][2] = (s1 - s3) * C14;
    m[3][2] = s2 * C15;
    m[4][2] = (s4 + s2) * C16;
    m[5][2] = s4 * C17;

    s1 = m[0][0] + m[1][0];
    e[0] = m[0][0];
    e[1] = s1 + m[2][0];
    e[2] = s1 - m[2][0];
    o[0] = m[3][0] - m[4][0];
    o[1] = m[4][0] + m[5][0];
    s1 = m[0][1] + m[1][1];
    e[3] = m[0][1];
    e[4] = s1 + m[2][1];
    e[5] = s1 - m[2][1];
    o[2] = m[3][1] - m[4][1];
    o[3] = m[4][1] + m[5][1];
    s1 = m[0][2] + m[1][2];
    o[4] = m[0][2];
    o[5] = s1 + m[2][2];
    o[6] = s1 - m[2][2];
    e[6] = m[3][2] - m[4][2];
    e[7] = m[4][2] + m[5][2];

    e[3] += e[0];
    e[4] += e[1];
    o[2] += o[0];
    x[2] = -(e[4] + e[6]);
    x[8] = e[4] - e[6];
    x[13] = o[2] + o[5];
    x[7] = o[2] - o[5];

    e[5] += e[2];
    o[3] += o[1];
    x[14] = -(e[5] + e[7]);
    x[4] = e[5] - e[7];
    x[1] = o[3] + o[6];
    x[11] = -o[3] + o[6];

    x[0] = e[0];
    x[3] = -o[0];
    x[6] = -e[2];
    x[10] = -e[3];
    x[12] = e[1];
    x[9] = -o[1];
    x[5] = -o[4];
}

/*****************
 *
 *  Short-length, in-place, DCT-IV modules for n1 = 2, 3, 4, 5, 7, 9 and 15.
 *  
 *  In most cases, we construct DCT-IV by using DCT-II and C.W.Kok's conversion method:
 *    C.W.Kok, "Fast Algorithm for Computing Discrete Cosine Transform", IEEE Trans. Signal Proc., vol.45, no.3, pp.757-760, Mar 1997
 *
 *  However, this is not the fastest method available!
 *  For several sizes, more efficient short-length DCT-IVs have been recently reported in
 *    A. Cariow, L. Lesiecki, Small-Size Algorithms for Type-IV Discrete Cosine Transform with Reduced Multiplicative Complexity,
 *    Radioelectron.Commun.Syst. 63, 465-487 (2020).
 */

/*!
 *  \brief Length 2 DCT-VI 
 *
 *   Direct factorization with 2 multiplications and 2 adds
 */
static void dct4_2(float* x)
{
    const float C1 = 0.9238795325112867f; /* cos(pi/8) */ 
    const float C3 = 0.3826834323650898f; /* cos(3pi/8) */ 
    
    float x0 = x[0]; 
    float x1 = x[1]; 

    x[0] = C1 * x0 + C3 * x1; 
    x[1] = C3 * x0 - C1 * x1; 
}

/*!
 *  \brief Length 3 DCT-VI 
 * 
 *   Direct factorization with 4 multiplications and 6 adds
 */
static void dct4_3(float* x)
{
    const float C31 = 0.9659258262890683f; /* cos(pi/12) */ 
    const float C32 = 0.7071067811865475f; /* cos(pi/4) */ 
    const float C33 = 0.2588190451025207f; /* cos(5*pi/12) */ 
    const float C34 = 0.5f * (C31 + C33); /* (C31 + C33)/2 */ 
    const float C35 = 0.5f * (C31 - C33); /* (C31 - C33)/2 */ 
    
    float s, d, t1, t2, t3, x1, bx1;
    
    s = x[0] + x[2]; 
    d = x[0] - x[2]; 
    t1 = d - x[1]; 
    x1 = C32 * t1;
    t2 = C34 * s; 
    t3 = C35 * d; 
    bx1 = t3 + C32 * x[1]; 
    x[0] = t2 + bx1;
    x[1] = x1; 
    x[2] = t2 - bx1;
};

/*!
 *  \brief Length 4 DCT-IV:
 */
static void dct4_4(float* x)
{
    x[0] *= 1.9615705608065f;  /* 2*cos(M_PI*1/16); */
    x[1] *= 1.6629392246051f;  /* 2*cos(M_PI*3/16); */
    x[2] *= 1.1111404660392f;  /* 2*cos(M_PI*5/16); */
    x[3] *= .39018064403232f;  /* 2*cos(M_PI*7/16); */

    dct2_4(x);

    x[0] *= 0.5;
    x[1] -= x[0];
    x[2] -= x[1];
    x[3] -= x[2];
}

/*!
 *  \brief Length-5 DCT-IV
 */
static void dct4_5(float* x)
{
    x[0] *= 1.9753766811903f;    /* 2*cos(pi*1/20); */
    x[1] *= 1.7820130483767f;    /* 2*cos(pi*3/20); */
    x[2] *= 1.4142135623731f;    /* 2*cos(pi*5/20); */
    x[3] *= 0.90798099947914f;   /* 2*cos(pi*7/20); */
    x[4] *= 0.31286893008048f;   /* 2*cos(pi*9/20); */

    dct2_5(x);

    x[0] *= 0.5f;
    x[1] -= x[0];
    x[2] -= x[1];
    x[3] -= x[2];
    x[4] -= x[3];
}

/*!
 *  \brief Length-7 DCT-IV
 */
static void dct4_7(float* x)
{
    x[0] *= 1.98742441979f;      /* 2*cos(pi*1,28); */
    x[1] *= 1.88776666062f;      /* 2*cos(pi*3,28); */
    x[2] *= 1.69344839846f;      /* 2*cos(pi*5,28); */
    x[3] *= 1.41421356237f;      /* 2*cos(pi*7,28); */
    x[4] *= 1.06406415302f;      /* 2*cos(pi*9,28); */
    x[5] *= 0.660558123910f;     /* 2*cos(pi*11,28); */
    x[6] *= 0.223928952200f;     /* 2*cos(pi*13,28); */

    dct2_7(x);

    x[0] *= 0.5f;
    x[1] -= x[0];
    x[2] -= x[1];
    x[3] -= x[2];
    x[4] -= x[3];
    x[5] -= x[4];
    x[6] -= x[5];
}

/*!
 *  \brief Length-9 DCT-IV
 */
static void dct4_9(float* x)
{
    x[0] *= 1.99238939618f;     /* 2*cos(pi*1/36); */
    x[1] *= 1.93185165258f;     /* 2*cos(pi*3/36); */
    x[2] *= 1.81261557407f;     /* 2*cos(pi*5/36); */
    x[3] *= 1.63830408858f;     /* 2*cos(pi*7/36); */
    x[4] *= 1.41421356237f;     /* 2*cos(pi*9/36); */
    x[5] *= 1.14715287270f;     /* 2*cos(pi*11/36); */
    x[6] *= 0.845236523474f;    /* 2*cos(pi*13/36); */
    x[7] *= 0.517638090196f;    /* 2*cos(pi*15/36); */
    x[8] *= 0.174311485506f;    /* 2*cos(pi*17/36); */
    
    dct2_9(x);

    x[0] *= 0.5f;
    x[1] -= x[0];
    x[2] -= x[1];
    x[3] -= x[2];
    x[4] -= x[3];
    x[5] -= x[4];
    x[6] -= x[5];
    x[7] -= x[6];
    x[8] -= x[7];
}

/*!
 *  \brief Length-15, DCT-IV
 */
static void dct4_15(float* x)
{
    x[0] *= 1.9972590695091f;    /* 2*cos(pi*1/60); */
    x[1] *= 1.9753766811903f;    /* 2*cos(pi*3/60); */
    x[2] *= 1.9318516525781f;    /* 2*cos(pi*5/60); */
    x[3] *= 1.8671608529944f;    /* 2*cos(pi*7/60); */
    x[4] *= 1.7820130483767f;    /* 2*cos(pi*9/60); */
    x[5] *= 1.6773411358909f;    /* 2*cos(pi*11/60); */
    x[6] *= 1.5542919229139f;    /* 2*cos(pi*13/60); */
    x[7] *= 1.4142135623731f;    /* 2*cos(pi*15/60); */
    x[8] *= 1.2586407820997f;    /* 2*cos(pi*17/60); */
    x[9] *= 1.0892780700300f;    /* 2*cos(pi*19/60); */
    x[10] *= 0.90798099947914f;  /* 2*cos(pi*21/60); */
    x[11] *= 0.71673589909058f;  /* 2*cos(pi*23/60); */
    x[12] *= 0.51763809020494f;  /* 2*cos(pi*25/60); */
    x[13] *= 0.31286893008048f;  /* 2*cos(pi*27/60); */
    x[14] *= 0.10467191248582f;  /* 2*cos(pi*29/60); */

    dct2_15(x);

    x[0] *= 0.5f;
    x[1] -= x[0];
    x[2] -= x[1];
    x[3] -= x[2];
    x[4] -= x[3];
    x[5] -= x[4];
    x[6] -= x[5];
    x[7] -= x[6];
    x[8] -= x[7];
    x[9] -= x[8];
    x[10] -= x[9];
    x[11] -= x[10];
    x[12] -= x[11];
    x[13] -= x[12];
    x[14] -= x[13];
}


/*****************
 *
 *  Recursive factorization for DCT-II and DCT-IV transforms of lengths n = 2^m * n1, where n1 = 3,4,5,7,9,15.
 *
 *    dct4_rec(float* x, float* y, int n, int n_min);
 *    dct2_rec(float* x, float* y, int n, int n_min);
 *
 *  This recursion is similar to the one described in:
 *      G.Plonka and M.Tashe, "Fast and numerically stable algorithms for discrete cosine transforms",
 *      Linear Algebra and Applications, 394 (1) 309-345 (2005).
 *
 *  This version is generalized to support non-dyadic transform lengths. The resulting transforms are unnormalized.  
 *  This factorization is numerically stable, and allows computation of large-size transforms.
 */

int dct2_rec(float *x, float *y, int n);
int dct4_rec(float *x, float *y, int n);

/*!
 *  \brief Recursive DCT-II of lengths n = 2 ^ m * n1, where n1 = 3, 4, 5, 7, 9, 15.
 *
 *  \param[in, out]  x - vector to be transformed
 *  \param[out]      y - temp buffer
 *  \param[in]       n - transform size
 *
 *  \returns 0 - success, 1 - unsupported transform size
 */
int dct2_rec(float *x, float* y, int n)
{

    /* check if split is possible: */
    if (!(n & 1) && n > 4)
    {
        int n1 = n / 2, i, err = 0;

        /* input butterfly: */
        for (i = 0; i < n1; i++) {
            y[i] = x[i] + x[n - 1 - i];
            y[i + n1] = x[i] - x[n - 1 - i];
        }

        /* Chen-Smith-Fralick split in DCT-II and DCT-IV: */
        err += dct2_rec(y, x, n1);
        err += dct4_rec(&y[n1], x, n1);
        if (err) return err;

        /* re-order and store the results: */
        for (i = 0; i < n1; i++)
            x[2 * i] = y[i];
        for (i = 0; i < n1; i++)
            x[2 * i + 1] = y[n1 + i];
    }
    /* check if we can use one of our small-size transforms: */
    else if (n == 15)
    {
        dct2_15(x);
    }
    else if (n == 9)
    {
        dct2_9(x);
    }
    else if (n == 7)
    {
        dct2_7(x);
    }
    else if (n == 5)
    {
        dct2_5(x);
    }
    else if (n == 4)
    {
        dct2_4(x);
    }
    else if (n == 3)
    {
        dct2_3(x);
    }
    else
    {
        /* unsupported transform size */
        return 1;
    }

    /* success:*/
    return 0;
}

 /*!
  *  \brief Recursive DCT-IV of lengths n = 2^m * n1, where n1 = 3,4,5,7,9,15.
  *
  *  \param[in,out]  x - vector to be transformed
  *  \param[out]     y - temp buffer
  *  \param[in]      n - transform size
  * 
  *  \returns 0 - success, !0 - unsupported transform size
  */
int dct4_rec(float* x, float* y, int n)
{
    /* check if split is possible: */
    if (!(n & 1) && n > 4)
    {
        int n1 = n / 2, i, err = 0;
        float temp;

        /* Givens' rotations: */
        for (i = 0; i < n1; i++) {
            /* with proper implementation, all these factors should be pre-computed */
            y[i] = (float)(x[i] * cos(M_PI * (double)(2 * i + 1) / (double)(4 * n)) + x[n - 1 - i] * sin(M_PI * (double)(2 * i + 1)/ (double)(4 * n)));
            y[n - 1 - i] = (float)(-x[i] * sin(M_PI * (double)(2 * i + 1) / (double)(4 * n)) + x[n - 1 - i] * cos(M_PI * (double)(2 * i + 1)/ (double) (4 * n)));
        }
        /* sign inversions */
        for (i = n1 + 1; i < n; i += 2)
            y[i] = -y[i];

        /* compute half-length DCT-IIs: */
        err += dct2_rec(y, x, n1);
        err += dct2_rec(&y[n1], x, n1);
        if (err) return err;

        /* combination of resuls: */
        y += n1;
        for (i = 0; i < n1 / 2; i++) {
            temp = y[i];
            y[i] = y[n1 - 1 - i];
            y[n1 - 1 - i] = temp;
        }
        y -= n1;

        /* sign inversions */
        for (i = n1 + 1; i < n; i += 2)
            y[i] = -y[i];

        /* final buttefly: */
        y[n - 1] = -y[n - 1];
        for (i = 1; i < n1; i++) {
            temp = y[i] + y[n1 - 1 + i];
            y[n1 - 1 + i] = y[i] - y[n1 - 1 + i];
            y[i] = temp;
        }
        /* store outputs: */
        for (i = 0; i < n1; i++)
            x[2 * i] = y[i];
        for (i = 0; i < n1; i++)
            x[2 * i + 1] = y[n1 + i];
    }
    /* check if we can use one of our small-size transforms: */
    else if (n == 15)
    {
        dct4_15(x);
    }
    else if (n == 9)
    {
        dct4_9(x);
    }
    else if (n == 7)
    {
        dct4_7(x);
    }
    else if (n == 5)
    {
        dct4_5(x);
    }
    else if (n == 4)
    {
        dct4_4(x);
    }
    else if (n == 3)
    {
        dct4_3(x);
    }
    else
    {
        /* unsupported transform size */
        return 1;
    }

    /* success:*/
    return 0;
}


#if 0

/**************** 
 * 
 *  Test program and demo
 * 
 *   This is a simple unit test for all above defined functions. 
 * 
 *   It includes reference implementations of DCT-II and DCT-IV, random signal generator, 
 *   and main program that executes and compares the outputs of the above modules agains the reference.
 *   It reports absolute average errors as observed for each transform. 
 * 
 */

#define N_MAX 1024

/*! 
 *  \brief Reference unnormalized DCT-II:
 *
 *   X[k] = sum_{n=0}^{n-1} x[n] * cos( pi*(n+0.5)*k / N), k = 0..N-1
 */
static void dct2_ref(float *x, float *X, int N)
{
    int n, k;
    /* compute DCT-II */
    for (n = 0; n < N; n++) {
        double s = 0.;
        for (k = 0; k < N; k++)
            s += x[k] * cos(M_PI * (double)(n * (2 * k + 1)) / (double)(2 * N));
        X[n] = (float)s;
    }
}

/*!
 *  \brief Reference unnormalized DCT-IV:
 *
 *   X[k] = sum_{n=0}^{n-1} x[n] * cos( pi*(n + 0.5)*(k + 0.5) / N), k = 0..N-1
 */
static void dct4_ref(float *x, float *X, int N)
{
    int n, k;
    for (k = 0; k < N; k++) {
        double s = 0.;
        for (n = 0; n < N; n++)
            s += x[n] * cos(M_PI * (double)(2 * n + 1) * (2 * k + 1)/(double) (4 * N));
        X[k] = (float)s;
    }
}

/*!
 *  \brief Generate random vector for testing
 */
static void get_test_vector(float* x, int n)
{
    int k;
    for (k = 0; k < n; k++) {
        x[k] = (float)rand() / (float)RAND_MAX;
    }
}

/*!
 *  \brief Compute total absolute error between 2 vectors
 */
static double abs_err(float *x, float *y, int n)
{
    int k; double s = 0;
    for (k = 0; k < n; k++) {
        s += fabs(x[k] - y[k]);
    }
    return s;
}

/*! 
 *  Test/demo program: 
 */
int main()
{
    /* test data arrays */
    static float x[N_MAX];
    static float X[N_MAX];
    static float X_ref[N_MAX];
    static float Y[N_MAX];

    /* transform sizes & functions: */
    static void (*dct2_test[])(float*) = {dct2_2, dct2_3, dct2_4, dct2_5, dct2_7, dct2_9, dct2_15};
    static void (*dct4_test[])(float*) = {dct4_2, dct4_3, dct4_4, dct4_5, dct4_7, dct4_9, dct4_15};
    static int N[] = {2, 3, 4, 5, 7, 9, 15}; 
    static int M = sizeof(N) / sizeof(N[0]);
    static int N2[] = {192, 256, 320, 448, 576, 960};
    static int M2 = sizeof(N2) / sizeof(N2[0]);
    static int Q = 10;

    /* local variables: */
    int m, n, q; double s;

    /* test short-lenth transforms: */
    printf("Short-length transfom modules:\n");
    for (m = 0; m < M; m++)
    {
        /* pick DCT length */
        n = N[m];

        /* test DCT-II: */
        printf("Testing %s, n=%d: ", "dct2", n);
        for (q = 0, s = 0.; q < Q; q++) {
            get_test_vector(x, n);
            dct2_ref(x, X_ref, n);
            memcpy(X, x, n * sizeof(float));
            dct2_test[m](X);
            s += abs_err(X, X_ref, n);
        }
        printf("err=%g\n", s / (double)(n*Q));

        /* test DCT-IV: */
        printf("Testing %s, n=%d: ", "dct4", n);
        for (q = 0, s = 0.; q < Q; q++) {
            get_test_vector(x, n);
            dct4_ref(x, X_ref, n);
            memcpy(X, x, n * sizeof(float));
            dct4_test[m](X);
            s += abs_err(X, X_ref, n);
        }
        printf("err=%g\n", s / (double)(n * Q));
    }

    /* test composite-lenth transforms: */
    printf("Composite-length transfoms:\n");
    for (m = 0; m < M2; m++)
    {
        /* pick DCT length */
        n = N2[m];

        /* test DCT-II: */
        printf("Testing %s, n=%d: ", "dct2", n);
        for (q = 0, s = 0.; q < Q; q++) {
            get_test_vector(x, n);
            dct2_ref(x, X_ref, n);
            memcpy(X, x, n * sizeof(float));
            dct2_rec(X, Y, n);
            s += abs_err(X, X_ref, n);
        }
        printf("err=%g\n", s / (double)(n * Q));

        /* test DCT-IV: */
        printf("Testing %s, n=%d: ", "dct4", n);
        for (q = 0, s = 0.; q < Q; q++) {
            get_test_vector(x, n);
            dct4_ref(x, X_ref, n);
            memcpy(X, x, n * sizeof(float));
            dct4_rec(X, Y, n);
            s += abs_err(X, X_ref, n);
        }
        printf("err=%g\n", s / (double)(n * Q));
    }

    /* all done */
    return 0;
}

#endif

/* dfx_dct.c -- end of file */