/*  
 *  Copyright 2020 Francesco Gozzini < gozzini AT cpt.univ-mrs.fr >
 *
 *  This file is part of SL2CFOAM-NEXT.
 *
 *  SL2CFOAM-NEXT is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  SL2CFOAM-NEXT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SL2CFOAM-NEXT. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SL2CFOAM_DBG_H__
#define __SL2CFOAM_DBG_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpfr.h>
#include <quadmath.h>
#include <time.h>

#include "integration_qagp.h"
#include "timing.h"
#include "sl2cfoam_tensors.h"

////////////////////////////////////////////////////////////////////////
// Functions useful for debugging.
// The inclusion of this header FORCES DEBUG CODE ON for that file.
////////////////////////////////////////////////////////////////////////

// Useful for writing code only for debugging stage.
#define DEBUG_ON  1
#define DEBUG_OFF 0

#define VARNAME(x) #x

// Print the name and values of variables.
#define DPI(x) printf("DEBUG - "VARNAME(x)": %d\n", (int)x);
#define DPD(x) printf("DEBUG - "VARNAME(x)": %g\n", (double)x);
#define DPC(x) printf("DEBUG - "VARNAME(x)" = %g + I*%g\n", creal((double complex)(x)), cimag((double complex)(x)));
#define DPS(x) printf("DEBUG - "VARNAME(x)": %s\n", (char*)x);
#define DPMPFR(x) printf("DEBUG - "VARNAME(x)": %g\n", (double)mpfr_get_d(x, MPFR_RNDN));
#define DPMPC(x) printf("DEBUG - "VARNAME(x)" = %g + I*%g\n", mpfr_get_d(x->re, MPFR_RNDN), mpfr_get_d(x->im, MPFR_RNDN));
#define DPMPZ(x) printf("DEBUG - "VARNAME(x)" = %ld\n", mpz_get_si(x));

// Prints a matrix of doubles.
#define DPDMAT(m, d1, d2) \
    { \
    printf(VARNAME(m)" = ["); \
    for (int i = 0; i < d1; i++) { \
    for (int j = 0; j < d2; j++) { \
        printf(" %.3g ", matrix_get(m, d1, i, j)); \
    } \
    printf(";\n"); \
    } \
    printf("];\n"); \
    }

// Prints a vector of doubles.
#define DPDVEC(v, d) \
    { \
    printf(VARNAME(v)" = ["); \
    for (int i = 0; i < d; i++) { \
        printf(" %.3g ", vector_get(v, i)); \
    } \
    printf("];\n"); \
    }


// For timing.
#define DPELAPSED(s,...) { printf("DEBUG - "); printf(s, ##__VA_ARGS__); printf(" Time elapsed: %.6f sec.\n", ELAPSED()); }

// Prints a function for plotting in Julia.
static inline void DPFUNC(sl2cfoam_integration_function* f, int N) {

        __float128 dx = 1.0Q / N;
        __float128 x[N-1];

        for (int i = 1; i < N; i++) {
            x[i-1] = i * dx;
        }

        __complex128 d[N-1];
        __complex128 y[N-1];
        __complex128 m[N-1];

        (*(f->function))(y, x, N, f->params);
        (*(f->measure))(m, x, N, f->params);

        for (int i = 1; i < N; i++) {
            d[i-1] = y[i-1]*m[i-1];
        }

        printf("x = [ ");
        for (int i = 0; i < N-2; i++) {
            printf("%.4f, ", (double)(x[i]));
        }
        printf("%.4f ];\n", (double)(x[N-2]));

        printf("d = [ ");
        double dr, di;
        for (int i = 0; i < N-2; i++) {
            dr = creal((double complex)d[i]);
            di = cimag((double complex)d[i]);
            printf("%.6g + %.6gim, ", dr, di);
        }
        dr = creal((double complex)d[N-2]);
        di = cimag((double complex)d[N-2]);
        printf("%.6g + %.6gim ];\n", dr, di);

}

// For checking memory.
#ifdef USE_MKL
#include <mkl.h>
#define MKL_MEM_ENABLE() { mkl_peak_mem_usage(MKL_PEAK_MEM_ENABLE); }
#define MKL_MEM_PEAK() { printf("DEBUG - "); printf("MKL peak memory allocation: %lld bytes.\n", mkl_peak_mem_usage(MKL_PEAK_MEM_RESET)); }
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_DBG_H__*/