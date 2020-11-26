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

#ifndef __SL2CFOAM_JSYMBOLS_H__
#define __SL2CFOAM_JSYMBOLS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"
#include "jsymbols.h"
#include "utils.h"
#include "wigxjpf.h"
#include "fastwigxj.h"
#include "sl2cfoam.h"
#include "error.h"

///////////////////////////////////////////////////////////
// Functions for computing wigner symbols.
///////////////////////////////////////////////////////////

// tensors for internal 6j symbols in 15j
TENSOR_INIT(6j_3, 3) // 3 indices
TENSOR_INIT(6j_4, 4) // 4 indices
TENSOR_INIT(6j_5, 5) // 5 indices
TENSOR_INIT(6j_6, 6) // 6 indices

// macro for assembling internal 6j tensor
// NB: the 6j symmetries help a lot here (we can exchange
//     any columns and flip two columns for free)
#define FILL_6J(two_ka, two_ja, two_kb, two_jb, two_jc, tens) \
for (dspin two_kb = two_kb##_absmin; two_kb <= two_kb##_absmax; two_kb += 2) { \
for (dspin two_ka = two_ka##_absmin; two_ka <= two_ka##_absmax; two_ka += 2) { \
for (dspin two_x = two_x_absmin; two_x <= two_x_absmax; two_x += 2) { \
    double w6j = fw6jja(two_ka, two_ja, two_x, two_kb, two_jb, two_jc); \
    TENSOR_SET(w6j, tens, 3, DIV2(two_x-two_x_absmin), DIV2(two_ka-two_ka##_absmin), DIV2(two_kb-two_kb##_absmin)); \
} \
} \
} \

// Computes the five 6j tensors in the indices ls and ks and x. Look at the 15j function
// below for the order and the spins of the symbols.
// TODO: with MPI some time can be saved by just computing needed symbols...
//       but the time spent here is non-significant compared to vertex assembly
//       so it is probably not worth the effort
static inline void sl2cfoam_6j_tensors_vertex(dspin two_js[10], dspin two_i1_min, dspin two_i1_max, dspin two_Dl,
                                              tensor_ptr(6j_4)* sjA, tensor_ptr(6j_5)* sjB, tensor_ptr(6j_6)* sjC, 
                                              tensor_ptr(6j_5)* sjD, tensor_ptr(6j_4)* sjE, dspin bounds[10]) {

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    // compute ks range
    size_t i1dim, k2dim, k3dim, k4dim, k5dim;

    dspin two_k2_absmin, two_k2_absmax, two_k3_absmin, two_k3_absmax,
          two_k4_absmin, two_k4_absmax, two_k5_absmin, two_k5_absmax;

    dspin two_i1_absmin = two_i1_min;
    dspin two_i1_absmax = two_i1_max;
    i1dim = DIV2(two_i1_absmax - two_i1_absmin) + 1;

    find_k_absolute_bounds(&k2dim, &two_k2_absmin, &two_k2_absmax, two_j23, two_j24, two_j25, two_j12, two_Dl, 4);
    find_k_absolute_bounds(&k3dim, &two_k3_absmin, &two_k3_absmax, two_j34, two_j35, two_j13, two_j23, two_Dl, 3);
    find_k_absolute_bounds(&k4dim, &two_k4_absmin, &two_k4_absmax, two_j45, two_j14, two_j24, two_j34, two_Dl, 2);
    find_k_absolute_bounds(&k5dim, &two_k5_absmin, &two_k5_absmax, two_j15, two_j25, two_j35, two_j45, two_Dl, 1);

    // compute x range

    // x bounds good for all ks
    dspin two_x_absmin, two_x_absmax;

    dspin two_zero = is_integer(two_i1_absmin + two_j25) ? 0 : 1;

    two_x_absmin = max3(two_zero,     two_i1_absmin-two_j25-two_Dl, two_j25-two_i1_absmax);
    two_x_absmin = max3(two_x_absmin, two_k5_absmin-two_j14, two_j14-two_k5_absmax);
    two_x_absmin = max3(two_x_absmin, two_k4_absmin-two_j35-two_Dl, two_j35-two_k4_absmax);
    two_x_absmin = max3(two_x_absmin, two_k3_absmin-two_j24-two_Dl, two_j24-two_k3_absmax);
    two_x_absmin = max3(two_x_absmin, two_k2_absmin-two_j13, two_j13-two_k2_absmax);

    two_x_absmax = min(two_i1_absmax + two_j25 + two_Dl, two_j14 + two_k5_absmax);
    two_x_absmax = min(two_x_absmax, two_k4_absmax + two_j35 + two_Dl);
    two_x_absmax = min(two_x_absmax, two_k3_absmax + two_j24 + two_Dl);
    two_x_absmax = min(two_x_absmax, two_k2_absmax + two_j13);

    // set bounds to return
    bounds[0] = two_k2_absmin;
    bounds[1] = two_k2_absmax;
    bounds[2] = two_k3_absmin;
    bounds[3] = two_k3_absmax;
    bounds[4] = two_k4_absmin;
    bounds[5] = two_k4_absmax;
    bounds[6] = two_k5_absmin;
    bounds[7] = two_k5_absmax;
    bounds[8] = two_x_absmin;
    bounds[9] = two_x_absmax;

    // sometimes triangular inequalities leave nothing to compute
    // (e.g. in batching)
    if (two_x_absmax < two_x_absmin) {
        *sjA = NULL;
        *sjB = NULL;
        *sjC = NULL;
        *sjD = NULL;
        *sjE = NULL;
        return;
    }

    size_t xdim = DIV2(two_x_absmax - two_x_absmin) + 1;
    size_t ldim = DIV2(two_Dl) + 1;

    // tensors to return
    tensor_ptr(6j_4) tA;
    tensor_ptr(6j_5) tB;
    tensor_ptr(6j_6) tC;
    tensor_ptr(6j_5) tD;
    tensor_ptr(6j_4) tE;

    TENSOR_CREATE(6j_4, tA, 4, xdim, i1dim, k5dim, ldim);
    TENSOR_CREATE(6j_5, tB, 5, xdim, k5dim, k4dim, ldim, ldim);
    TENSOR_CREATE(6j_6, tC, 6, xdim, k4dim, k3dim, ldim, ldim, ldim);
    TENSOR_CREATE(6j_5, tD, 5, xdim, k3dim, k2dim, ldim, ldim);
    TENSOR_CREATE(6j_4, tE, 4, xdim, i1dim, k2dim, ldim);

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE)
    {
    #endif

    wig_thread_temp_init(CONFIG.max_two_spin);

    // tensor for 3 inner indices
    tensor_ptr(6j_3) tkx;

    ////////////////////////////////////////////////////////////////////
    // first tensor (x i1 k5 l25)
    ////////////////////////////////////////////////////////////////////

    TENSOR_CREATE(6j_3, tkx, 3, xdim, i1dim, k5dim);

    #ifdef USE_OMP
    #pragma omp for schedule(dynamic, 1)
    #endif
    for (dspin two_l25 = two_j25; two_l25 <= two_j25 + two_Dl; two_l25 += 2) {

        TENSOR_ZERO(tkx);
        FILL_6J(two_i1, two_l25, two_k5, two_j14, two_j15, tkx);

        double* dst = tA->d + TENSOR_INDEX(tA, 4, 0, 0, 0, DIV2(two_l25-two_j25));
        memcpy(dst, tkx->d, tkx->dim * sizeof(double));

    }

    TENSOR_FREE(tkx);

    ////////////////////////////////////////////////////////////////////
    // second tensor (x k5 k4 l35 l45)
    ////////////////////////////////////////////////////////////////////

    TENSOR_CREATE(6j_3, tkx, 3, xdim, k5dim, k4dim);

    #ifdef USE_OMP
    #pragma omp for collapse(2) schedule(dynamic, 1)
    #endif
    for (dspin two_l45 = two_j45; two_l45 <= two_j45 + two_Dl; two_l45 += 2) {
    for (dspin two_l35 = two_j35; two_l35 <= two_j35 + two_Dl; two_l35 += 2) {

        TENSOR_ZERO(tkx);
        FILL_6J(two_k5, two_j14, two_k4, two_l35, two_l45, tkx);

        double* dst = tB->d + TENSOR_INDEX(tB, 5, 0, 0, 0, DIV2(two_l35-two_j35), DIV2(two_l45-two_j45));
        memcpy(dst, tkx->d, tkx->dim * sizeof(double));

    }
    }

    TENSOR_FREE(tkx);

    ////////////////////////////////////////////////////////////////////
    // third tensor (x  k4 k3 l24 l34 l35)
    ////////////////////////////////////////////////////////////////////

    TENSOR_CREATE(6j_3, tkx, 3, xdim, k4dim, k3dim);

    #ifdef USE_OMP
    #pragma omp for  collapse(3) schedule(dynamic, 1)
    #endif
    for (dspin two_l35 = two_j35; two_l35 <= two_j35 + two_Dl; two_l35 += 2) {
    for (dspin two_l34 = two_j34; two_l34 <= two_j34 + two_Dl; two_l34 += 2) {
    for (dspin two_l24 = two_j24; two_l24 <= two_j24 + two_Dl; two_l24 += 2) {

        TENSOR_ZERO(tkx);
        FILL_6J(two_k4, two_l35, two_k3, two_l24, two_l34, tkx);

        double* dst = tC->d + TENSOR_INDEX(tC, 6, 0, 0, 0, DIV2(two_l24-two_j24), DIV2(two_l34-two_j34), DIV2(two_l35-two_j35));
        memcpy(dst, tkx->d, tkx->dim * sizeof(double));

    }
    }
    }

    TENSOR_FREE(tkx);

    ////////////////////////////////////////////////////////////////////
    // fourth tensor (x k3 k2 l23 l24)
    ////////////////////////////////////////////////////////////////////

    TENSOR_CREATE(6j_3, tkx, 3, xdim, k3dim, k2dim);

    #ifdef USE_OMP
    #pragma omp for collapse(2) schedule(dynamic, 1)
    #endif
    for (dspin two_l24 = two_j24; two_l24 <= two_j24 + two_Dl; two_l24 += 2) {
    for (dspin two_l23 = two_j23; two_l23 <= two_j23 + two_Dl; two_l23 += 2) {

        TENSOR_ZERO(tkx);
        FILL_6J(two_k3, two_l24, two_k2, two_j13, two_l23, tkx);

        double* dst = tD->d + TENSOR_INDEX(tD, 5, 0, 0, 0, DIV2(two_l23-two_j23), DIV2(two_l24-two_j24));
        memcpy(dst, tkx->d, tkx->dim * sizeof(double));

    }
    }

    TENSOR_FREE(tkx);

    ////////////////////////////////////////////////////////////////////
    // fifth tensor (x i1 k2 l25)
    ////////////////////////////////////////////////////////////////////

    TENSOR_CREATE(6j_3, tkx, 3, xdim, i1dim, k2dim);

    #ifdef USE_OMP
    #pragma omp for schedule(dynamic, 1)
    #endif
    for (dspin two_l25 = two_j25; two_l25 <= two_j25 + two_Dl; two_l25 += 2) {

        TENSOR_ZERO(tkx);
        FILL_6J(two_i1, two_l25, two_k2, two_j13, two_j12, tkx);

        double* dst = tE->d + TENSOR_INDEX(tE, 4, 0, 0, 0, DIV2(two_l25-two_j25));
        memcpy(dst, tkx->d, tkx->dim * sizeof(double));

    }

    TENSOR_FREE(tkx);

    wig_temp_free();

    #ifdef USE_OMP
    } // omp parallel
    #endif

    *sjA = tA;
    *sjB = tB;
    *sjC = tC;
    *sjD = tD;
    *sjE = tE;

}

// Computes an irreducbile 15j of first type (see Yutsis for convention).
static inline double sl2cfoam_w15j(dspin two_j12, dspin two_j13, dspin two_j14, dspin two_j15, dspin two_l23,
                                   dspin two_l24, dspin two_l25, dspin two_l34, dspin two_l35, dspin two_l45, 
                                   dspin two_i1, dspin two_k2, dspin two_k3, dspin two_k4, dspin two_k5) {

    dspin two_x, two_x_min, two_x_max;

    two_x_min = max(abs(two_i1-two_l25), abs(two_j14-two_k5));
    two_x_min = max(two_x_min, abs(two_k4-two_l35));
    two_x_min = max(two_x_min, abs(two_l24-two_k3));
    two_x_min = max(two_x_min, abs(two_k2-two_j13));

    two_x_max = min(two_i1+two_l25, two_j14+two_k5);
    two_x_max = min(two_x_max, two_k4+two_l35);
    two_x_max = min(two_x_max, two_l24+two_k3);
    two_x_max = min(two_x_max, two_k2+two_j13);

    if (two_x_max < two_x_min) {
        return 0;
    }

    long double res = 0;

    // TODO: continue if 0, maybe a cutoff is better?
    double w6j_1, w6j_2, w6j_3, w6j_4, w6j_5;

    for (two_x = two_x_min; two_x <= two_x_max; two_x += 2) {

        w6j_1 = fw6jja(two_i1, two_l25, two_x,
                       two_k5, two_j14, two_j15);
        if (w6j_1 == 0) continue;

        w6j_2 = fw6jja(two_j14, two_k5, two_x,
                       two_l35, two_k4, two_l45);
        if (w6j_2 == 0) continue;

        w6j_3 = fw6jja(two_k4, two_l35, two_x,
                       two_k3, two_l24, two_l34);
        if (w6j_3 == 0) continue;

        w6j_4 = fw6jja(two_l24, two_k3, two_x,
                       two_j13, two_k2, two_l23);
        if (w6j_4 == 0) continue;

        w6j_5 = fw6jja(two_k2, two_j13, two_x,
                       two_i1, two_l25, two_j12);
        if (w6j_5 == 0) continue;

        res += (long double)(DIM(two_x) * w6j_1 * w6j_2 * w6j_3 * w6j_4 * w6j_5);

    }

    int sign = real_negpow(two_j12 + two_j13 + two_j15 + two_j14 + two_l23 + 
                           two_l24 + two_l25 + two_l34 + two_l35 + two_l45 + 
                           two_i1 + two_k2 + two_k3 + two_k4 + two_k5);

    return (double)(sign * res);

}

// Computes a wigner 4j symbol (on-the-fly).
static inline double sl2cfoam_w4jm(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
			                       dspin two_m1, dspin two_m2, dspin two_m3, dspin two_m4,
				                   dspin two_i) {

	double w1, w2;
    w1 = fw3jja6(two_j1, two_j2, two_i, two_m1, two_m2, -two_m1-two_m2);
    w2 = fw3jja6(two_i, two_j3, two_j4, two_m1+two_m2, two_m3, two_m4);

    return real_negpow(two_i + two_m1 + two_m2) * w1 * w2;
    
}

#ifdef SL2CFOAM_SHORT_NAMES
#define init_wigxjpf_thread(...) sl2cfoam_init_wigxjpf_thread(__VA_ARGS__)
#define clear_wigxjpf_thread(...) sl2cfoam_clear_wigxjpf_thread(__VA_ARGS__)
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_JSYMBOLS_H__*/
