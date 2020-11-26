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

#include <complex.h>
#include <quadmath.h>
#include <omp.h>

#define MPFR_WANT_FLOAT128
#include <mpfr.h>
#include <mpc.h>

#include "integration_qagp.h"
#include "common.h"
#include "utils.h"
#include "error.h"
#include "verb.h"
#include "cgamma.h"
#include "sl2cfoam_tensors.h"
#include "jsymbols.h"
#include "dsmall.h"
#include "blas_wrapper.h"

TENSOR_INIT(dsmall_integral, 4);
TENSOR_INIT(boost_wig, 6);

typedef struct __dsmall_params {
    int prec;
    double rho;
    dspin two_j;
    dspin two_l;
    dspin two_p;
    mpc_ptr* Ym;
    mpc_ptr* Yn;
    __float128 dx_min;
} __dsmall_params;

static void dsmall_prod_qagp(__complex128 ds[], __float128 xs[], size_t N, void* params) {

    __dsmall_params* p1 = ((__dsmall_params**)params)[0];
    __dsmall_params* p2 = ((__dsmall_params**)params)[1];
    __dsmall_params* p3 = ((__dsmall_params**)params)[2];
    __dsmall_params* p4 = ((__dsmall_params**)params)[3];

    __complex128 d1[N];
    __complex128 d2[N];
    __complex128 d3[N];
    __complex128 d4[N];

    sl2cfoam_dsmall(d1, xs, N, p1->prec, p1->Ym, p1->Yn, NULL,
                    p1->rho, p1->two_j, p1->two_j, p1->two_l, p1->two_p);

    sl2cfoam_dsmall(d2, xs, N, p2->prec, p2->Ym, p2->Yn, NULL,
                    p2->rho, p2->two_j, p2->two_j, p2->two_l, p2->two_p);

    sl2cfoam_dsmall(d3, xs, N, p3->prec, p3->Ym, p3->Yn, NULL,
                    p3->rho, p3->two_j, p3->two_j, p3->two_l, p3->two_p);

    sl2cfoam_dsmall(d4, xs, N, p4->prec, p4->Ym, p4->Yn, NULL,
                    p4->rho, p4->two_j, p4->two_j, p4->two_l, p4->two_p);

    for (int i = 0; i < N; i++) {
        ds[i] = d1[i] * d2[i] * d3[i] * d4[i];
    }

}

// TODO: consider parallelize this loop (but it's probably faster than parallel latency)
static void dsmall_measure_qagp(__complex128 ds[], __float128 xs[], size_t N, void* params) {

    const __float128 k =  M_1_PIq / 16.0Q;

    __float128 x, invx, mx;
    for (int i = 0; i < N; i++) {

        x = xs[i];
        invx = 1.0Q / x;

        mx = invx;
        mx *= SQ(invx - 1.0Q);
        mx *= SQ(x + 1.0Q);
        mx *= k;

        ds[i] = mx;

    }

}

// computes all the parameterd and the Y coefficients for a dsmall call
static inline void dsmall_fill_params(__dsmall_params* p, __float128 dx_min,
                                      double rho, dspin two_j, dspin two_l, dspin two_p) {
    
    // compute precision
    int prec_base, prec_add;

    // there should be enough significant digits to cancel the prefactor 
    // (1-(1-dx)^2)^-(j1+j2+1) ~ (2*dx)^-(j1+j2+1) ...
    prec_add = - (int)((DIV2(two_j+two_l)+1) * log2q(2 * dx_min));

    // and some more of course
    // base precision is found by "trial and error"
    if (two_l <= 50) {
        prec_base = 64;
    } else if (two_l <= 100) {
        prec_base = 128;
    } else {
        prec_base = 256;
    }

    // total precision;
    int precision = prec_base + prec_add;

    // compute the Y coefficients
    int Jmp = DIV2(abs(two_j - two_p)); // |J-p|
    int Jpp = DIV2(abs(two_j + two_p)); // |J+p|
    int m_max = DIV2(two_j + two_l) - Jmp;
    int n_max = DIV2(two_j + two_l) - Jpp;

    mpc_ptr* Ym;
    mpc_ptr* Yn;

    Ym = (mpc_ptr*) malloc((m_max+1) * sizeof(mpc_ptr));
    for (int m = 0; m <= m_max; m++) {
        Ym[m] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(Ym[m], precision);
    }

    Yn = (mpc_ptr*) malloc((n_max+1) * sizeof(mpc_ptr));
    for (int n = 0; n <= n_max; n++) {
        Yn[n] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(Yn[n], precision);
    }

    // NOTICE the order of the arguments
    sl2cfoam_dsmall_Yc(Ym, precision,  rho, two_j, two_l, two_j,  two_p);
    sl2cfoam_dsmall_Yc(Yn, precision, -rho, two_j, two_j, two_l, -two_p);

    p->prec = precision;
    p->rho = rho;
    p->two_j = two_j;
    p->two_l = two_l;
    p->two_p = two_p;
    p->Ym = Ym;
    p->Yn = Yn;
    p->dx_min = dx_min;

}

static inline void dsmall_free_params(__dsmall_params* p) {

    dspin two_j, two_l, two_p;
    two_j = p->two_j;
    two_l = p->two_l;
    two_p = p->two_p;

    int Jmp = DIV2(abs(two_j - two_p)); // |J-p|
    int Jpp = DIV2(abs(two_j + two_p)); // |J+p|
    int m_max = DIV2(two_j + two_l) - Jmp;
    int n_max = DIV2(two_j + two_l) - Jpp;

    // clear the Y coefficients
    for (int m = 0; m <= m_max; m++) {
        mpc_clear(p->Ym[m]);
        free(p->Ym[m]);
    }
    free(p->Ym);

    for (int n = 0; n <= n_max; n++) {
        mpc_clear(p->Yn[n]);
        free(p->Yn[n]);
    }
    free(p->Yn);

}

static inline __complex128 dsmall_phase(dspin two_j, dspin two_l, cgamma_lanczos* lanczos) {

    spin j = SPIN(two_j);
    spin l = SPIN(two_l);
    double rho = RHO(j);

    // compute Speziale's phase
    __complex128 sph;

    mpc_t cg, z1, z2;
    mpfr_t gabs;
    mpc_init2(cg, MPBITS);
    mpc_init2(z1, MPBITS);
    mpc_init2(z2, MPBITS);
    mpfr_init2(gabs, MPBITS);

    mpc_set_d_d(z1, j + 1.0, rho, MPC_RNDNN);
    sl2cfoam_complex_gamma(z2, z1, lanczos);

    mpc_abs(gabs, z2, MPFR_RNDN);
    mpc_div_fr(z2, z2, gabs, MPC_RNDNN);

    mpc_set_d_d(z1, l + 1.0, -rho, MPC_RNDNN);
    sl2cfoam_complex_gamma(cg, z1, lanczos);

    mpc_set_d_d(z1, l + 1.0, +rho, MPC_RNDNN);
    sl2cfoam_complex_gamma(z1, z1, lanczos);

    mpc_abs(gabs, z1, MPFR_RNDN);
    mpc_div_fr(cg, cg, gabs, MPC_RNDNN);

    mpc_mul(z2, z2, cg, MPC_RNDNN);

    sph = mpfr_get_float128(z2->re, MPFR_RNDN) + I * mpfr_get_float128(z2->im, MPFR_RNDN);
    sph *= complex_negpow(DIV2(two_l-two_j));

    mpc_clear(cg);
    mpc_clear(z1);
    mpc_clear(z2);
    mpfr_clear(gabs);

    return sph;
    
}

sl2cfoam_dmatrix sl2cfoam_b4_accurate(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
                                      dspin two_l1, dspin two_l2, dspin two_l3, dspin two_l4,
                                      dspin two_i_min, dspin two_i_max,
                                      dspin two_k_min, dspin two_k_max) {

    dspin two_jis[4] = { two_j1, two_j2, two_j3, two_j4 };
    dspin two_lis[4] = { two_l1, two_l2, two_l3, two_l4 };

    size_t dimp1 = DIM(two_j1);
    size_t dimp2 = DIM(two_j2);
    size_t dimp3 = DIM(two_j3);
    size_t dimp4 = DIM(two_j4);

    // parameters for integration
    int qagp_max_nevals;
    double qagp_tolerance;

    switch (ACCURACY)
    {
    case SL2CFOAM_ACCURACY_NORMAL:
        qagp_max_nevals = 1000;
        qagp_tolerance = 1e-4;
        break;

    case SL2CFOAM_ACCURACY_HIGH:
        qagp_max_nevals = 4000;
        qagp_tolerance = 1e-6;
        break;

    case SL2CFOAM_ACCURACY_VERYHIGH:
        qagp_max_nevals = 10000;
        qagp_tolerance = 1e-9;
        break;
    
    default:
        error("wrong accuracy value");
    }

    // tensor for dsmall integrals
    tensor_ptr(dsmall_integral) dtens;
    TENSOR_CREATE(dsmall_integral, dtens, 4, dimp1, dimp2, dimp3, dimp4);

    // compute al necessary coefficients and phases
    // for each of the 4 dsmall loop over all possible ps
    __complex128 phases[4];
    __dsmall_params dp1[dimp1];
    __dsmall_params dp2[dimp2];
    __dsmall_params dp3[dimp3];
    __dsmall_params dp4[dimp4];

    __dsmall_params* dp_all[4] = { dp1, dp2, dp3, dp4 };

    cgamma_lanczos lanczos;
    sl2cfoam_cgamma_lanczos_fill(&lanczos);

    dspin two_ji, two_li;
    spin ji;
    double rho;
    for (int dsmall_index = 0; dsmall_index < 4; dsmall_index++) {

        two_ji = two_jis[dsmall_index];
        two_li = two_lis[dsmall_index];

        ji = SPIN(two_ji);
        rho = RHO(ji);

        phases[dsmall_index] = dsmall_phase(two_ji, two_li, &lanczos);

        #ifdef USE_OMP
        #pragma omp parallel for if(OMP_PARALLELIZE)
        #endif
        for (dspin two_p = -two_ji; two_p <= two_ji; two_p +=2) {

            // compute dx_min for this p
            __float128 dx_min;
            if (abs(two_p) <= two_ji/3) {
                dx_min = 1e-6;
            } else if (abs(two_p) <= 2*two_ji/3) {
                dx_min = 1e-8;
            } else {
                dx_min = 1e-9;
            }

            dsmall_fill_params(&(dp_all[dsmall_index][DIV2(two_p+two_ji)]), dx_min, rho, two_ji, two_li, two_p);

        }

    }

    sl2cfoam_cgamma_lanczos_free(&lanczos);

    // global phase
    __complex128 phase_p = phases[0] * phases[1] * phases[2] * phases[3];

    // compute the integrals

    #ifdef USE_OMP
    #pragma omp parallel for collapse(3) schedule(dynamic, 2) if(OMP_PARALLELIZE)
    #endif
    for (dspin two_p4 = -two_j4; two_p4 <= two_j4; two_p4 += 2) {
    for (dspin two_p3 = -two_j3; two_p3 <= two_j3; two_p3 += 2) {
    for (dspin two_p2 = -two_j2; two_p2 <= two_j2; two_p2 += 2) {

        dspin two_p1 = - two_p4 - two_p3 - two_p2;
        if (two_p1 < -two_j1 || two_p1 > two_j1) {
            continue;
        }

        sl2cfoam_integration_function f_qagp;
        f_qagp.function = &dsmall_prod_qagp;
        f_qagp.measure = &dsmall_measure_qagp;

        __dsmall_params* dp1 = &(dp_all[0][DIV2(two_p1+two_j1)]);
        __dsmall_params* dp2 = &(dp_all[1][DIV2(two_p2+two_j2)]);
        __dsmall_params* dp3 = &(dp_all[2][DIV2(two_p3+two_j3)]);
        __dsmall_params* dp4 = &(dp_all[3][DIV2(two_p4+two_j4)]);

        __dsmall_params* dp[4] = { dp1, dp2, dp3, dp4 };
        f_qagp.params = dp;

        __float128 dx_min;
        dx_min = fminq(dp1->dx_min, dp2->dx_min);
        dx_min = fminq(dx_min, dp3->dx_min);
        dx_min = fminq(dx_min, dp4->dx_min);

        __complex128 res_qagp;
        __float128 abserr_qagp;
        size_t nevals_qagp;
        sl2cfoam_qagp_retcode rcode;
        rcode = sl2cfoam_qagp(&f_qagp, qagp_tolerance, &res_qagp, &abserr_qagp, &nevals_qagp, dx_min, qagp_max_nevals);

        if (rcode != QAGP_SUCCESS) {
            verb(SL2CFOAM_VERBOSE_LOW, "adaptive integral (p1, p2, p3, p4) = (%d, %d, %d, %d) above tolerance, code = %d\n", two_p1, two_p2, two_p3, two_p4, rcode);
        }

        if (abserr_qagp >= cabsq(res_qagp)) {
            verb(SL2CFOAM_VERBOSE_LOW, "adaptive integral (p1, p2, p3, p4) = (%d, %d, %d, %d) absolute error >= |result|, code = %d\n", two_p1, two_p2, two_p3, two_p4, rcode);
        }

        // multiply global phases
        __complex128 integ_p = res_qagp * phase_p;

        // now the integrals should be real, truncate to double
        double dip = (double)crealq(integ_p);

        // set the value in the integral tensor
        TENSOR_SET(dip, dtens, 4, DIV2(two_p1+two_j1), DIV2(two_p2+two_j2), DIV2(two_p3+two_j3), DIV2(two_p4+two_j4));

    } // p2
    } // p3
    } // p4

    // free the memory of the coefficients
    for (int dsmall_index = 0; dsmall_index < 4; dsmall_index++) {

        two_ji = two_jis[dsmall_index];
        for (dspin two_p = -two_ji; two_p <= two_ji; two_p +=2) {
            dsmall_free_params(&(dp_all[dsmall_index][DIV2(two_p+two_ji)]));
        }

    }

    // tensor for wigner symbols

    int dimi = DIV2(two_i_max-two_i_min) + 1;
    int dimk = DIV2(two_k_max-two_k_min) + 1;

    tensor_ptr(boost_wig) wtens;
    TENSOR_CREATE(boost_wig, wtens, 6, dimi, dimk, dimp1, dimp2, dimp3, dimp4);

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE)
    {
    #endif

    wig_thread_temp_init(CONFIG.max_two_spin);

    #ifdef USE_OMP
    #pragma omp for collapse(3)
    #endif
    for (dspin two_p4 = -two_j4; two_p4 <= two_j4; two_p4 += 2) {
    for (dspin two_p3 = -two_j3; two_p3 <= two_j3; two_p3 += 2) {
    for (dspin two_p2 = -two_j2; two_p2 <= two_j2; two_p2 += 2) {

        dspin two_p1 = - two_p4 - two_p3 - two_p2;
        if (two_p1 < -two_j1 || two_p1 > two_j1) {
            continue;
        }

        for (dspin two_k = two_k_min; two_k <= two_k_max; two_k += 2) {
        for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {
        
            int ii = DIV2(two_i-two_i_min);
            int ki = DIV2(two_k-two_k_min);

                double w4jj = sl2cfoam_w4jm(two_j1, two_j2, two_j3, two_j4,
                                            two_p1, two_p2, two_p3, two_p4, two_i);

                double w4jl = sl2cfoam_w4jm(two_l1, two_l2, two_l3, two_l4,
                                            two_p1, two_p2, two_p3, two_p4, two_k);

                TENSOR_SET(w4jj * w4jl, wtens, 6, ii, ki, DIV2(two_p1+two_j1), DIV2(two_p2+two_j2), DIV2(two_p3+two_j3), DIV2(two_p4+two_j4));

        } // i
        } // k

    } // p2
    } // p3
    } // p4

    wig_temp_free();

    #ifdef USE_OMP
    } // parallel
    #endif

    // matrix in indices (i, k)
    sl2cfoam_dmatrix b4 = dmatrix_alloc(dimi, dimk);

    // contract the dsmall ad 4js tensors along indices p1, p2, p3, p4
    long dgemv_m = dimi * dimk;
    long dgemv_n = dimp1 * dimp2 * dimp3 * dimp4;
    BLASW_DGEMV(dgemv_m, dgemv_n, wtens->d, dtens->d, b4);

    TENSOR_FREE(dtens);
    TENSOR_FREE(wtens);

    return b4;

}
