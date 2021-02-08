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

#include <math.h>
#include <complex.h>
#include <quadmath.h>
#include <omp.h>

#define MPFR_WANT_FLOAT128
#include <mpfr.h>
#include <mpc.h>

#include "common.h"
#include "utils.h"
#include "error.h"
#include "verb.h"
#include "cgamma.h"
#include "sl2cfoam_tensors.h"
#include "dsmall.h"
#include "integration_gk.h"
#include "blas_wrapper.h"
#include "wigxjpf.h"
#include "fastwigxj.h"

TENSOR_INIT(dsmall_integral, 4);
TENSOR_INIT(boost_wig, 3);

static void dsmall_measure(double ds[], double xs[], size_t N) {

    const double k =  M_1_PI / 16.0;

    double x, invx, mx;
    for (int i = 0; i < N; i++) {

        x = xs[i];
        invx = 1.0 / x;

        mx = invx;
        mx *= SQ(invx - 1.0);
        mx *= SQ(x + 1.0);
        mx *= k;

        ds[i] = mx;

    }

}

static inline __complex128 dsmall_phase(dspin two_j, dspin two_l, cgamma_lanczos* lanczos) {

    spin j = SPIN(two_j);
    spin l = SPIN(two_l);
    double rho = RHO(j);

    // compute Speziale's phase
    __complex128 sph;

    mpc_t z1, z2;
    mpfr_t gabs;
    mpc_init2(z1, MPBITS);
    mpc_init2(z2, MPBITS);
    mpfr_init2(gabs, MPBITS);

    mpc_set_d_d(z1, j + 1.0, rho, MPC_RNDNN);
    sl2cfoam_complex_gamma(z2, z1, lanczos);

    mpc_abs(gabs, z2, MPFR_RNDN);
    mpc_div_fr(z2, z2, gabs, MPC_RNDNN);

    mpc_set_d_d(z1, l + 1.0, -rho, MPC_RNDNN);
    sl2cfoam_complex_gamma(z1, z1, lanczos);

    mpc_abs(gabs, z1, MPFR_RNDN);
    mpc_div_fr(z1, z1, gabs, MPC_RNDNN);

    mpc_mul(z2, z2, z1, MPC_RNDNN);

    sph = mpfr_get_float128(z2->re, MPFR_RNDN) + I * mpfr_get_float128(z2->im, MPFR_RNDN);
    sph *= complex_negpow(DIV2(two_l-two_j));

    mpc_clear(z1);
    mpc_clear(z2);
    mpfr_clear(gabs);

    return sph;
        
}

static inline void dsmall_prefactors(mpc_ptr* rop, __float128* xs, int nxs, int precision,
                                     dspin two_j, dspin two_l, double rho) {

    mpc_t emi_r_rho_ab;
    mpc_t pref;

    mpfr_t emr_ab;
    mpfr_t mr_rho_ab;

    mpc_init2(emi_r_rho_ab, precision);
    mpc_init2(pref, precision);

    mpfr_init2(emr_ab, precision);
    mpfr_init2(mr_rho_ab, precision);

    for (int i = 0; i < nxs; i++) {

        // exp(-r)
        mpfr_set_float128(emr_ab, xs[i], MPFR_RNDN);

        // get exp(-i*r*rho)
        mpfr_log(mr_rho_ab, emr_ab, MPFR_RNDN);
        mpfr_mul_d(mr_rho_ab, mr_rho_ab, rho, MPFR_RNDN);
        mpfr_sin_cos(emi_r_rho_ab->im, emi_r_rho_ab->re, mr_rho_ab, MPFR_RNDN);

        // prefactor
        mpfr_sqr(pref->re, emr_ab, MPFR_RNDN);
        mpfr_set_ui(pref->im, 0, MPFR_RNDN);
        mpfr_ui_sub(pref->re, 1, pref->re, MPFR_RNDN);
        mpfr_pow_si(pref->re, pref->re, -DIV2(two_j+two_l) - 1, MPFR_RNDN);

        mpc_mul(rop[i], pref, emi_r_rho_ab, MPC_RNDNN);

    }

    mpc_clear(emi_r_rho_ab);
    mpc_clear(pref);

    mpfr_clear(mr_rho_ab);
    mpfr_clear(emr_ab);

}

// IF there is no parallelization above (and nested parallelization is disabled)
// then the parallelization strategy is: the computation of dsmall is parallelized
// inside the dsmall function. Then the rest is parallelized over the p indices.
sl2cfoam_dmatrix sl2cfoam_b4(dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
                             dspin two_l1, dspin two_l2, dspin two_l3, dspin two_l4) {

    dspin two_jis[4] = { two_j1, two_j2, two_j3, two_j4 };
    dspin two_lis[4] = { two_l1, two_l2, two_l3, two_l4 };

    size_t dimp1 = DIM(two_j1);
    size_t dimp2 = DIM(two_j2);
    size_t dimp3 = DIM(two_j3);
    size_t dimp4 = DIM(two_j4);

    // fix global grid

    dspin two_l_max = max4(two_l1, two_l2, two_l3, two_l4);
    
    // set the number of intervals for integration
    // (next variable looks ok even as low as 1 for most cases...)

    double interval_mult;

    switch (ACCURACY)
    {
    case SL2CFOAM_ACCURACY_NORMAL:
        interval_mult = 1.5;
        break;

    case SL2CFOAM_ACCURACY_HIGH:
        interval_mult = 2.5;
        break;

    case SL2CFOAM_ACCURACY_VERYHIGH:
        interval_mult = 4.0;
        break;
    
    default:
        error("wrong accuracy value");
    }

    // TODO: better study of this criterion
    //       maybe no less than 2 intervals (~120 points) to begin with?
    double glmax = SPIN(two_l_max) * fmax(1.0, sqrt(IMMIRZI));
    int intervals = max(1, floor(interval_mult * sqrt(glmax)));
    int nxs = GK_POINTS * intervals;

    __float128* grid = sl2cfoam_grid_harmonic(intervals);
    __float128* qxs = sl2cfoam_gk_grid_abscissae(intervals, grid);
    __float128 dx_min = qxs[1];

    double xs[nxs];
    for (int i = 0; i < nxs; i++) {
        xs[i] = (double)qxs[i];
    }

    // relative error threshold in integration for prompting a warning
    const double gk_tol = 0.01;

    // the dsmall integral can be exactly zero by certain symmetries
    // so if we find that the error is large and the value is very small
    // then it is probably just numerical noise
    const double integral_zero = 1e-16;

    // compute measure
    double measure[nxs];
    dsmall_measure(measure, xs, nxs);

    // compute weights
    double* wgks = sl2cfoam_gk_grid_weights(intervals, grid);
    double* wgs = sl2cfoam_gk_grid_weights_gauss(intervals, grid);

    // arrays with dsmall values for each p and r
    sl2cfoam_cmatrix dp1 = cmatrix_alloc(nxs, dimp1);
    sl2cfoam_cmatrix dp2 = cmatrix_alloc(nxs, dimp2);
    sl2cfoam_cmatrix dp3 = cmatrix_alloc(nxs, dimp3);
    sl2cfoam_cmatrix dp4 = cmatrix_alloc(nxs, dimp4);

    sl2cfoam_cmatrix dpall[4] = { dp1, dp2, dp3, dp4 };

    cgamma_lanczos lanczos;
    sl2cfoam_cgamma_lanczos_fill(&lanczos);

    dspin two_ji, two_li;
    for (int dsmall_index = 0; dsmall_index < 4; dsmall_index++) {

        two_ji = two_jis[dsmall_index];
        two_li = two_lis[dsmall_index];

        spin ji = SPIN(two_ji);
        double rho = RHO(ji);

        // fix precision for this dsmall
        int prec_base, prec_add;

        // there should be enough significant digits to cancel the prefactor 
        // (1-(1-dx)^2)^-(j1+j2+1) ~ (2*dx)^-(j1+j2+1) ...
        prec_add = - (int)((DIV2(two_ji+two_li)+1) * log2q(2 * dx_min));

        // and some more of course
        // base precision is found by "trial and error"
        if (two_li <= 50) {
            prec_base = 64;
        } else if (two_li <= 100) {
            prec_base = 128;
        } else {
            prec_base = 256;
        }

        // total precision;
        int precision = prec_base + prec_add;

        // precompute expensive prefactors
        mpc_ptr prefactors[nxs];
        for (int i = 0; i < nxs; i++) {
            prefactors[i] = (mpc_ptr)malloc(sizeof(mpc_t));
            mpc_init2(prefactors[i], precision);
        }

        dsmall_prefactors(prefactors, qxs, nxs, precision, two_ji, two_li, rho);
        
        // compute Speziale's phase
        __complex128 sph = dsmall_phase(two_ji, two_li, &lanczos);
        int mph = real_negpow(two_ji-two_li);

        sl2cfoam_cmatrix dpi = dpall[dsmall_index];
        sl2cfoam_cvector dpip;
        sl2cfoam_cvector dpipm;

        __complex128 ds[nxs];

        for (dspin two_p = is_integer(two_ji) ? 0 : 1; two_p <= two_ji; two_p += 2) {

            // compute the Y coefficients
            int Jmp = DIV2(abs(two_ji - two_p)); // |J-p|
            int Jpp = DIV2(abs(two_ji + two_p)); // |J+p|
            int m_max = DIV2(two_ji + two_li) - Jmp;
            int n_max = DIV2(two_ji + two_li) - Jpp;

            mpc_ptr Ym[m_max+1];
            mpc_ptr Yn[n_max+1];

            for (int m = 0; m <= m_max; m++) {
                Ym[m] = (mpc_ptr)malloc(sizeof(mpc_t));
                mpc_init2(Ym[m], precision);
            }

            for (int n = 0; n <= n_max; n++) {
                Yn[n] = (mpc_ptr)malloc(sizeof(mpc_t));
                mpc_init2(Yn[n], precision);
            }

            sl2cfoam_dsmall_Yc(Ym, precision,  rho, two_ji, two_li, two_ji,  two_p);
            sl2cfoam_dsmall_Yc(Yn, precision, -rho, two_ji, two_ji, two_li, -two_p);

            // compute dsmall at all points
            sl2cfoam_dsmall(ds, qxs, nxs, precision, Ym, Yn, prefactors, 
                            rho, two_ji, two_ji, two_li, two_p);

            // multiply phase for p and set
            dpip = matrix_column(dpi, nxs, DIV2(two_p+two_ji));
            for (int i = 0; i < nxs; i++) {
                dpip[i] = (double complex)(ds[i] * sph);
            }

            // set the values for -p
            dpipm = matrix_column(dpi, nxs, DIV2(-two_p+two_ji));
            for (int i = 0; i < nxs; i++) {
                dpipm[i] = (double complex)(mph * conjq(dpip[i]));
            }

            // clear the Y coefficients
            for (int m = 0; m <= m_max; m++) {
                mpc_clear(Ym[m]);
                free(Ym[m]);
            }

            for (int n = 0; n <= n_max; n++) {
                mpc_clear(Yn[n]);
                free(Yn[n]);
            }

        } // p

        // clear prefactors
        for (int i = 0; i < nxs; i++) {
            mpc_clear(prefactors[i]);
            free(prefactors[i]);
        }

    } // loop over 4 dsmalls

    sl2cfoam_cgamma_lanczos_free(&lanczos);

    // tensor for dsmall integrals
    tensor_ptr(dsmall_integral) dtens;
    TENSOR_CREATE(dsmall_integral, dtens, 4, dimp1, dimp2, dimp3, dimp4);

    int wcount = 0;
    const int warn_max_per_thread = 10;

    #ifdef USE_OMP
    #pragma omp parallel for collapse(3) private(wcount) if(OMP_PARALLELIZE)
    #endif
    for (dspin two_p4 = -two_j4; two_p4 <= two_j4; two_p4 += 2) {
    for (dspin two_p3 = -two_j3; two_p3 <= two_j3; two_p3 += 2) {
    for (dspin two_p2 = -two_j2; two_p2 <= two_j2; two_p2 += 2) {

        dspin two_p1 = - two_p4 - two_p3 - two_p2;
        if (two_p1 < -two_j1 || two_p1 > two_j1) {
            continue;
        }

        int p1i, p2i, p3i, p4i;
        p1i = DIV2(two_p1+two_j1);
        p2i = DIV2(two_p2+two_j2);
        p3i = DIV2(two_p3+two_j3);
        p4i = DIV2(two_p4+two_j4);

        // multiply and then take real part
        double prod[nxs];
        double complex cprod;
        for (int i = 0; i < nxs; i++) {

            cprod = (matrix_column(dp1, nxs, p1i))[i] * 
                    (matrix_column(dp2, nxs, p2i))[i] *
                    (matrix_column(dp3, nxs, p3i))[i] *
                    (matrix_column(dp4, nxs, p4i))[i];
            prod[i] = creal(cprod);

        }

        double res, abserr, absres;

        res = sl2cfoam_gk_grid(intervals, prod, measure, wgks, wgs, &abserr);

        absres = fabs(res);
        if (abserr >= gk_tol * absres) {

            if (absres < integral_zero) {

                // this is numerical noise
                // do not bother and leave the integral to 0 exactly
                continue;

            }

            wcount++;
            if (wcount <= warn_max_per_thread) {

                // probably worth a warning
                if (wcount < warn_max_per_thread) {
                    warning("gk integral (pi) = (%d, %d, %d, %d) relative error = %.3g >= %.3g (integral = %.3g)",
                        two_p1, two_p2, two_p3, two_p4, abserr/absres, gk_tol, res);
                } else {
                    warning("gk integral (pi) = (%d, %d, %d, %d) relative error = %.3g >= %.3g (integral = %.3g)\n"
                            "(further warnings for integration on this thread not shown...)",
                            two_p1, two_p2, two_p3, two_p4, abserr/absres, gk_tol, res);
                }
                
            }
            
        }

        // set real part (should be) in dsmall tensor
        TENSOR_SET(res, dtens, 4, p1i, p2i, p3i, p4i);

    } // p2
    } // p3
    } // p4

    // compute 4jm tensors

    dspin two_i_min = max(abs(two_j1-two_j2), abs(two_j3-two_j4));
    dspin two_i_max = min(two_j1+two_j2, two_j3+two_j4);
    dspin two_k_min = max(abs(two_l1-two_l2), abs(two_l3-two_l4));
    dspin two_k_max = min(two_l1+two_l2, two_l3+two_l4);

    int dimi = DIV2(two_i_max-two_i_min) + 1;
    int dimk = DIV2(two_k_max-two_k_min) + 1;

    // tensors for wigner symbols
    tensor_ptr(boost_wig) wt_ip1p2;
    tensor_ptr(boost_wig) wt_ip3p4;
    tensor_ptr(boost_wig) wt_kp1p2;
    tensor_ptr(boost_wig) wt_kp3p4;
    TENSOR_CREATE(boost_wig, wt_ip1p2, 3, dimi, dimp1, dimp2);
    TENSOR_CREATE(boost_wig, wt_ip3p4, 3, dimi, dimp3, dimp4);
    TENSOR_CREATE(boost_wig, wt_kp1p2, 3, dimk, dimp1, dimp2);
    TENSOR_CREATE(boost_wig, wt_kp3p4, 3, dimk, dimp3, dimp4);

    // final b4 matrix in indices (i, k)
    sl2cfoam_dmatrix b4 = dmatrix_alloc(dimi, dimk);

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE)
    {
    #endif

    wig_thread_temp_init(CONFIG.max_two_spin);

    #ifdef USE_OMP
    #pragma omp for collapse(3)
    #endif
    for (dspin two_p2 = -two_j2; two_p2 <= two_j2; two_p2 += 2) {
    for (dspin two_p1 = -two_j1; two_p1 <= two_j1; two_p1 += 2) {
    for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {

        int ii = DIV2(two_i-two_i_min);

        double w3j = fw3jja6(two_j1, two_j2, two_i,
                            two_p1, two_p2, -two_p1-two_p2);
        TENSOR_SET(w3j, wt_ip1p2, 3, ii, DIV2(two_p1+two_j1), DIV2(two_p2+two_j2));

    } // i
    } // p1
    } // p2

    #ifdef USE_OMP
    #pragma omp for collapse(3)
    #endif
    for (dspin two_p4 = -two_j4; two_p4 <= two_j4; two_p4 += 2) {
    for (dspin two_p3 = -two_j3; two_p3 <= two_j3; two_p3 += 2) {
    for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {

        int ii = DIV2(two_i-two_i_min);

        double w3j = fw3jja6(two_i, two_j3, two_j4,
                                 -two_p3-two_p4, two_p3, two_p4);
        TENSOR_SET(w3j, wt_ip3p4, 3, ii, DIV2(two_p3+two_j3), DIV2(two_p4+two_j4));

    } // i
    } // p3
    } // p4

    #ifdef USE_OMP
    #pragma omp for collapse(3)
    #endif
    for (dspin two_p2 = -two_j2; two_p2 <= two_j2; two_p2 += 2) {   
    for (dspin two_p1 = -two_j1; two_p1 <= two_j1; two_p1 += 2) {
    for (dspin two_k = two_k_min; two_k <= two_k_max; two_k += 2) {

        int ki = DIV2(two_k-two_k_min);

        double w3j = fw3jja6(two_l1, two_l2, two_k,
                            two_p1, two_p2, -two_p1-two_p2);
        TENSOR_SET(w3j, wt_kp1p2, 3, ki, DIV2(two_p1+two_j1), DIV2(two_p2+two_j2));

    } // k
    } // p1
    } // p2

    #ifdef USE_OMP
    #pragma omp for collapse(3)
    #endif
    for (dspin two_p4 = -two_j4; two_p4 <= two_j4; two_p4 += 2) {
    for (dspin two_p3 = -two_j3; two_p3 <= two_j3; two_p3 += 2) {
    for (dspin two_k = two_k_min; two_k <= two_k_max; two_k += 2) {

        int ki = DIV2(two_k-two_k_min);

        double w3j = fw3jja6(two_k, two_l3, two_l4,
                            -two_p3-two_p4, two_p3, two_p4);
        TENSOR_SET(w3j, wt_kp3p4, 3, ki, DIV2(two_p3+two_j3), DIV2(two_p4+two_j4));

    } // k
    } // p3
    } // p4

    // now loop over all possible values
    // loop over ps and then over i,k
    // I need a temporary array for each thread for reduction
    // it's faster than reversing the loops
    // (probably due to NUMA access and first-touch policy with dtens tensor)

    sl2cfoam_dmatrix b4_thread = dmatrix_alloc(dimi, dimk);

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

        int p1i, p2i, p3i, p4i;
        p1i = DIV2(two_p1+two_j1);
        p2i = DIV2(two_p2+two_j2);
        p3i = DIV2(two_p3+two_j3);
        p4i = DIV2(two_p4+two_j4);
        
        double integral;

        integral = TENSOR_GET(dtens, 4, p1i, p2i, p3i, p4i);

        if (integral == 0.0) continue;

        for (dspin two_k = two_k_min; two_k <= two_k_max; two_k += 2) {
        for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {

            int ii = DIV2(two_i-two_i_min);
            int ki = DIV2(two_k-two_k_min);

            // check if boost is 0 by symmetries of the 3js
            if (((two_j1 == two_j2 && two_l2 == two_l1) || (two_j3 == two_j4  && two_l4 == two_l3))
                && (two_k - two_i) % 4 != 0) {

                    continue;

            }

            double w3j1, w3j2, w3j3, w3j4;
            w3j1 = TENSOR_GET(wt_ip1p2, 3, ii, p1i, p2i);
            w3j2 = TENSOR_GET(wt_ip3p4, 3, ii, p3i, p4i);
            w3j3 = TENSOR_GET(wt_kp1p2, 3, ki, p1i, p2i);
            w3j4 = TENSOR_GET(wt_kp3p4, 3, ki, p3i, p4i);

            double kisum;
            kisum = real_negpow(two_i + two_k + 2*(two_p1 + two_p2)) 
                    * w3j1 * w3j2 * w3j3 * w3j4 * integral;
            
            matrix_get(b4_thread, dimi, ii, ki) += kisum;

        } // i
        } // k

    } // p2
    } // p3
    } // p4

    // reduce over all threads
    #ifdef USE_OMP
    #pragma omp critical
    {
    #endif

    for (dspin two_k = two_k_min; two_k <= two_k_max; two_k += 2) {
    for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {

        int ii = DIV2(two_i-two_i_min);
        int ki = DIV2(two_k-two_k_min);

        matrix_get(b4, dimi, ii, ki) += sqrt(DIM(two_i) * DIM(two_k)) 
                                        * matrix_get(b4_thread, dimi, ii, ki);

    } // i
    } // k

    #ifdef USE_OMP
    }
    #endif
    
    wig_temp_free();

    #ifdef USE_OMP
    } // omp parallel
    #endif

    TENSOR_FREE(dtens);
    TENSOR_FREE(wt_ip1p2);
    TENSOR_FREE(wt_ip3p4);
    TENSOR_FREE(wt_kp1p2);
    TENSOR_FREE(wt_kp3p4);
    matrix_free(dp1);
    matrix_free(dp2);
    matrix_free(dp3);
    matrix_free(dp4);
    free(grid);
    free(qxs);
    free(wgks);
    free(wgs);

    return b4;

}
