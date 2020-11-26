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

#include "common.h"
#include "utils.h"
#include "error.h"

// Conventions: the map to Speziale is
//
//  C  |  S
// ---------
//  J  |  k
// rho | rho
//  j' |  l
//  j  |  j
//  p  |  p
//
// j2 is j', j1 is j.

static inline void alpha(mpc_t a, int prec, double rho, 
                         dspin two_J, dspin two_j2, dspin two_j1, dspin two_p,
                         int s2, int s1, int m) {

    double j1mm = SPIN(two_j1) - m;
    int Jmp = DIV2(abs(two_J-two_p)); // |J-p|                                         

    mpc_t prod;
    mpc_init2(prod, prec);
    mpc_set_d_d(prod, j1mm, -rho, MPC_RNDNN);

    mpc_t pc;
    mpc_init2(pc, prec);

    // complex product
    int k_max = Jmp + s1 + s2;
    for (int k = 1; k <= k_max; k++) {
        mpc_set_d_d(pc, j1mm - k, -rho, MPC_RNDNN);
        mpc_mul(prod, prod, pc, MPC_RNDNN);
    }

    // multiply binomial
    mpz_t bz;
    mpz_init(bz);

    mpfr_t bf;
    mpfr_init2(bf, prec);

    mpz_set_si(bz, DIV2(two_j1 + two_j2) - Jmp - s1 - s2);
    mpz_bin_ui(bz, bz, m);
    mpfr_set_z(bf, bz, MPFR_RNDN);

    mpc_fr_div(a, bf, prod, MPC_RNDNN);

    // clear
    mpc_clear(prod);
    mpc_clear(pc);
    mpfr_clear(bf);
    mpz_clear(bz);

}

void sl2cfoam_dsmall_Yc(mpc_ptr Ys[], int prec, double rho, 
                        dspin two_k, dspin two_j, dspin two_l, dspin two_p) {

    // map names to Collet conventions
    dspin two_J = two_k;
    dspin two_j1 = two_j;
    dspin two_j2 = two_l;


    int Jmp = DIV2(abs(two_J-two_p)); // |J-p|
    int Jpp = DIV2(abs(two_J+two_p)); // |J+p|

    // maximum index m
    int m_max = DIV2(two_j1 + two_j2) - Jmp;

    int s1_max = DIV2(two_j1 - Jmp - Jpp);
    int s2_max = DIV2(two_j2 - Jmp - Jpp);

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE)
    {
    #endif

    mpc_t sum;
    mpc_init2(sum, prec);

    mpc_t psum;
    mpc_init2(psum, prec);

    mpfr_t af;
    mpfr_init2(af, prec);

    mpz_t az, bz;
    mpz_inits(az, bz, NULL);

    int deltam_down, deltam_up;
    int B;

    // loop over index m
    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int m = 0; m <= m_max; m++) {

        mpc_set_ui(sum, 0, MPC_RNDNN);

        for (int s1 = 0; s1 <= s1_max; s1++) {

            deltam_up = DIV2(two_j1 + two_j2) - Jmp - s1;
            if (m > deltam_up) continue;

        for (int s2 = 0; s2 <= s2_max; s2++) {

            deltam_down = s2;
            if (m < deltam_down) continue;

            // factors f and N (with binomials)
            B = Jmp + s1 + s2;

            mpz_set_si(bz, B);
            mpz_bin_ui(az, bz, s1);
            mpfr_set_z(af, az, MPFR_RNDN);

            mpz_bin_ui(az, bz, s2);
            mpfr_mul_z(af, af, az, MPFR_RNDN);

            mpz_fac_ui(bz, B);

            mpz_fac_ui(az, DIV2(two_j1 + Jpp - Jmp) - s1);
            mpz_mul(bz, bz, az);

            mpz_fac_ui(az, DIV2(two_j1 - Jpp - Jmp) - s1);
            mpz_mul(bz, bz, az);

            mpz_fac_ui(az, DIV2(two_j2 + Jpp - Jmp) - s2);
            mpz_mul(bz, bz, az);

            mpz_fac_ui(az, DIV2(two_j2 - Jpp - Jmp) - s2);
            mpz_mul(bz, bz, az);

            mpfr_div_z(af, af, bz, MPFR_RNDN);

            // multiply alpha and sum

            alpha(psum, prec, rho,
                  two_J, two_j2, two_j1, two_p, s2, s1, m-s2);
            mpc_mul_si(psum, psum, real_negpow(2 * (m - s1)), MPC_RNDNN);
            mpc_mul_fr(psum, psum, af, MPC_RNDNN);

            mpc_add(sum, sum, psum, MPC_RNDNN);

        }
        }

        // add normalization prefactor

        mpz_set_si(bz, DIM(two_j1) * DIM(two_j2));

        mpz_fac_ui(az, DIV2(two_j1 + two_J));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j1 - two_J));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j1 + two_p));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j1 - two_p));
        mpz_mul(bz, bz, az);

        mpz_fac_ui(az, DIV2(two_j2 + two_J));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j2 - two_J));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j2 + two_p));
        mpz_mul(bz, bz, az);
        mpz_fac_ui(az, DIV2(two_j2 - two_p));
        mpz_mul(bz, bz, az);

        mpfr_set_z(af, bz, MPFR_RNDN);
        mpfr_sqrt(af, af, MPFR_RNDN);

        mpc_mul_fr(Ys[m], sum, af, MPC_RNDNN);

    }

    // clear
    mpc_clear(sum);
    mpc_clear(psum);
    mpfr_clear(af);
    mpz_clears(az, bz, NULL);

    #ifdef USE_OMP
    } // omp parallel
    #endif

}

void sl2cfoam_dsmall(__complex128 ds[], __float128 xs[], size_t N,
                     int prec, mpc_ptr Ym[], mpc_ptr Yn[], mpc_ptr prefactor[],
                     double rho, dspin two_k, dspin two_j, dspin two_l, dspin two_p) {

    if (rho == 0) {
        error("dsmall with rho == 0 not implemented");
    }

    // map names to Collet conventions
    dspin two_J = two_k;
    dspin two_j1 = two_j;
    dspin two_j2 = two_l;

    int Jmp = DIV2(abs(two_J-two_p)); // |J-p|
    int Jpp = DIV2(abs(two_J+two_p)); // |J+p|

    int m_max = DIV2(two_j1 + two_j2) - Jmp;
    int n_max = DIV2(two_j1 + two_j2) - Jpp;

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE)
    {
    #endif

    mpc_t emi_r_rho_ab;
    mpc_t sum_m;
    mpc_t sum_n;
    mpc_t sum;
    mpc_t psum;
    mpc_t pref;

    mpfr_t emr_ab;
    mpfr_t mr_rho_ab;
    mpfr_t emr_sq_ab;
    mpfr_t emul;

    mpc_init2(sum_m, prec);
    mpc_init2(sum_n, prec);
    mpc_init2(sum, prec);
    mpc_init2(emi_r_rho_ab, prec);
    mpc_init2(psum, prec);
    mpc_init2(pref, prec);

    mpfr_init2(emr_ab, prec);
    mpfr_init2(mr_rho_ab, prec);
    mpfr_init2(emr_sq_ab, prec);
    mpfr_init2(emul, prec);

    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int i = 0; i < N; i++) {

        // exp(-r)
        mpfr_set_float128(emr_ab, xs[i], MPFR_RNDN);

        // exp(-2r)
        mpfr_sqr(emr_sq_ab, emr_ab, MPFR_RNDN);

        if (prefactor == NULL) {

            // get exp(-i*r*rho)
            mpfr_log(mr_rho_ab, emr_ab, MPFR_RNDN);
            mpfr_mul_d(mr_rho_ab, mr_rho_ab, rho, MPFR_RNDN);
            mpfr_sin_cos(emi_r_rho_ab->im, emi_r_rho_ab->re, mr_rho_ab, MPFR_RNDN);

            // prefactor
            mpc_set_fr(pref, emr_sq_ab, MPFR_RNDN);
            mpfr_ui_sub(pref->re, 1, pref->re, MPFR_RNDN);
            mpfr_pow_si(pref->re, pref->re, -DIV2(two_j1+two_j2) - 1, MPFR_RNDN);

            mpc_mul(pref, pref, emi_r_rho_ab, MPC_RNDNN);

        } else {
            mpc_set(pref, prefactor[i], MPC_RNDNN);
        }

        // sum over m
        mpc_set_ui(sum_m, 0, MPC_RNDNN);
        mpfr_pow_ui(emul, emr_ab, 1 + Jmp, MPFR_RNDN); // emul = exp(-r*(1 + |J-p|))
        for (int m = 0; m <= m_max; m++) {

            mpfr_set_fr(psum->re, emul, MPFR_RNDN);
            mpfr_set_ui(psum->im, 0, MPFR_RNDN);
            mpc_mul(psum, psum, Ym[m], MPC_RNDNN);
            mpc_add(sum_m, sum_m, psum, MPC_RNDNN);
            mpfr_mul(emul, emul, emr_sq_ab, MPFR_RNDN); // emul *= exp(-2r)

        }
        mpc_mul(sum_m, sum_m, pref, MPC_RNDNN);

        // sum over n
        mpc_set_ui(sum_n, 0, MPC_RNDNN);
        mpfr_pow_ui(emul, emr_ab, 1 + Jpp, MPFR_RNDN);
        for (int n = 0; n <= n_max; n++) {

            mpfr_set_fr(psum->re, emul, MPFR_RNDN);
            mpfr_set_ui(psum->im, 0, MPFR_RNDN);
            mpc_mul(psum, psum, Yn[n], MPC_RNDNN);
            mpc_add(sum_n, sum_n, psum, MPC_RNDNN);
            mpfr_mul(emul, emul, emr_sq_ab, MPFR_RNDN);

        }
        mpc_conj(pref, pref, MPC_RNDNN);
        mpc_mul(sum_n, sum_n, pref, MPC_RNDNN);
        mpc_mul_si(sum_n, sum_n, real_negpow(two_j2 - two_j1), MPC_RNDNN);

        mpc_add(sum, sum_m, sum_n, MPC_RNDNN);

        ds[i] =       mpfr_get_float128(sum->re, MPFR_RNDN)
                + I * mpfr_get_float128(sum->im, MPFR_RNDN);

    }

    mpc_clear(emi_r_rho_ab);
    mpc_clear(sum_m);
    mpc_clear(sum_n);
    mpc_clear(sum);
    mpc_clear(psum);
    mpc_clear(pref);

    mpfr_clear(mr_rho_ab);
    mpfr_clear(emr_ab);
    mpfr_clear(emr_sq_ab);
    mpfr_clear(emul);

    #ifdef USE_OMP
    } // omp parallel
    #endif

}