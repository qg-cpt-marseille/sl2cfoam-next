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
#include <mpfr.h>
#include <mpc.h>

#include "utils.h"
#include "wigxjpf.h"
#include "fastwigxj.h"
#include "sl2cfoam.h"

// computes the wigner matrix D^j_(jm) (phi, theta, -phi)
// in arbitrary precision
static void wigner_dmatel(mpc_t rop, int precision, 
                          dspin two_j, int two_m, double theta, double phi) {

    int Jpm, Jmm;
    Jpm = DIV2(two_j + two_m);
    Jmm = DIV2(two_j - two_m);

    spin j = SPIN(two_j);
    spin m = SPIN(two_m);

    mpz_t rop1;
    mpz_t rop2;
    mpz_t rop3;

    mpz_init(rop1);
    mpz_init(rop2);
    mpz_init(rop3);

    mpz_fac_ui(rop1, two_j);
    mpz_fac_ui(rop2, Jpm);
    mpz_fac_ui(rop3, Jmm);

    mpz_mul(rop2, rop2, rop3);

    mpfr_t sqf;
    mpfr_init(sqf);

    // compute sqrt prefactor
    mpfr_set_z(sqf, rop1, MPFR_RNDN);
    mpfr_div_z(sqf, sqf, rop2, MPFR_RNDN);
    mpfr_sqrt(sqf, sqf, MPFR_RNDN);

    mpz_clear(rop1);
    mpz_clear(rop2);
    mpz_clear(rop3);

    mpfr_t x;
    mpfr_t c;
    mpfr_t s;
    mpfr_t pow1;
    mpfr_t pow2;

    mpfr_init2(x, precision);
    mpfr_init2(s, precision);
    mpfr_init2(c, precision);

    mpfr_set_d(x, theta * 0.5, MPFR_RNDN);

    // compute sine and cosine of theta/2
    mpfr_sin_cos(s, c, x, MPFR_RNDN);

    mpfr_init2(pow1, precision);
    mpfr_init2(pow2, precision);

    mpfr_pow_ui(pow1, c, Jpm, MPFR_RNDN);
    mpfr_pow_ui(pow2, s, Jmm, MPFR_RNDN);

    // create the D small function
    mpfr_mul(pow1, pow1, pow2, MPFR_RNDN);
    mpfr_mul(pow1, pow1, sqf, MPFR_RNDN);
    
    mpc_t mel;
    mpc_init2(mel, precision);
    mpfr_set_d(x, (j-m) * (+phi), MPC_RNDNN);

    // exp(I * x)
    mpfr_sin_cos(mel->im, mel->re, x, MPC_RNDNN);

    mpc_mul_fr(mel, mel, pow1, MPC_RNDNN);

    mpc_set(rop, mel, MPC_RNDNN);

    mpfr_clear(sqf);
    mpfr_clear(x);
    mpfr_clear(c);
    mpfr_clear(s);
    mpfr_clear(pow1);
    mpfr_clear(pow2);
    mpc_clear(mel);

}

// computes a coherent state in arbitrary precision
static double complex coherent_state(int precision, dspin two_i, 
                                     dspin two_j1, dspin two_j2, dspin two_j3, dspin two_j4,
                                     double theta_phi[4][2]) {

    size_t d1_size, d2_size, d3_size, d4_size;
    d1_size = DIM(two_j1);
    d2_size = DIM(two_j2);
    d3_size = DIM(two_j3);
    d4_size = DIM(two_j4);

    mpc_ptr wd1[d1_size];
    mpc_ptr wd2[d2_size];
    mpc_ptr wd3[d3_size];
    mpc_ptr wd4[d4_size];

    for (int i = 0; i < d1_size; i++) {
        wd1[i] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(wd1[i], precision);
    }

    for (int i = 0; i < d2_size; i++) {
        wd2[i] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(wd2[i], precision);
    }

    for (int i = 0; i < d3_size; i++) {
        wd3[i] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(wd3[i], precision);
    }

    for (int i = 0; i < d4_size; i++) {
        wd4[i] = (mpc_ptr)malloc(sizeof(mpc_t));
        mpc_init2(wd4[i], precision);
    }

    double w3j12[d1_size][d2_size];
    double w3j34[d3_size][d4_size];

    mpc_t csm;
    mpc_init2(csm, precision);
    mpc_set_ui(csm, 0, MPC_RNDNN);

    #ifdef USE_OMP
    #pragma omp parallel shared(csm) if(OMP_PARALLELIZE)
    {
    #endif

    wig_thread_temp_init(CONFIG.max_two_spin);

    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
        wigner_dmatel(wd1[DIV2(two_m1+two_j1)], precision, two_j1, two_m1, theta_phi[0][0], theta_phi[0][1]);
    }

    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
        wigner_dmatel(wd2[DIV2(two_m2+two_j2)], precision, two_j2, two_m2, theta_phi[1][0], theta_phi[1][1]);
    }

    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int two_m3 = -two_j3; two_m3 <= two_j3; two_m3 += 2) {
        wigner_dmatel(wd3[DIV2(two_m3+two_j3)], precision, two_j3, two_m3, theta_phi[2][0], theta_phi[2][1]);
    }

    #ifdef USE_OMP
    #pragma omp for
    #endif
    for (int two_m4 = -two_j4; two_m4 <= two_j4; two_m4 += 2) {
        wigner_dmatel(wd4[DIV2(two_m4+two_j4)], precision, two_j4, two_m4, theta_phi[3][0], theta_phi[3][1]);
    }

    #ifdef USE_OMP
    #pragma omp for collapse(2)
    #endif
    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
    for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
        w3j12[DIV2(two_m1+two_j1)][DIV2(two_m2+two_j2)] = 
            fw3jja6(two_j1, two_j2, two_i, two_m1, two_m2, -two_m1 - two_m2);
    }
    }

    #ifdef USE_OMP
    #pragma omp for collapse(2)
    #endif
    for (int two_m3 = -two_j3; two_m3 <= two_j3; two_m3 += 2) {
    for (int two_m4 = -two_j4; two_m4 <= two_j4; two_m4 += 2) {
        w3j34[DIV2(two_m3+two_j3)][DIV2(two_m4+two_j4)] = 
            fw3jja6(two_j3, two_j4, two_i, two_m3, two_m4, -two_m3 - two_m4);
    }
    }

    mpc_t thsum, prod;
    mpc_init2(thsum, precision);
    mpc_init2(prod, precision);

    mpc_set_ui(thsum, 0, MPC_RNDNN);

    #ifdef USE_OMP
    #pragma omp for collapse(3) schedule(dynamic, 4)
    #endif
    for (int two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2) {
    for (int two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2) {
    for (int two_m3 = -two_j3; two_m3 <= two_j3; two_m3 += 2) {

        dspin two_m4 = -two_m1-two_m2-two_m3;
        if (two_m4 < -two_j4 || two_m4 > two_j4) {
            continue;
        }

        int m1i, m2i, m3i, m4i;
        m1i = DIV2(two_m1+two_j1);
        m2i = DIV2(two_m2+two_j2);
        m3i = DIV2(two_m3+two_j3);
        m4i = DIV2(two_m4+two_j4);

        double w4j;        
        w4j = real_negpow(two_i -(-two_m1 -two_m2)) * w3j12[m1i][m2i] * w3j34[m3i][m4i];

        mpc_set_d(prod, w4j, MPC_RNDNN);
        mpc_mul(prod, prod, wd1[m1i], MPC_RNDNN);
        mpc_mul(prod, prod, wd2[m2i], MPC_RNDNN);
        mpc_mul(prod, prod, wd3[m3i], MPC_RNDNN);
        mpc_mul(prod, prod, wd4[m4i], MPC_RNDNN);
        
        mpc_add(thsum, thsum, prod, MPC_RNDNN);
                
    }
    }
    }

    #ifdef USE_OMP
    #pragma omp critical
    #endif
    mpc_add(csm, csm, thsum, MPC_RNDNN);

    mpc_clear(thsum);
    mpc_clear(prod);

    wig_temp_free();

    #ifdef USE_OMP
    } // omp parallel
    #endif

    for (int i = 0; i < d1_size; i++) {
        mpc_clear(wd1[i]);
        free(wd1[i]);
    }

    for (int i = 0; i < d2_size; i++) {
        mpc_clear(wd2[i]);
        free(wd2[i]);
    }

    for (int i = 0; i < d3_size; i++) {
        mpc_clear(wd3[i]);
        free(wd3[i]);
    }

    for (int i = 0; i < d4_size; i++) {
        mpc_clear(wd4[i]);
        free(wd4[i]);
    }

    double complex ret;
    ret = sqrt(DIM(two_i)) * mpc_get_dc(csm, MPC_RNDNN);

    mpc_clear(csm);

    return ret;

}

sl2cfoam_cvector sl2cfoam_coherentstate_range(sl2cfoam_dspin two_js[4],
                                              sl2cfoam_dspin two_i_min, sl2cfoam_dspin two_i_max,
                                              double theta_phi[4][2]) {

    // check intertwiners range
    dspin two_i_min_allowed, two_i_max_allowed;                                                               
    two_i_min_allowed = (dspin) max(abs(two_js[0]-two_js[1]), abs(two_js[2]-two_js[3]));
    two_i_max_allowed = (dspin) min(two_js[0]+two_js[1], two_js[2]+two_js[3]);

    if (two_i_min < two_i_min_allowed || two_i_max > two_i_max_allowed) {
        warning("intertwiner range must be in [%d %d]", DIV2(two_i_min_allowed), DIV2(two_i_max_allowed));
        return NULL;
    }

    // compute and put in vector
    size_t isize = DIV2(two_i_max-two_i_min) + 1;
    sl2cfoam_cvector csv = cvector_alloc(isize);

    // compute precision for mpc
    int precision = 128 + (int)(0.5 * max4(two_js[0], two_js[1], two_js[2], two_js[3]));

    for (dspin two_i = two_i_min; two_i <= two_i_max; two_i += 2) {

        double complex cs = 0. + 0.*I;

        cs = coherent_state(precision, two_i, 
                            two_js[0], two_js[1], two_js[2], two_js[3],
                            theta_phi);

        vector_set(cs, csv, DIV2(two_i-two_i_min));

    }

    return csv;

}

sl2cfoam_cvector sl2cfoam_coherentstate_fullrange(sl2cfoam_dspin two_js[4], 
                                                  double theta_phi[4][2]) {

    dspin two_i_min, two_i_max;                                                               
    two_i_min = (dspin) max(abs(two_js[0]-two_js[1]), abs(two_js[2]-two_js[3]));
    two_i_max = (dspin) min(two_js[0]+two_js[1], two_js[2]+two_js[3]);

    return sl2cfoam_coherentstate_range(two_js, two_i_min, two_i_max, theta_phi);

}
