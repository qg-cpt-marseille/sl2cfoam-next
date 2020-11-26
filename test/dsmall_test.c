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

#include <quadmath.h>
#include <complex.h>

#include "utils.h"
#include "common.h"
#include "dsmall.h"
#include "integration_qagp.h"

#define check(name, cond) \
    if (!(cond)) { fprintf(stderr, "test '%s' failed\n", name); exit(EXIT_FAILURE); }

typedef struct __dsmall_params {
    int prec;
    double rho;
    dspin two_ji;
    dspin two_li;
    dspin two_pi;
    mpc_ptr* Ym;
    mpc_ptr* Yn;
} __dsmall_params;

static void dsmall_qagp(__complex128 ds[], __float128 xs[], size_t N, void* params) {

    __dsmall_params* p = (__dsmall_params*)params;
    sl2cfoam_dsmall(ds, xs, N, p->prec, p->Ym, p->Yn, NULL,
                    p->rho, p->two_ji, p->two_ji, p->two_li, p->two_pi);

}

int main(int argc, char** argv) {

    dspin two_ji = (dspin)atoi(argv[1]);
    dspin two_li = (dspin)atoi(argv[2]);
    dspin two_pi = (dspin)atoi(argv[3]);

    double IMMIRZI;
    sscanf(argv[4], "%lf", &IMMIRZI);

    dspin two_j_max = two_li;

    // fix precision
    // strategy: fix the minimum allowed interval for the integration
    //           then bit precision follows
    __float128 dx_min = 1e-6;

    int prec_base, prec_add;

    // there should be enough significant digits to cancel the prefactor 
    // (1-(1-dx)^2)^-(j1+j2+1) ~ (2*dx)^-(j1+j2+1) ...
    prec_add = - (int)((DIV2(two_ji+two_li)+1) * log2q(2 * dx_min));

    // and some more of course
    // base precision is found by "trial and error"
    if (two_j_max <= 50) {
        prec_base = 64;
    } else if (two_j_max <= 100) {
        prec_base = 128;
    } else {
        prec_base = 256;
    }

    // total precision;
    int precision = prec_base + prec_add;
    printf("Working with precision of %d bits...\n", precision);

    spin ji = SPIN(two_ji);
    double rho = RHO(ji);

    // compute the Y coefficients
    int Jmp = DIV2(abs(two_ji - two_pi)); // |J-p|
    int Jpp = DIV2(abs(two_ji + two_pi)); // |J+p|
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

    // NOTICE the order of the arguments
    sl2cfoam_dsmall_Yc(Ym, precision,  rho, two_ji, two_li, two_ji,  two_pi);
    sl2cfoam_dsmall_Yc(Yn, precision, -rho, two_ji, two_ji, two_li, -two_pi);

    // compute dsmall if argument N (number of points) specified
    // prints points for importing in Julia
    if (argc > 5) {

        int N = atoi(argv[5]);

        __float128 x_max = 1.0Q;
        __float128 dx = x_max / N;
        __float128 x[N-1]; // exp(-r)

        for (int i = 1; i < N; i++) {
            x[i-1] = i * dx;
        }

        __complex128 d[N-1];
        sl2cfoam_dsmall(d, x, N-1, precision, Ym, Yn, NULL, rho, two_ji, two_ji, two_li, two_pi);

        // print for plotting in Julia
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

        // compute integral
        // (simple trapezoidal)
        __complex128 sum = 0.0Q;
        double sumr, sumi;

        //sum += 0.5 * d[0];
        for (int i = 1; i < N-1; i++) {
            sum += d[i];
        }
        sum += 0.5 * d[N-1];
        sum *= dx;

        sumr = creal((double complex)sum);
        sumi = cimag((double complex)sum);
        printf("integral (trapz) = %.6g + I * %.6g\n", sumr, sumi);

    }

    // integrate
    __dsmall_params p;
    p.prec = precision;
    p.rho = rho;
    p.two_ji = two_ji;
    p.two_li = two_li;
    p.two_pi = two_pi;
    p.Ym = Ym;
    p.Yn = Yn;

    printf("Testing adaptive QAGP integration (dx_min = %.6g)...\n", (double)dx_min);

    sl2cfoam_integration_function f_qagp;
    f_qagp.function = &dsmall_qagp;
    f_qagp.measure = NULL;
    f_qagp.params = &p;

    __complex128 res_qagp;
    __float128 abserr_qagp;
    size_t nevals_qagp;
    sl2cfoam_qagp_retcode rcode;
    rcode = sl2cfoam_qagp(&f_qagp, 1e-6, &res_qagp, &abserr_qagp, &nevals_qagp, dx_min, 4000);

    printf("integral (QAGP) = %.6g + I * %.6g (err: %.6g | nevals: %d | retcode = %d)\n", 
            (double)crealq(res_qagp), (double)cimagq(res_qagp), (double)abserr_qagp, (int)nevals_qagp, rcode);

    for (int m = 0; m <= m_max; m++) {
        mpc_clear(Ym[m]);
        free(Ym[m]);
    }

    for (int n = 0; n <= n_max; n++) {
        mpc_clear(Yn[n]);
        free(Yn[n]);
    }

    return EXIT_SUCCESS;

}