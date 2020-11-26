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
#include <math.h>

#include "common.h"
#include "dsmall.h"
#include "integration_gk.h"
#include "integration_qagp.h"

#define check(name, cond) \
    if (!(cond)) { fprintf(stderr, "test '%s' failed\n", name); exit(EXIT_FAILURE); }

void f_qagp(__complex128 ys[], __float128 xs[], size_t N, void* params) {

    __float128 x, y;
    for (int i = 0; i < N; i++) {

        x = xs[i];
        y = x * cosq(1/x);
        ys[i] = y + I * 0.0Q;

    }

}

void f_gk(double ys[], __float128 xs[], size_t N) {

    double x, y;
    for (int i = 0; i < N; i++) {

        x = (double)xs[i];
        y = x * cos(1/x);
        ys[i] = y;

    }

}

int main(int argc, char** argv) {

    // test QAGP integration

    sl2cfoam_integration_function F2;
    F2.function = &f_qagp;
    F2.measure = NULL;
    F2.params = NULL;

    __complex128 r2;
    __float128 ae2;
    size_t neval;
    sl2cfoam_qagp_retcode rc;
    rc = sl2cfoam_qagp(&F2, 1e-9, &r2, &ae2, &neval, 1e-9, 5000);

    printf("integral (QAGP) = %.12g, abserr = %.3g, neval = %d, ret = %d\n", (double)crealq(r2), (double)ae2, (int)neval, (int)rc);

    // test fixed HARMONIC integration

    int intervals = 10;
    __float128* grid_harmonic = sl2cfoam_grid_harmonic(intervals);
    __float128* grid_uniform = sl2cfoam_grid_uniform(intervals);
    int nxs = intervals * GK_POINTS;

    double ms[nxs];
    for (int i = 0; i < nxs; i++) {
        ms[i] = 1.0;
    }

    // compute weights
    double* wgks_harmonic = sl2cfoam_gk_grid_weights(intervals, grid_harmonic);
    double* wgs_harmonic = sl2cfoam_gk_grid_weights_gauss(intervals, grid_harmonic);
    double* wgks_uniform = sl2cfoam_gk_grid_weights(intervals, grid_uniform);
    double* wgs_uniform = sl2cfoam_gk_grid_weights_gauss(intervals, grid_uniform);

    __float128* xs_harmonic = sl2cfoam_gk_grid_abscissae(intervals, grid_harmonic);
    __float128* xs_uniform = sl2cfoam_gk_grid_abscissae(intervals, grid_uniform);

    double ysh[nxs];
    double ysu[nxs];

    double ae3;
    double r3;

    f_gk(ysh, xs_harmonic, nxs);
    r3 = sl2cfoam_gk_grid(intervals, ysh, ms, wgks_harmonic, wgs_harmonic, &ae3);
    printf("integral (HARM) = %.12g, abserr = %.3g, nevals = %d\n", r3, ae3, nxs);

    f_gk(ysu, xs_uniform, nxs);
    r3 = sl2cfoam_gk_grid(intervals, ysu, ms, wgks_uniform, wgs_uniform, &ae3);
    printf("integral (UNIF) = %.12g, abserr = %.3g, nevals = %d\n", r3, ae3, nxs);

    free(grid_harmonic);
    free(grid_uniform);
    free(xs_harmonic);
    free(xs_uniform);
    free(wgs_harmonic);
    free(wgks_harmonic);
    free(wgs_uniform);
    free(wgks_uniform);

    return EXIT_SUCCESS;

}