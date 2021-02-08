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

#include "utils.h"
#include "common.h"
#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"
#include "cgamma.h"
#include "error.h"

#define check(name, cond) \
    if (!(cond)) { fprintf(stderr, "test '%s' failed\n", name); exit(EXIT_FAILURE); }

int main(int argc, char** argv) {

    dspin two_j1;
    dspin two_j2;
    dspin two_j3;
    dspin two_j4;
    dspin two_l1;
    dspin two_l2;
    dspin two_l3;
    dspin two_l4;

    double IMMIRZI = 1.2;

    bool test_acc = false;

    if (argc == 1) {

        two_j1 = 20;
        two_j2 = 20;
        two_j3 = 20;
        two_j4 = 20;
        two_l1 = 20;
        two_l2 = 20;
        two_l3 = 20;
        two_l4 = 20;
        printf("Setting j = l = 10...\n");

    } else if (argc == 10 || argc == 11) {

        two_j1 = (dspin)atoi(argv[1]);
        two_j2 = (dspin)atoi(argv[2]);
        two_j3 = (dspin)atoi(argv[3]);
        two_j4 = (dspin)atoi(argv[4]);
        two_l1 = (dspin)atoi(argv[5]);
        two_l2 = (dspin)atoi(argv[6]);
        two_l3 = (dspin)atoi(argv[7]);
        two_l4 = (dspin)atoi(argv[8]);

        sscanf(argv[9], "%lf", &IMMIRZI);

        if (argc == 11) test_acc = true;

    } else {
        error("Usage: %s [two_j1 ... two_l1 ...] [immirzi] [acc]", argv[0]);
    }

    struct sl2cfoam_config libconf;
	libconf.verbosity = SL2CFOAM_VERBOSE_HIGH;

    dspin two_l_max = max4(two_l1, two_l2, two_l3, two_l4);

    // set accuracy
    // and init library

    if (two_l_max < 2 * 50) {
        libconf.accuracy = SL2CFOAM_ACCURACY_NORMAL;
        printf("Using NORMAL accuracy ... \n");
    } else {
        libconf.accuracy = SL2CFOAM_ACCURACY_HIGH;
        printf("Using HIGH accuracy ... \n");
    }

    libconf.max_two_spin = two_l_max + 20;

    sl2cfoam_init_conf("test/test_data/", IMMIRZI, &libconf);

    // compute

    sl2cfoam_dmatrix b4 = sl2cfoam_b4(two_j1, two_j2, two_j3, two_j4,
                                      two_l1, two_l2, two_l3, two_l4);

    const int CUT = 10;

    dspin two_i_min = max(abs(two_j1-two_j2), abs(two_j3-two_j4));
    dspin two_i_max = min(two_j1+two_j2, two_j3+two_j4);
    dspin two_k_min = max(abs(two_l1-two_l2), abs(two_l3-two_l4));
    dspin two_k_max = min(two_l1+two_l2, two_l3+two_l4);

    int dimi = DIV2(two_i_max-two_i_min) + 1;

    printf("\nPrinting (at most) the first %d values...\n", (CUT+2)*(CUT+2)/4);

    for (dspin two_i = two_i_min; two_i <= min(two_i_min + CUT, two_i_max); two_i += 2) {
    for (dspin two_k = two_k_min; two_k <= min(two_k_min + CUT, two_k_max); two_k += 2) {

        double b4v = matrix_get(b4, dimi, DIV2(two_i-two_i_min), DIV2(two_k-two_k_min));

        printf("i = %.1f, k = %.1f: %.6g\n", SPIN(two_i), SPIN(two_k), b4v);

    }
    }

    printf("...\n\n");

    if (!test_acc) goto exit;

    // test the first b4s using accurate integration

    printf("Testing the non-zero values shown against accurate integration...\n");
    double eps = 1e-6;

    int dimia = DIV2(min(two_i_min + CUT, two_i_max)-two_i_min) + 1;

    sl2cfoam_dmatrix b4_acc = sl2cfoam_b4_accurate(two_j1, two_j2, two_j3, two_j4,
                                                   two_l1, two_l2, two_l3, two_l4,
                                                   two_i_min, min(two_i_min + CUT, two_i_max),
                                                   two_k_min, min(two_k_min + CUT, two_k_max));

    for (dspin two_i = two_i_min; two_i <= min(two_i_min + CUT, two_i_max); two_i += 2) {
    for (dspin two_k = two_k_min; two_k <= min(two_k_min + CUT, two_k_max); two_k += 2) {

        double b4f = matrix_get(b4, dimi, DIV2(two_i-two_i_min), DIV2(two_k-two_k_min));
        if (b4f == 0.0) continue;

        double b4a = matrix_get(b4_acc, dimia, DIV2(two_i-two_i_min), DIV2(two_k-two_k_min));

        printf("i = %.1f, k = %.1f: b4 fast = %.12g, b4 accurate = %.12g, relative error = %g\n", 
               SPIN(two_i), SPIN(two_k), b4f, b4a, fabs((b4a-b4f)/b4a));

    }
    }

exit:

    sl2cfoam_free();

    return EXIT_SUCCESS;

}