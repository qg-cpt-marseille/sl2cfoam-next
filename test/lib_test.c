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

#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"

#include "common.h"
#include "utils.h"
#include "mpi_utils.h"

#define check(name, cond) \
    if (!(cond)) { fprintf(stderr, "test '%s' failed\n", name); passed = false; goto clean; }

#define check_within_eps(name, v, vref) \
    check(name, fabs((v-vref)/vref) < 1e-6)

#define NTEST 3

int main(int argc, char **argv) {

    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    #endif

    MPI_FUNC_INIT();

    // parameters
    double Immirzi = 1.2;
    dspin two_j = 2;
    int Dl = 1;

    dspin two_js[10] = { two_j, two_j, two_j, two_j, two_j,
                         two_j, two_j, two_j, two_j, two_j };
    
    tensor_ptr(vertex) t;

    // test some values
    int is[NTEST][5] = { {0, 0, 0, 2, 0}, {0, 0, 1, 0, 1}, {2, 0, 1, 2, 1} };

    // values depend on the Y-map

    #ifdef RHO_GJ
    double ampls[NTEST] = { 6.17259170394e-09, 5.33694074629e-09, -7.26116697167e-10 };
    #endif

    #ifdef RHO_GJP1
    double ampls[NTEST] = { 8.76263107783e-11, 7.19679864345e-11, -1.0176125168e-11 };
    #endif
    
    // library conf
    struct sl2cfoam_config libconf;
	libconf.verbosity = SL2CFOAM_VERBOSE_HIGH;
    libconf.accuracy = SL2CFOAM_ACCURACY_NORMAL;
    libconf.max_two_spin = 20;
    libconf.max_MB_mem_per_thread = 0;

    ///////////////////////////////////////////////////////

    sl2cfoam_init_conf("test/test_data/", Immirzi, &libconf);

    MPI_MASTERONLY_DO printf("Computing full range tensor all 2j = %d, Dl = %d...\n", two_j, Dl);

    t = sl2cfoam_vertex_fullrange(two_js, Dl, TENSOR_RESULT_RETURN);

MPI_MASTERONLY_START

    printf("Testing values...\n");

    bool passed = true;

    double a;
    char msg[64];
    for (int s = 0; s < NTEST; s++) {
        a = TENSOR_GET(t, 5, is[s][4], is[s][3], is[s][2], is[s][1], is[s][0]);
        sprintf(msg, "eq %d", s+1);
        check_within_eps(msg, a, ampls[s]);
    }

clean:

    TENSOR_FREE(t);
    
    // remove boost tensors
    remove("test/test_data/vertex/immirzi_1.200/boosters/b4__2-2-2-2__gf-1__imm-1.200__dl-1.sl2t");
    remove("test/test_data/vertex/immirzi_1.200/boosters/b4__2-2-2-2__gf-2__imm-1.200__dl-1.sl2t");
    remove("test/test_data/vertex/immirzi_1.200/boosters/b4__2-2-2-2__gf-3__imm-1.200__dl-1.sl2t");
    remove("test/test_data/vertex/immirzi_1.200/boosters/b4__2-2-2-2__gf-4__imm-1.200__dl-1.sl2t");

    // remove directories
    remove("test/test_data/vertex/immirzi_1.200/amplitudes");
    remove("test/test_data/vertex/immirzi_1.200/boosters");
    remove("test/test_data/vertex/immirzi_1.200");
    remove("test/test_data/vertex");

    if (passed) {
        printf("... test PASSED.\n");
    } else {
        printf("... test FAILED.\n");
    }

MPI_MASTERONLY_END

    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    
}