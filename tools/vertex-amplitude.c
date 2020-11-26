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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"
#include "utils.h"
#include "mpi_utils.h"

typedef sl2cfoam_dspin dspin;
typedef sl2cfoam_spin spin;

int main(int argc, char **argv) {

    if (argc < 4) {
        fprintf(stderr, "Usage: %s folder Immirzi two_j1,two_j2,...,two_j10 two_i1,...,two_i5 Dl\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    #endif

    MPI_FUNC_INIT();

    // parse folder
    char folder[1024];
    strcpy(folder, argv[1]);

    // parse Immirzi
    double Immirzi;
    sscanf(argv[2], "%lf", &Immirzi);

    // parse spins
    dspin two_js[10];
    sscanf(argv[3], "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", 
                    &two_js[0], &two_js[1], &two_js[2], &two_js[3], &two_js[4],
		 		    &two_js[5], &two_js[6], &two_js[7], &two_js[8], &two_js[9]);

    // parse intertwiners
    dspin two_is[5];
    sscanf(argv[4], "%d,%d,%d,%d,%d", 
                    &two_is[0], &two_is[1], &two_is[2], &two_is[3], &two_is[4]);

    // parse Dl
    int Dl;
    sscanf(argv[5], "%d", &Dl);
    dspin two_Dl = 2 * Dl;

    // initialize the library
    struct sl2cfoam_config libconf;
	libconf.verbosity = SL2CFOAM_VERBOSE_LOW;

    // find maximum spin
    dspin two_j_max = 0;
    for (int i = 0; i < 10; i++) {
        if (two_js[i] > two_j_max) two_j_max = two_js[i];
    }

    // TODO: following criterion is not tested
    libconf.max_two_spin = 3 * (two_j_max + two_Dl);

    if (two_j_max + two_Dl < 2 * 50) {
        libconf.accuracy = SL2CFOAM_ACCURACY_NORMAL;
    } else {
        libconf.accuracy = SL2CFOAM_ACCURACY_HIGH;
    }

    libconf.max_MB_mem_per_thread = 2000;

	sl2cfoam_init_conf(folder, Immirzi, &libconf);

    MPI_MASTERONLY_DO printf("Computing a single amplitude...\n");

    double a = sl2cfoam_vertex_amplitude(two_js, two_is, Dl);

    MPI_MASTERONLY_DO printf("done. Result = %.12g\n", a);

    sl2cfoam_free();

    #ifdef USE_MPI
    MPI_Finalize();
    #endif

}