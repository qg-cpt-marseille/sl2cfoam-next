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

    const char* usage = "Usage: %s [options] folder Immirzi 2j1,...,2j10 Dl\n";

    if (argc < 4) {
        fprintf(stderr, usage, argv[0]);
        exit(EXIT_FAILURE);
    }

    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    #endif

    MPI_FUNC_INIT();

    /////////////////////////////////////////////////////////////////////
    // read command line args
    /////////////////////////////////////////////////////////////////////

    int verbosity = SL2CFOAM_VERBOSE_OFF;
    int accuracy = SL2CFOAM_ACCURACY_NORMAL;

    dspin max_two_spin = 0;
    size_t MB_per_thread = 0;

    bool store_batches = false;

    int opt;
    while ((opt = getopt(argc, argv, "vVhHBJ:m:")) != -1) {
        switch (opt) {

        case 'v': verbosity = SL2CFOAM_VERBOSE_LOW; break;
        case 'V': verbosity = SL2CFOAM_VERBOSE_HIGH; break;

        case 'h': accuracy = SL2CFOAM_ACCURACY_HIGH; break;
        case 'H': accuracy = SL2CFOAM_ACCURACY_VERYHIGH; break;

        case 'B': store_batches = true; break;

        case 'J': 
            sscanf(optarg, "%d", &max_two_spin);
            break;

        case 'm':
            sscanf(optarg, "%lu", &MB_per_thread);
			break;

        default:
            fprintf(stderr, usage, argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    int argind = optind;

    // parse folder
    char folder[1024];
    strcpy(folder, argv[argind++]);

    // parse Immirzi
    double Immirzi;
    sscanf(argv[argind++], "%lf", &Immirzi);

    // parse spins
    dspin two_js[10];
    sscanf(argv[argind++], "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", 
                           &two_js[0], &two_js[1], &two_js[2], &two_js[3], &two_js[4],
		 		           &two_js[5], &two_js[6], &two_js[7], &two_js[8], &two_js[9]);

    // parse Dl
    int Dl;
    sscanf(argv[argind++], "%d", &Dl);
    dspin two_Dl = 2 * Dl;

    /////////////////////////////////////////////////////////////////////
    // initialize the library
    /////////////////////////////////////////////////////////////////////

    struct sl2cfoam_config libconf;
	libconf.verbosity = verbosity;
    libconf.accuracy = accuracy;
    libconf.max_MB_mem_per_thread = MB_per_thread;

    if (max_two_spin == 0) {

        // find maximum spin
        dspin two_j_max = 0;

        for (int i = 0; i < 10; i++) {
            if (two_js[i] > two_j_max) two_j_max = two_js[i];
        }

        // TODO: following criterion is not tested
        libconf.max_two_spin = 3 * (two_j_max + two_Dl);

    } else {
        libconf.max_two_spin = max_two_spin;
    }

	sl2cfoam_init_conf(folder, Immirzi, &libconf);

    sl2cfoam_tensor_result tresult = TENSOR_RESULT_STORE;
    if (store_batches) tresult |= TENSOR_RESULT_STORE_BATCHES;

    sl2cfoam_vertex_fullrange(two_js, Dl, tresult);

    sl2cfoam_free();

    #ifdef USE_MPI
    MPI_Finalize();
    #endif 

}