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

///////////////////////////////////////////////////////////////
// Setup functions.
// Initialize various configurataions at start.
///////////////////////////////////////////////////////////////

#include <dirent.h> 
#include <string.h>
#include <sys/stat.h>
#include <mpfr.h>
#include <omp.h>

#include "common.h"
#include "sl2cfoam.h"
#include "wigxjpf.h"
#include "fastwigxj.h"
#include "utils.h"
#include "error.h"
#include "mpi_utils.h"

// path of found fastwig tables
char* FASTWIG_3J_TABLE_PATH;
char* FASTWIG_6J_TABLE_PATH;

#ifdef USE_MPI
// stores if MPI was found to be initialized outside the library
static bool __mpi_managed_outside;
#endif

// searches for all files in root folder 
// with extension ".3j" or ".6j". 
// The largest tables found are set for loading.
static void find_fastwig_tables(char* folder) {
  
    DIR *d;
    struct dirent *dir;
    d = opendir(folder);

    if (d == NULL) error("error opening root folder");

    char table_3j_path[strlen(folder) + 256];
    char table_3j_path_toload[strlen(folder) + 256];
    long table_3j_size = 0L;
    
    char table_6j_path[strlen(folder) + 256];
    char table_6j_path_toload[strlen(folder) + 256];
    long table_6j_size = 0L;

    char* fname;
    while ((dir = readdir(d)) != NULL) {

        fname = dir->d_name;
        char *dot = strrchr(fname, '.');

        // 3j
        if (dot && !strcmp(dot, ".3j")) {

            strcpy(table_3j_path, folder);
            strcat(table_3j_path, "/");
            strcat(table_3j_path, fname);

            // check size
            struct stat st;
            stat(table_3j_path, &st);

            if (st.st_size > table_3j_size) {

                // (larger) found, set path and size
                strcpy(table_3j_path_toload, folder);
                strcat(table_3j_path_toload, "/");
                strcat(table_3j_path_toload, fname);

                table_3j_size = st.st_size;

            }

        }

        // 6j
        if (dot && !strcmp(dot, ".6j")) {

            strcpy(table_6j_path, folder);
            strcat(table_6j_path, "/");
            strcat(table_6j_path, fname);

            // check size
            struct stat st;
            stat(table_6j_path, &st);

            if (st.st_size > table_6j_size) {

                // (larger) found, set path and size
                strcpy(table_6j_path_toload, folder);
                strcat(table_6j_path_toload, "/");
                strcat(table_6j_path_toload, fname);

                table_6j_size = st.st_size;

            }

        }

    }
    closedir(d);

    if (table_3j_size == 0) error("no 3j table found in root folder");
    if (table_6j_size == 0) error("no 6j table found in root folder");

    FASTWIG_3J_TABLE_PATH = strdup(table_3j_path_toload);
    FASTWIG_6J_TABLE_PATH = strdup(table_6j_path_toload);

}

static inline void load_fastwig_tables() {
    fastwigxj_load(FASTWIG_3J_TABLE_PATH, 3, NULL);
    fastwigxj_load(FASTWIG_6J_TABLE_PATH, 6, NULL);
}

static inline void unload_fastwig_tables() {
    fastwigxj_unload(3);
    fastwigxj_unload(6);
}

void sl2cfoam_init_conf(char* root_folder, double Immirzi, struct sl2cfoam_config* conf) {

    // check root folder is accessible
    DIR *d = opendir(root_folder);
    if (d == NULL) error("error opening root folder");
    closedir(d);

    DATA_ROOT = strdup(root_folder);

    memcpy(&CONFIG, conf, sizeof(struct sl2cfoam_config));

    // set Immirzi and create folder structure
    sl2cfoam_set_Immirzi(Immirzi);

    sl2cfoam_set_verbosity(conf->verbosity);
    sl2cfoam_set_accuracy(conf->accuracy);
    
    // initialize wigxjpf
    wig_table_init(conf->max_two_spin, 6);

    #ifndef NO_IO

    // load fastwig tables
    find_fastwig_tables(DATA_ROOT);
    load_fastwig_tables();

    #endif

    // enable OMP parallelization by default
    OMP_PARALLELIZE = true;

    // no nested parallelism
    omp_set_max_active_levels(1);

    // setup BLAS libraries
    #ifdef USE_MKL

    // 64bit interface
    mkl_set_interface_layer(MKL_INTERFACE_ILP64);

    // disable builtin threading
    mkl_set_threading_layer(MKL_THREADING_SEQUENTIAL);

    #elif USE_SYSBLAS

    #ifdef OPENBLAS_THREAD
    openblas_set_num_threads(1)
    #endif

    #endif

    // if library is compiled for MPI
    // initialize MPI if it's not already initialized
    // needed for using as shared library
    // runs on a single node if not called with mpirun
    #ifdef USE_MPI

    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(NULL, NULL);
        __mpi_managed_outside = false;
    } else {
        __mpi_managed_outside = true;
    }

    #endif

}

void sl2cfoam_init(char* root_folder, double Immirzi) {
	
    struct sl2cfoam_config def;
    def.verbosity = SL2CFOAM_VERBOSE_OFF;
    def.accuracy = SL2CFOAM_ACCURACY_NORMAL;
    def.max_two_spin = 3 * 2 * 50;

    sl2cfoam_init_conf(root_folder, Immirzi, &def);

}

void sl2cfoam_free() {
	
    // wigxjpf
    wig_table_free();

    #ifndef NO_IO
    unload_fastwig_tables();
    #endif

    // free paths
    free(DATA_ROOT);
    free(DIR_BOOSTERS);
    free(DIR_AMPLS);
    free(FASTWIG_6J_TABLE_PATH);
    free(FASTWIG_3J_TABLE_PATH);

    #ifdef USE_MPI

    if (!__mpi_managed_outside) {
        MPI_Finalize();
    }
    
    #endif

}
