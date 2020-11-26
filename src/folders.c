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

#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#include "common.h"
#include "folders.h"
#include "error.h"
#include "utils.h"
#include "mpi_utils.h"

// assembles the paths for all needed folders
static void build_folders(struct sl2cfoam_data_folders* df) {

    int len = strlen(DATA_ROOT) + 256;

    df->vertex = (char*)malloc(len*sizeof(char));

    df->vertex_imm = (char*)malloc(len*sizeof(char));
    df->vertex_imm_boost = (char*)malloc(len*sizeof(char));
    df->vertex_imm_ampl = (char*)malloc(len*sizeof(char));
    
    strcpy(df->vertex, DATA_ROOT);
    strcpy(df->vertex_imm, DATA_ROOT);
    strcpy(df->vertex_imm_boost, DATA_ROOT);
    strcpy(df->vertex_imm_ampl, DATA_ROOT);

    char tmp[len];

    strcat(df->vertex, "/vertex/");

    sprintf(tmp, "/vertex/immirzi_%.3f/", IMMIRZI);
    strcat(df->vertex_imm, tmp);

    sprintf(tmp, "/vertex/immirzi_%.3f/boosters/", IMMIRZI);
    strcat(df->vertex_imm_boost, tmp);

    sprintf(tmp, "/vertex/immirzi_%.3f/amplitudes/", IMMIRZI);
    strcat(df->vertex_imm_ampl, tmp);

}

struct sl2cfoam_data_folders* sl2cfoam_get_data_folders() {

    struct sl2cfoam_data_folders* df;
    df = (struct sl2cfoam_data_folders*) malloc(sizeof(struct sl2cfoam_data_folders));

    // build
    build_folders(df);

    return df;

}

void sl2cfoam_free_data_folders(struct sl2cfoam_data_folders* df) {

    free(df->vertex);
    free(df->vertex_imm);
    free(df->vertex_imm_boost);
    free(df->vertex_imm_ampl);
    free(df);

}

void sl2cfoam_check_data_folders(struct sl2cfoam_data_folders* df) {

    MPI_FUNC_INIT();

MPI_MASTERONLY_START

    if (mkdir(DATA_ROOT, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create root directory"); }
    }

    if (mkdir(df->vertex, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(df->vertex_imm, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(df->vertex_imm_boost, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

    if (mkdir(df->vertex_imm_ampl, 0755) == -1) {
        if (errno != EEXIST) { error("cannot create directory"); }
    }

MPI_MASTERONLY_END

}
