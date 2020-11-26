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

#include <omp.h>

#include "common.h"
#include "folders.h"
#include "error.h"
#include "utils.h"
#include "mpi_utils.h"
#include "boosters.h"
#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"

// default filename for booster tensors
static const char* boosters_fn = "b4__%d-%d-%d-%d__gf-%d__imm-%.3f__dl-%d.sl2t";

// maximum allowed number of shells
#define DL_MAX 50

// macros to map the 10 spins of the 15j to the 4 booster symbols

#define MAP_SPINS_2(two_ja, two_jb, two_jc, two_jd, gf) \
    two_ja = two_j23; two_jb = two_j24; two_jc = two_j25; two_jd = two_j12; gf = 4;

#define MAP_SPINS_3(two_ja, two_jb, two_jc, two_jd, gf) \
    two_ja = two_j34; two_jb = two_j35; two_jc = two_j13; two_jd = two_j23; gf = 3;

#define MAP_SPINS_4(two_ja, two_jb, two_jc, two_jd, gf) \
    two_ja = two_j45; two_jb = two_j14; two_jc = two_j24; two_jd = two_j34; gf = 2;

#define MAP_SPINS_5(two_ja, two_jb, two_jc, two_jd, gf) \
    two_ja = two_j15; two_jb = two_j25; two_jc = two_j35; two_jd = two_j45; gf = 1;

sl2cfoam_tensor_boosters* sl2cfoam_boosters(int gf,
                                            dspin two_ja, dspin two_jb, dspin two_jc,  dspin two_jd, 
                                            int Dl, bool store) {

    MPI_FUNC_INIT();

    if (Dl > DL_MAX)
        error("too many shells requested, maximum is %d shells", DL_MAX);

    dspin two_Dl = (dspin)(2 * Dl);

    // get folder from config
    struct sl2cfoam_data_folders* df = get_data_folders();
    check_data_folders(df);

    // build filename
    char filename[256];
    sprintf(filename, boosters_fn, two_ja, two_jb, two_jc, two_jd, gf, IMMIRZI, Dl);

    const char* folder = df->vertex_imm_boost;

    // build path
    char path[strlen(folder) + 256];
    char path_found[strlen(folder) + 256];

    strcpy(path, folder);
    strcat(path, "/");
    strcat(path, filename);

    // result tensor
    tensor_ptr(boosters) b4t = NULL;

    // check if tensor already exists
    // if yes return it
    if (file_exist(path)) {
        TENSOR_LOAD(boosters, b4t, 6, path);
        goto tensor_return;
    }

    // check for tensors with same spin but different Dl
    bool found = false;
    dspin two_Dl_found;

    // go "up" from Dl+1 to DL_MAX then reset to 0 up to Dl (then stop)

    int dli = Dl + 1;
    do {

        // build filename
        sprintf(filename, boosters_fn, two_ja, two_jb, two_jc, two_jd, gf, IMMIRZI, dli);

        strcpy(path_found, folder);
        strcat(path_found, "/");
        strcat(path_found, filename);

        if (file_exist(path_found)) {
            found = true;
            two_Dl_found = (dspin)(2 * dli);
        }

        dli++;

        if (dli > DL_MAX) dli = 0;

    } while (found != true && dli != Dl);

    // must compute (at least partially)

    tensor_ptr(boosters) b4t_found;

    size_t ldim, idim, kdim_absmax;

    dspin two_i_min = max(abs(two_ja-two_jb), abs(two_jc-two_jd));
    dspin two_i_max = min(two_ja+two_jb, two_jc+two_jd);

    ldim = Dl + 1;
    idim = DIV2(two_i_max - two_i_min) + 1;

    // maximum dimension for intw space, k index
    dspin two_k_absmin, two_k_absmax;
    find_k_absolute_bounds(&kdim_absmax, &two_k_absmin, &two_k_absmax, 
                           two_ja, two_jb, two_jc, two_jd, two_Dl, gf);

    // do nothing for degenerate cases
    if (kdim_absmax == 0) goto tensor_return;

    // create tensor and set all entries to 0
    // the gauge-fixed index has dimension 1
    switch (gf)
    {
        case 1:
            TENSOR_CREATE(boosters, b4t, 6, idim, kdim_absmax, 1, ldim, ldim, ldim);
            break;
        
        case 2:
            TENSOR_CREATE(boosters, b4t, 6, idim, kdim_absmax, ldim, 1, ldim, ldim);
            break;

        case 3:
            TENSOR_CREATE(boosters, b4t, 6, idim, kdim_absmax, ldim, ldim, 1, ldim);
            break;

        case 4:
            TENSOR_CREATE(boosters, b4t, 6, idim, kdim_absmax, ldim, ldim, ldim, 1);
            break;
        
        default:
            error("gauge-fixed index must be 1 to 4");
    }

    // if a previous tensor with same spins has been found, load it
    if (found) {
        TENSOR_LOAD(boosters, b4t_found, 6, path_found);
    }

    // trick: collapse the for loop for the gauge fixed index
    dspin two_la_max, two_lb_max, two_lc_max, two_ld_max;
    switch (gf)
    {
        case 1:
            two_la_max = two_ja;
            two_lb_max = two_jb + two_Dl;
            two_lc_max = two_jc + two_Dl;
            two_ld_max = two_jd + two_Dl;
            
            break;
        
        case 2:
            two_la_max = two_ja + two_Dl;
            two_lb_max = two_jb;
            two_lc_max = two_jc + two_Dl;
            two_ld_max = two_jd + two_Dl;
            break;

        case 3:
            two_la_max = two_ja + two_Dl;
            two_lb_max = two_jb + two_Dl;
            two_lc_max = two_jc;
            two_ld_max = two_jd + two_Dl;
            break;

        case 4:
            two_la_max = two_ja + two_Dl;
            two_lb_max = two_jb + two_Dl;
            two_lc_max = two_jc + two_Dl;
            two_ld_max = two_jd;
            break;
        
        default:
            error("gauge-fixed index must be 1 to 4");
    }

    // this is needed in case of reusing found tensor
    dspin two_la_max_found, two_lb_max_found, two_lc_max_found, two_ld_max_found;
    if (found) {
        switch (gf)
        {
            case 1:
                two_la_max_found = two_ja;
                two_lb_max_found = two_jb + two_Dl_found;
                two_lc_max_found = two_jc + two_Dl_found;
                two_ld_max_found = two_jd + two_Dl_found;
                
                break;
            
            case 2:
                two_la_max_found = two_ja + two_Dl_found;
                two_lb_max_found = two_jb;
                two_lc_max_found = two_jc + two_Dl_found;
                two_ld_max_found = two_jd + two_Dl_found;
                break;

            case 3:
                two_la_max_found = two_ja + two_Dl_found;
                two_lb_max_found = two_jb + two_Dl_found;
                two_lc_max_found = two_jc;
                two_ld_max_found = two_jd + two_Dl_found;
                break;

            case 4:
                two_la_max_found = two_ja + two_Dl_found;
                two_lb_max_found = two_jb + two_Dl_found;
                two_lc_max_found = two_jc + two_Dl_found;
                two_ld_max_found = two_jd;
                break;
            
            default:
                error("gauge-fixed index must be 1 to 4");
        }
    }

    //////////////////////////////////////////////////////////////////////
    // compute the virtual spins ls to be summed over
    // - in MPI the ls are interleaved across the available nodes and then
    //   each node parallelize over his own ls
    //////////////////////////////////////////////////////////////////////

    // TODO: this can  be done more efficiently, but for a reasonable number
    //       of shells (~10 or little more) the loops are ok..
    long l_loops = CUBE(DIV2(two_Dl) + 1);
    dspin* ls_todo = (dspin*)calloc(l_loops, 4 * sizeof(dspin));
    size_t ls_todo_size = 0;

    long lcount = -1;
    for (dspin two_ld = two_jd; two_ld <= two_ld_max; two_ld += 2) {
    for (dspin two_lc = two_jc; two_lc <= two_lc_max; two_lc += 2) {
    for (dspin two_lb = two_jb; two_lb <= two_lb_max; two_lb += 2) {
    for (dspin two_la = two_ja; two_la <= two_la_max; two_la += 2) {

        lcount++;

        #ifdef USE_MPI

        // interleave across nodes
        if ((lcount % mpi_size) != mpi_rank) continue;

        #endif 

        ls_todo[ls_todo_size * 4 + 0] = two_la;
        ls_todo[ls_todo_size * 4 + 1] = two_lb;
        ls_todo[ls_todo_size * 4 + 2] = two_lc;
        ls_todo[ls_todo_size * 4 + 3] = two_ld;
        ls_todo_size++;

    } // la
    } // lb
    } // lc
    } // ld

    // check if to enable parallelization or not at this level
    // do NOT if the ls to compute (on this node if MPI)
    // are less than the number of cores
    bool go_parallel = true;
    if (ls_todo_size < omp_get_max_threads()) go_parallel = false;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(dynamic, 1) if(OMP_PARALLELIZE && go_parallel)
    #endif
    for (size_t lind = 0; lind < ls_todo_size; lind++) {

        dspin two_la = ls_todo[lind * 4 + 0];
        dspin two_lb = ls_todo[lind * 4 + 1];
        dspin two_lc = ls_todo[lind * 4 + 2];
        dspin two_ld = ls_todo[lind * 4 + 3];

        dspin two_k_min = max(abs(two_la-two_lb), abs(two_lc-two_ld));
        dspin two_k_max = min(two_la+two_lb, two_lc+two_ld);

        // check if there are no allowed intertwiner for this ls
        if (two_k_max < two_k_min) continue;

        dspin kdim = DIV2(two_k_max - two_k_min) + 1;

        // matrices with found values
        // NB: the matrices have same number of rows but different number of columns
        //     can copy whole matrix (thanks to column-major)
        sl2cfoam_dmatrix dst = b4t->d + TENSOR_INDEX(b4t, 6, 0, 0, DIV2(two_la-two_ja), DIV2(two_lb-two_jb), DIV2(two_lc-two_jc), DIV2(two_ld-two_jd));
        sl2cfoam_dmatrix src;

        if (found && two_Dl < two_Dl_found) {

            // here two_Dl_found > two_Dl so
            // just take the values from the tensor already computed
            src = b4t_found->d + TENSOR_INDEX(b4t_found, 6, 0, 0, DIV2(two_la-two_ja), DIV2(two_lb-two_jb), DIV2(two_lc-two_jc), DIV2(two_ld-two_jd));
            
            memcpy(dst, src, idim * min(kdim_absmax, b4t_found->dims[1]) * sizeof(double));
            continue;

        } 

        if (found && two_la <= two_la_max_found
                  && two_lb <= two_lb_max_found
                  && two_lc <= two_lc_max_found
                  && two_ld <= two_ld_max_found) {

            // here two_Dl_found < two_Dl but case covered by found tensor,
            // again just take the values from the tensor already computed
            src = b4t_found->d + TENSOR_INDEX(b4t_found, 6, 0, 0, DIV2(two_la-two_ja), DIV2(two_lb-two_jb), DIV2(two_lc-two_jc), DIV2(two_ld-two_jd));
            
            memcpy(dst, src, idim * min(kdim_absmax, b4t_found->dims[1]) * sizeof(double));
            continue;

        }

        // must compute

        sl2cfoam_dmatrix b4_ik = sl2cfoam_b4(two_ja, two_jb, two_jc, two_jd,
                                             two_la, two_lb, two_lc, two_ld);

        memcpy(dst, b4_ik, idim * kdim * sizeof(double));

        matrix_free(b4_ik);

    }

    #ifdef USE_MPI

    // reduce tensors over all nodes to master
    if (mpi_rank == MPI_MASTER) {
        MPI_Reduce(MPI_IN_PLACE, b4t->d, b4t->dim, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(b4t->d, b4t->d, b4t->dim, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    }

    #endif

    if (store) {
        MPI_MASTERONLY_DO TENSOR_STORE(b4t, path);
    }

    if (found) {
       TENSOR_FREE(b4t_found); 
    }

    #ifdef USE_MPI

    // synchronize nodes to ensure the stored tensor
    // is seen by all from now on
    MPI_Barrier(MPI_COMM_WORLD);

    #endif

tensor_return:

    free_data_folders(df);

    return b4t;
    
}

void sl2cfoam_boosters_tensors_vertex(dspin two_js[10], int Dl,
                                      tensor_ptr(boosters)* b2, tensor_ptr(boosters)* b3,
                                      tensor_ptr(boosters)* b4, tensor_ptr(boosters)* b5,
                                      dspin b_two_i_mins[4]) {

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    two_j12 = two_js[0];
    two_j13 = two_js[1];
    two_j14 = two_js[2];
    two_j15 = two_js[3];
    two_j23 = two_js[4];
    two_j24 = two_js[5];
    two_j25 = two_js[6];
    two_j34 = two_js[7];
    two_j35 = two_js[8];
    two_j45 = two_js[9];

    int gf;
    dspin two_ja, two_jb, two_jc, two_jd;

    // compute all boosters and store the tensors

    // booster 2
    MAP_SPINS_2(two_ja, two_jb, two_jc, two_jd, gf);

    *b2 = sl2cfoam_boosters(gf, two_ja, two_jb, two_jc, two_jd, Dl, true);
    b_two_i_mins[0] = max(abs(two_ja-two_jb), abs(two_jc-two_jd));

    // booster 3
    MAP_SPINS_3(two_ja, two_jb, two_jc, two_jd, gf);

    *b3 = sl2cfoam_boosters(gf, two_ja, two_jb, two_jc, two_jd, Dl, true);
    b_two_i_mins[1] = max(abs(two_ja-two_jb), abs(two_jc-two_jd));

    // booster 4
    MAP_SPINS_4(two_ja, two_jb, two_jc, two_jd, gf);

    *b4 = sl2cfoam_boosters(gf, two_ja, two_jb, two_jc, two_jd, Dl, true);
    b_two_i_mins[2] = max(abs(two_ja-two_jb), abs(two_jc-two_jd));

    // booster 5
    MAP_SPINS_5(two_ja, two_jb, two_jc, two_jd, gf);

    *b5 = sl2cfoam_boosters(gf, two_ja, two_jb, two_jc, two_jd, Dl, true);
    b_two_i_mins[3] = max(abs(two_ja-two_jb), abs(two_jc-two_jd));

}

sl2cfoam_tensor_boosters* sl2cfoam_boosters_load(int gf,
                                                 dspin two_ja, dspin two_jb, dspin two_jc,  dspin two_jd, 
                                                 int Dl) {

    // get folder from config
    struct sl2cfoam_data_folders* df = get_data_folders();
    check_data_folders(df);

    // build filename
    char filename[256];
    sprintf(filename, boosters_fn, two_ja, two_jb, two_jc, two_jd, gf, IMMIRZI, Dl);

    const char* folder = df->vertex_imm_boost;

    // build path
    char path[strlen(folder) + 256];

    strcpy(path, folder);
    strcat(path, "/");
    strcat(path, filename);

    tensor_ptr(boosters) t;
    TENSOR_LOAD(boosters, t, 6, path);

    free_data_folders(df);

    return t;

}

void sl2cfoam_boosters_free(sl2cfoam_tensor_boosters* t) {
    TENSOR_FREE(t);
}
