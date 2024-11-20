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
#include <omp.h>

#include "common.h"
#include "utils.h"
#include "mpi_utils.h"
#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"
#include "boosters.h"
#include "jsymbols.h"
#include "verb.h"
#include "error.h"
#include "blas_wrapper.h"
#include "timing.h"


// filenames for storing the tensors
const char* fn_template = "vertex__%d-%d-%d-%d-%d-%d-%d-%d-%d-%d__i1_%d-%d__i2_%d-%d__i3_%d-%d__i4_%d-%d__i5_%d-%d__imm-%.3f__dl-%d.sl2t";
const char* fn_template_fullrange = "vertex__%d-%d-%d-%d-%d-%d-%d-%d-%d-%d__fullrange__imm-%.3f__dl-%d.sl2t";

// build the file path of a vertex tensor (checking for fullrange case)
static void build_vertex_path(char* path,  sl2cfoam_dspin two_js[10],
                                          sl2cfoam_dspin two_i1_min, sl2cfoam_dspin two_i1_max, 
                                          sl2cfoam_dspin two_i2_min, sl2cfoam_dspin two_i2_max, 
                                          sl2cfoam_dspin two_i3_min, sl2cfoam_dspin two_i3_max, 
                                          sl2cfoam_dspin two_i4_min, sl2cfoam_dspin two_i4_max, 
                                          sl2cfoam_dspin two_i5_min, sl2cfoam_dspin two_i5_max, 
                                          int Dl) {

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    // check if given intertwiners are allowed
    dspin two_i1_min_allowed, two_i1_max_allowed;
    dspin two_i2_min_allowed, two_i2_max_allowed;
    dspin two_i3_min_allowed, two_i3_max_allowed;
    dspin two_i4_min_allowed, two_i4_max_allowed;
    dspin two_i5_min_allowed, two_i5_max_allowed;

    ASSIGN_INTW_RANGES(allowed);
    
    char filename[512];

    if (two_i1_min == two_i1_min_allowed && two_i1_max == two_i1_max_allowed &&
        two_i2_min == two_i2_min_allowed && two_i2_max == two_i2_max_allowed &&
        two_i3_min == two_i3_min_allowed && two_i3_max == two_i3_max_allowed &&
        two_i4_min == two_i4_min_allowed && two_i4_max == two_i4_max_allowed &&
        two_i5_min == two_i5_min_allowed && two_i5_max == two_i5_max_allowed) {

        sprintf(filename, fn_template_fullrange, two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                                                 two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                                                 IMMIRZI, Dl);

    } else {

        sprintf(filename, fn_template, two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                                       two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                                       two_i1_min, two_i1_max, 
                                       two_i2_min, two_i2_max, 
                                       two_i3_min, two_i3_max, 
                                       two_i4_min, two_i4_max, 
                                       two_i5_min, two_i5_max,
                                       IMMIRZI, Dl);

    }

    strcpy(path, DIR_AMPLS);
    strcat(path, "/");
    strcat(path, filename);

}

// build the file path of a vertex tensor (assuming fullrange)
static void build_vertex_path_fullrange(char* path, sl2cfoam_dspin two_js[10], int Dl) {

    char filename[512];

    sprintf(filename, fn_template_fullrange, two_js[0], two_js[1], two_js[2], two_js[3], two_js[4],
                                             two_js[5], two_js[6], two_js[7], two_js[8], two_js[9],
                                             IMMIRZI, Dl);

    strcpy(path, DIR_AMPLS);
    strcat(path, "/");
    strcat(path, filename);

}

// computes the approximate memory needed by the matrices used in a single thread
// returns in MBytes
static size_t max_mem_needed_per_thread(dspin two_js[10], 
                                        dspin two_i1_min, dspin two_i1_max, 
                                        dspin two_i2_min, dspin two_i2_max, 
                                        dspin two_i3_min, dspin two_i3_max, 
                                        dspin two_i4_min, dspin two_i4_max, 
                                        dspin two_i5_min, dspin two_i5_max,
                                        int Dl) {

    dspin two_Dl = (dspin)(2 * Dl);

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    size_t i1_size, i2_size, i3_size, i4_size, i5_size;
    i1_size = DIV2(two_i1_max-two_i1_min) + 1;
    i2_size = DIV2(two_i2_max-two_i2_min) + 1;
    i3_size = DIV2(two_i3_max-two_i3_min) + 1;
    i4_size = DIV2(two_i4_max-two_i4_min) + 1;
    i5_size = DIV2(two_i5_max-two_i5_min) + 1;

    size_t k2_size, k3_size, k4_size, k5_size;
    dspin two_k2_absmin, two_k2_absmax, two_k3_absmin, two_k3_absmax,
          two_k4_absmin, two_k4_absmax, two_k5_absmin, two_k5_absmax;

    find_k_absolute_bounds(&k2_size, &two_k2_absmin, &two_k2_absmax, two_j23, two_j24, two_j25, two_j12, two_Dl, 4);
    find_k_absolute_bounds(&k3_size, &two_k3_absmin, &two_k3_absmax, two_j34, two_j35, two_j13, two_j23, two_Dl, 3);
    find_k_absolute_bounds(&k4_size, &two_k4_absmin, &two_k4_absmax, two_j45, two_j14, two_j24, two_j34, two_Dl, 2);
    find_k_absolute_bounds(&k5_size, &two_k5_absmin, &two_k5_absmax, two_j15, two_j25, two_j35, two_j45, two_Dl, 1);

    size_t tot_bytes = 0;

    tot_bytes += sizeof(double) * k2_size * k3_size * k4_size * k5_size * i1_size; // internal 15j tensor
    tot_bytes += sizeof(double) * k3_size * k4_size * k5_size * i1_size * i2_size; // temp contraction matrix (pampl_tmp)
    tot_bytes += sizeof(double) * i1_size * i2_size * i3_size * i4_size * i5_size; // final tensor (pampl)

    size_t tot_MB = tot_bytes / (1 << 20); // integer truncated division
    return tot_MB;

}

double sl2cfoam_vertex_amplitude(dspin two_js[10], dspin two_is[5], int Dl) {

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    dspin two_i1, two_i2, two_i3, two_i4, two_i5;
    two_i1 = two_is[0];
    two_i2 = two_is[1];
    two_i3 = two_is[2];
    two_i4 = two_is[3];
    two_i5 = two_is[4];

    // check if given intertwiners are allowed
    dspin two_i1_min_allowed, two_i1_max_allowed;
    dspin two_i2_min_allowed, two_i2_max_allowed;
    dspin two_i3_min_allowed, two_i3_max_allowed;
    dspin two_i4_min_allowed, two_i4_max_allowed;
    dspin two_i5_min_allowed, two_i5_max_allowed;

    ASSIGN_INTW_RANGES(allowed);

    // checks
    bool check_failed = false;

    CHECK_NULL_INTW(i1, check_failed);
    CHECK_NULL_INTW(i2, check_failed);
    CHECK_NULL_INTW(i3, check_failed);
    CHECK_NULL_INTW(i4, check_failed);
    CHECK_NULL_INTW(i5, check_failed);

    if (check_failed) return 0.0;

    CHECK_ALLOWED_INTW(i1, check_failed);
    CHECK_ALLOWED_INTW(i2, check_failed);
    CHECK_ALLOWED_INTW(i3, check_failed);
    CHECK_ALLOWED_INTW(i4, check_failed);
    CHECK_ALLOWED_INTW(i5, check_failed);

    if (check_failed) return 0.0;

    sl2cfoam_tensor_vertex* t = sl2cfoam_vertex_range(two_js,
                                                      two_is[0], two_is[0],
                                                      two_is[1], two_is[1], 
                                                      two_is[2], two_is[2],  
                                                      two_is[3], two_is[3],  
                                                      two_is[4], two_is[4],  
                                                      Dl, TENSOR_RESULT_RETURN);

    // this code covers also MPI version
    // junk (INT32_MAX) is returned by nodes different from master
    double ampl = (double)INT32_MAX;

    if (t != NULL) {

        ampl = TENSOR_GET(t, 5, 0, 0, 0, 0, 0);
        TENSOR_FREE(t);

    }

    return ampl;

}

sl2cfoam_tensor_vertex* sl2cfoam_vertex_fullrange(dspin two_js[10], int Dl, sl2cfoam_tensor_result tresult) {

    MPI_FUNC_INIT();

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    // compute the intertwiners ranges
    dspin two_i1_min_allowed, two_i1_max_allowed;
    dspin two_i2_min_allowed, two_i2_max_allowed;
    dspin two_i3_min_allowed, two_i3_max_allowed;
    dspin two_i4_min_allowed, two_i4_max_allowed;
    dspin two_i5_min_allowed, two_i5_max_allowed;

    ASSIGN_INTW_RANGES(allowed);

    // full tensor
    tensor_ptr(vertex) t_full = NULL;

    // look for previously computed tensor

#ifndef NO_IO
    
    char path[strlen(DIR_AMPLS) + 512];
    build_vertex_path_fullrange(path, two_js, Dl);

    if (file_exist(path)) {
        
        MPI_MASTERONLY_START

        verb(SL2CFOAM_VERBOSE_HIGH, "[ previously computed full tensor found ]\n");

        if (bitmask_has(tresult, TENSOR_RESULT_RETURN)) {
                
            TENSOR_LOAD(vertex, t_full, 5, path);

            if (t_full == NULL) {
                verb(SL2CFOAM_VERBOSE_LOW, "Error loading previously computed full tensor\n");
            }

            return t_full;

        }
            
        MPI_MASTERONLY_END

        return NULL;

    }

#endif

    // no stored tensor found, computing

    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_LOW, "Computing full tensor...\n");
    MPI_MASTERONLY_DO TIC();  

    // check if the max memory is capped
    // if yes the computation is batched
    bool batch_computation = false;

    size_t needed_MB;
    size_t max_MB = CONFIG.max_MB_mem_per_thread;
    if (max_MB > 0) {

        needed_MB = max_mem_needed_per_thread(two_js,
                                              two_i1_min_allowed, two_i1_max_allowed, 
                                              two_i2_min_allowed, two_i2_max_allowed, 
                                              two_i3_min_allowed, two_i3_max_allowed, 
                                              two_i4_min_allowed, two_i4_max_allowed, 
                                              two_i5_min_allowed, two_i5_max_allowed, 
                                              Dl);

        // batch if memory is not enough (required memory is more than 90% available)
        if ((double)needed_MB > (0.9 * max_MB)) {

            batch_computation = true;

        }

    }

    if (batch_computation) {

        bool build_full = bitmask_has(tresult, TENSOR_RESULT_RETURN) || bitmask_has(tresult, TENSOR_RESULT_STORE);

        size_t i1_size, i2_size, i3_size, i4_size, i5_size;
        i1_size = DIV2(two_i1_max_allowed-two_i1_min_allowed) + 1;
        i2_size = DIV2(two_i2_max_allowed-two_i2_min_allowed) + 1;
        i3_size = DIV2(two_i3_max_allowed-two_i3_min_allowed) + 1;
        i4_size = DIV2(two_i4_max_allowed-two_i4_min_allowed) + 1;
        i5_size = DIV2(two_i5_max_allowed-two_i5_min_allowed) + 1;

        // restrict the range of i1 to fit memory
        // the memory needed is a multiple of i1_size
        size_t i1_size_per_thread;
        i1_size_per_thread = (size_t)round((0.9 * max_MB) / (needed_MB / (double)i1_size));

        if (i1_size_per_thread < 1) error("cannot fit memory limit per thread");

        // compute the needed intervals
        int num_batches = i1_size / i1_size_per_thread + (i1_size % i1_size_per_thread == 0 ? 0 : 1);
        MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_LOW, "[ batching full tensor into %d small tensors ]\n", num_batches);

        // flags for batch tensor computation
        sl2cfoam_tensor_result bres = 0;
        if (build_full) bitmask_set(bres, TENSOR_RESULT_RETURN);
        if (bitmask_has(tresult, TENSOR_RESULT_STORE_BATCHES)) bitmask_set(bres, TENSOR_RESULT_STORE);

        dspin two_i1_min_batch, two_i1_max_batch;
        two_i1_min_batch = two_i1_min_allowed;

MPI_MASTERONLY_START

        if (build_full) TENSOR_CREATE(vertex, t_full, 5, i5_size, i4_size, i3_size, i2_size, i1_size);

MPI_MASTERONLY_END

        tensor_ptr(vertex) tbatch;
        size_t skip = 0;
        while (true) {

            two_i1_max_batch = min(two_i1_min_batch + 2*(i1_size_per_thread - 1), two_i1_max_allowed);
            
            MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "[ batching: computing two_i1 from %d to %d... ]\n", two_i1_min_batch, two_i1_max_batch);
            tbatch = sl2cfoam_vertex_range(two_js,
                                           two_i1_min_batch, two_i1_max_batch, 
                                           two_i2_min_allowed, two_i2_max_allowed, 
                                           two_i3_min_allowed, two_i3_max_allowed, 
                                           two_i4_min_allowed, two_i4_max_allowed, 
                                           two_i5_min_allowed, two_i5_max_allowed, 
                                           Dl, bres);

MPI_MASTERONLY_START

            // accumulate if full tensor must be built
            if (build_full) {

                memcpy(t_full->d + skip, tbatch->d, tbatch->dim * sizeof(double));
                skip += tbatch->dim;
                TENSOR_FREE(tbatch);

            }

MPI_MASTERONLY_END            

            if (two_i1_max_batch == two_i1_max_allowed) break;

            two_i1_min_batch = two_i1_max_batch + 2;

        }

    } else {

        // do not batch
        t_full = sl2cfoam_vertex_range(two_js,
                                       two_i1_min_allowed, two_i1_max_allowed, 
                                       two_i2_min_allowed, two_i2_max_allowed, 
                                       two_i3_min_allowed, two_i3_max_allowed, 
                                       two_i4_min_allowed, two_i4_max_allowed, 
                                       two_i5_min_allowed, two_i5_max_allowed, 
                                       Dl, TENSOR_RESULT_RETURN);

    }

    MPI_MASTERONLY_DO TOC();
    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_LOW, "... done in %.3f seconds.\n", ELAPSED());

MPI_MASTERONLY_START

    // add tag
    if (t_full != NULL) {

        sl2cfoam_vertex_tag tag;

        tag.Immirzi = IMMIRZI;
        memcpy(tag.two_js, two_js, 10 * sizeof(dspin));
        tag.two_i1_range[0] = two_i1_min_allowed;
        tag.two_i1_range[1] = two_i1_max_allowed;
        tag.two_i2_range[0] = two_i2_min_allowed;
        tag.two_i2_range[1] = two_i2_max_allowed;
        tag.two_i3_range[0] = two_i3_min_allowed;
        tag.two_i3_range[1] = two_i3_max_allowed;
        tag.two_i4_range[0] = two_i4_min_allowed;
        tag.two_i4_range[1] = two_i4_max_allowed;
        tag.two_i5_range[0] = two_i5_min_allowed;
        tag.two_i5_range[1] = two_i5_max_allowed;
        tag.Dl = Dl;

        TENSOR_TAG(t_full, &tag, sizeof(tag));

    }

#ifndef NO_IO

    // store if configured
    if (t_full != NULL && bitmask_has(tresult, TENSOR_RESULT_STORE)) {
        TENSOR_STORE(t_full, path);
    }

#endif

    // return the result if configured (and master)
    if (bitmask_has(tresult, TENSOR_RESULT_RETURN))
        return t_full;

MPI_MASTERONLY_END
        
    if (t_full != NULL) TENSOR_FREE(t_full);
    return NULL;

}

// tensor for internal 15j in the intertwiner indices
// indices: (i1 k5 k4 k3 k2)
TENSOR_INIT(15j, 5)

sl2cfoam_tensor_vertex* sl2cfoam_vertex_range(dspin two_js[10],
                                              dspin two_i1_min, dspin two_i1_max, 
                                              dspin two_i2_min, dspin two_i2_max, 
                                              dspin two_i3_min, dspin two_i3_max, 
                                              dspin two_i4_min, dspin two_i4_max, 
                                              dspin two_i5_min, dspin two_i5_max, 
                                              int Dl, sl2cfoam_tensor_result tresult) {

    MPI_FUNC_INIT();

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

    dspin two_Dl = (dspin)(2 * Dl);

    // check if given range is allowed
    dspin two_i1_min_allowed, two_i1_max_allowed;
    dspin two_i2_min_allowed, two_i2_max_allowed;
    dspin two_i3_min_allowed, two_i3_max_allowed;
    dspin two_i4_min_allowed, two_i4_max_allowed;
    dspin two_i5_min_allowed, two_i5_max_allowed;

    ASSIGN_INTW_RANGES(allowed);

    // checks
    bool check_failed = false;

    CHECK_NULL_INTW(i1, check_failed);
    CHECK_NULL_INTW(i2, check_failed);
    CHECK_NULL_INTW(i3, check_failed);
    CHECK_NULL_INTW(i4, check_failed);
    CHECK_NULL_INTW(i5, check_failed);

    if (check_failed) return NULL;

    CHECK_ALLOWED_INTW_RANGE(i1, check_failed);
    CHECK_ALLOWED_INTW_RANGE(i2, check_failed);
    CHECK_ALLOWED_INTW_RANGE(i3, check_failed);
    CHECK_ALLOWED_INTW_RANGE(i4, check_failed);
    CHECK_ALLOWED_INTW_RANGE(i5, check_failed);

    if (check_failed) return NULL;

    // spins look good

    // result tensor
    tensor_ptr(vertex) t_range = NULL;

    // found tensor (if found)
    tensor_ptr(vertex) t_found = NULL;

    ///////////////////////////////////////////////////////////////
    // search for stored tensors with lower (or equal) Dl
    ///////////////////////////////////////////////////////////////

    bool found = false;
    dspin two_Dl_found = 0;

MPI_MASTERONLY_START
#ifndef NO_IO

    // check for computed tensors with lower or equal Dl
    
    char path_found[strlen(DIR_AMPLS) + 512];

    for (dspin two_dli = two_Dl; two_dli >= 0; two_dli -= 2) {

        int Dl_found = DIV2(two_dli);

        build_vertex_path(path_found, two_js, 
                                      two_i1_min, two_i1_max, 
                                      two_i2_min, two_i2_max, 
                                      two_i3_min, two_i3_max, 
                                      two_i4_min, two_i4_max, 
                                      two_i5_min, two_i5_max,
                                      Dl_found);

        // check if requested tensor has already been computed
        if (file_exist(path_found)) {

            verb(SL2CFOAM_VERBOSE_HIGH, "[ previously computed partial tensor found ]\n");

            TENSOR_LOAD(vertex, t_found, 5, path_found);

            if (t_found == NULL) {
                verb(SL2CFOAM_VERBOSE_LOW, "Error loading previously computed partial tensor\n");
            }

            if (t_found != NULL) {
                
                found = true;
                two_Dl_found = two_dli;
                break;

            }

        }

    }

#endif
MPI_MASTERONLY_END

    #ifdef USE_MPI

    // if MPI broadcast the found values to all nodes
    MPI_Bcast(&found, 1, MPI_C_BOOL, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&two_Dl_found, 1, MPI_INT, MPI_MASTER, MPI_COMM_WORLD);

    #endif

    if (found) {

        if (two_Dl_found == two_Dl) {
            t_range = t_found;
            goto tensor_range_return;
        }

        MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, 
            "[ found vertex tensor with lower Dl = %d, accumulating... ]\n", DIV2(two_Dl_found));
        
    } 

    // compute dimensions
    size_t i1_size, i2_size, i3_size, i4_size, i5_size;
    i1_size = DIV2(two_i1_max-two_i1_min) + 1;
    i2_size = DIV2(two_i2_max-two_i2_min) + 1;
    i3_size = DIV2(two_i3_max-two_i3_min) + 1;
    i4_size = DIV2(two_i4_max-two_i4_min) + 1;
    i5_size = DIV2(two_i5_max-two_i5_min) + 1;

    // create result tensor
    TENSOR_CREATE(vertex, t_range, 5, i5_size, i4_size, i3_size, i2_size, i1_size);

    // in MPI mode each node owns his own copy of the tensor
    // then finally reduce over MPI_MASTER
    size_t tot_size = i1_size * i2_size * i3_size * i4_size * i5_size;

    ///////////////////////////////////////////////////////////////
    // compute 6js tensors
    ///////////////////////////////////////////////////////////////

    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "Computing 6js tensors...\n");
    MPI_MASTERONLY_DO TIC();

    tensor_ptr(6j_4) sjA;
    tensor_ptr(6j_5) sjB;
    tensor_ptr(6j_6) sjC;
    tensor_ptr(6j_5) sjD;
    tensor_ptr(6j_4) sjE;

    dspin bounds[10];

    sl2cfoam_6j_tensors_vertex(two_js, two_i1_min, two_i1_max, two_Dl, &sjA, &sjB, &sjC, &sjD, &sjE, bounds);

    MPI_MASTERONLY_DO TOC();
    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "... done in %.3f seconds.\n", ELAPSED());

    // compute k range for internal 15j assembly
    dspin two_k2_absmin, two_k3_absmin, two_k4_absmin, two_k5_absmin;

    two_k2_absmin = bounds[0];
    two_k3_absmin = bounds[2];
    two_k4_absmin = bounds[4];
    two_k5_absmin = bounds[6];

    dspin two_i1_absmin = two_i1_min;

    dspin two_x_absmin = bounds[8];
    dspin two_x_absmax = bounds[9];

    // sometimes triangular inequalities leave nothing to compute
    // (e.g. in batching)
    if (two_x_absmax < two_x_absmin) goto tensor_range_return;

    size_t x_size = DIV2(two_x_absmax - two_x_absmin) + 1;

    sl2cfoam_dvector xdims = dvector_alloc(x_size);
    for (dspin two_x = two_x_absmin; two_x <= two_x_absmax; two_x += 2) {
        xdims[DIV2(two_x-two_x_absmin)] = (double)DIM(two_x);
    }

    ///////////////////////////////////////////////////////////////
    // precompute boosters tensors
    ///////////////////////////////////////////////////////////////

    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "Computing booster tensors...\n");
    MPI_MASTERONLY_DO TIC();

    tensor_ptr(boosters) booster_2;
    tensor_ptr(boosters) booster_3;
    tensor_ptr(boosters) booster_4;
    tensor_ptr(boosters) booster_5;

    dspin boosters_two_i_min[4];

    sl2cfoam_boosters_tensors_vertex(two_js, Dl, &booster_2, &booster_3, &booster_4, &booster_5, boosters_two_i_min);

    dspin b2_two_i_min, b3_two_i_min, b4_two_i_min, b5_two_i_min;

    b2_two_i_min = boosters_two_i_min[0];
    b3_two_i_min = boosters_two_i_min[1];
    b4_two_i_min = boosters_two_i_min[2];
    b5_two_i_min = boosters_two_i_min[3];

    MPI_MASTERONLY_DO TOC();
    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "... done in %.3f seconds.\n", ELAPSED());    

    //////////////////////////////////////////////////////////////////////
    // compute the virtual spins ls to be summed over
    // - consider found tensors with lower Dl
    // - in MPI the ls are interleaved across the available nodes and then
    //   each node parallelize over his own ls
    //////////////////////////////////////////////////////////////////////

    // TODO: this can  be done more efficiently, but for a reasonable number
    //       of shells (~10 or little more) the loops are ok...
    long l_loops = SQ( (long)CUBE(Dl+1) );
    dspin* ls_todo = (dspin*)calloc(l_loops, 6 * sizeof(dspin));
    size_t ls_todo_size = 0;

    long lcount = -1;
    for (dspin two_l45 = two_j45; two_l45 <= two_j45 + two_Dl; two_l45 += 2) {
    for (dspin two_l35 = two_j35; two_l35 <= two_j35 + two_Dl; two_l35 += 2) {
    for (dspin two_l34 = two_j34; two_l34 <= two_j34 + two_Dl; two_l34 += 2) {
    for (dspin two_l25 = two_j25; two_l25 <= two_j25 + two_Dl; two_l25 += 2) {
    for (dspin two_l24 = two_j24; two_l24 <= two_j24 + two_Dl; two_l24 += 2) {
    for (dspin two_l23 = two_j23; two_l23 <= two_j23 + two_Dl; two_l23 += 2) {

        lcount++;

        // check if all ls have been computed in found tensor
        // if yes skip
        if (found && two_l45 <= two_j45 + two_Dl_found
                  && two_l35 <= two_j35 + two_Dl_found
                  && two_l34 <= two_j34 + two_Dl_found
                  && two_l25 <= two_j25 + two_Dl_found
                  && two_l24 <= two_j24 + two_Dl_found
                  && two_l23 <= two_j23 + two_Dl_found) continue;

        #ifdef USE_MPI

        // interleave across nodes
        if ((lcount % mpi_size) != mpi_rank) continue;

        #endif 

        ls_todo[ls_todo_size * 6 + 0] = two_l23;
        ls_todo[ls_todo_size * 6 + 1] = two_l24;
        ls_todo[ls_todo_size * 6 + 2] = two_l25;
        ls_todo[ls_todo_size * 6 + 3] = two_l34;
        ls_todo[ls_todo_size * 6 + 4] = two_l35;
        ls_todo[ls_todo_size * 6 + 5] = two_l45;
        ls_todo_size++;

    } // l23
    } // l24
    } // l25
    } // l34
    } // l35
    } // l45

    ///////////////////////////////////////////////////////////////
    // assemble boosters and 6js (via 15j)
    // start of parallel OMP code
    ///////////////////////////////////////////////////////////////

    // if there are enough enough shells to compute then
    // parallelize over the ls
    // otherwise parallelize internal 15j and DGEMMs
    // TODO: MKL is set to SEQUENTIAL in any case because of strange
    //       bug with the threaded version
    bool parallel_ls = true;
    int nthreads = omp_get_max_threads();
    if (ls_todo_size < nthreads) parallel_ls = false;

    MPI_MASTERONLY_DO verb(SL2CFOAM_VERBOSE_HIGH, "Assembling amplitudes...\n");
    MPI_MASTERONLY_DO TIC();

    #ifdef USE_OMP
    #pragma omp parallel if(OMP_PARALLELIZE && parallel_ls)
    {
    #endif

    size_t pampl_cols = i1_size * i2_size * i3_size * i4_size;
    sl2cfoam_dvector pampl = dvector_alloc(tot_size);

    // start (parallel) loop over the needed ls

    #ifdef USE_OMP
    #pragma omp for schedule(dynamic, 1)
    #endif
    for (size_t lind = 0; lind < ls_todo_size; lind++) {

        dspin two_l23 = ls_todo[lind * 6 + 0];
        dspin two_l24 = ls_todo[lind * 6 + 1];
        dspin two_l25 = ls_todo[lind * 6 + 2];
        dspin two_l34 = ls_todo[lind * 6 + 3];
        dspin two_l35 = ls_todo[lind * 6 + 4];
        dspin two_l45 = ls_todo[lind * 6 + 5];

        int il23, il24, il25, il34, il35, il45;

        il23 = DIV2(two_l23 - two_j23);
        il24 = DIV2(two_l24 - two_j24);
        il25 = DIV2(two_l25 - two_j25);
        il34 = DIV2(two_l34 - two_j34);
        il35 = DIV2(two_l35 - two_j35);
        il45 = DIV2(two_l45 - two_j45);
        
        dspin two_k2_min, two_k2_max, two_k3_min, two_k3_max,
              two_k4_min, two_k4_max, two_k5_min, two_k5_max;

        two_k2_min = max(abs(two_l23-two_l24), abs(two_l25-two_j12));
        two_k2_max = min(two_l23+two_l24, two_l25+two_j12);
        if (two_k2_max < two_k2_min) continue;

        two_k3_min = max(abs(two_l34-two_l35), abs(two_j13-two_l23));
        two_k3_max = min(two_l34+two_l35, two_j13+two_l23);
        if (two_k3_max < two_k3_min) continue;

        two_k4_min = max(abs(two_l45-two_j14), abs(two_l24-two_l34));
        two_k4_max = min(two_l45+two_j14, two_l24+two_l34);
        if (two_k4_max < two_k4_min) continue;

        two_k5_min = max(abs(two_j15-two_l25), abs(two_l35-two_l45));
        two_k5_max = min(two_j15+two_l25, two_l35+two_l45);
        if (two_k5_max < two_k5_min) continue;

        size_t k2_size_l, k3_size_l, k4_size_l, k5_size_l;
        k2_size_l = DIV2(two_k2_max - two_k2_min) + 1;
        k3_size_l = DIV2(two_k3_max - two_k3_min) + 1;
        k4_size_l = DIV2(two_k4_max - two_k4_min) + 1;
        k5_size_l = DIV2(two_k5_max - two_k5_min) + 1;

        dspin jl_sum = two_j12 + two_j13 + two_j14 + two_j15 + two_l23 +  
                       two_l24 + two_l25 + two_l34 + two_l35 + two_l45;

        // compute stricter bounds for x
        dspin two_x_absmin_l, two_x_absmax_l;

        dspin two_zero = is_integer(two_i1_min + two_j25) ? 0 : 1;

        two_x_absmin_l = max3(two_zero,       two_i1_min-two_l25, two_l25-two_i1_max);
        two_x_absmin_l = max3(two_x_absmin_l, two_k5_min-two_j14, two_j14-two_k5_max);
        two_x_absmin_l = max3(two_x_absmin_l, two_k4_min-two_l35, two_l35-two_k4_max);
        two_x_absmin_l = max3(two_x_absmin_l, two_k3_min-two_l24, two_l24-two_k3_max);
        two_x_absmin_l = max3(two_x_absmin_l, two_k2_min-two_j13, two_j13-two_k2_max);

        two_x_absmax_l = min(two_i1_max + two_l25, two_k5_max + two_j14);
        two_x_absmax_l = min(two_x_absmax_l, two_k4_max + two_l35);
        two_x_absmax_l = min(two_x_absmax_l, two_k3_max + two_l24);
        two_x_absmax_l = min(two_x_absmax_l, two_k2_max + two_j13);

        if (two_x_absmax_l < two_x_absmin_l) continue;

        int x_size_l = DIV2(two_x_absmax_l-two_x_absmin_l) + 1;

        // first index in larger x arrays
        int ix0 = DIV2(two_x_absmin_l - two_x_absmin);

        // build 15j

        // internal 15j symbol tensor
        tensor_ptr(15j) w15j_tensor;
        TENSOR_CREATE(15j, w15j_tensor, 5, i1_size, k5_size_l, k4_size_l, k3_size_l, k2_size_l);

        #ifdef USE_OMP
        #pragma omp parallel for collapse(2) if(OMP_PARALLELIZE && !parallel_ls)
        #endif
        for (dspin two_k2 = two_k2_min; two_k2 <= two_k2_max; two_k2 += 2) {
        for (dspin two_k3 = two_k3_min; two_k3 <= two_k3_max; two_k3 += 2) {

            int ik2 = DIV2(two_k2-two_k2_absmin);
            int ik3 = DIV2(two_k3-two_k3_absmin);
            const sl2cfoam_dvector w6j4 = sjD->d + TENSOR_INDEX(sjD, 5, 0, ik3, ik2, il23, il24);

        for (dspin two_k4 = two_k4_min; two_k4 <= two_k4_max; two_k4 += 2) {

            int ik4 = DIV2(two_k4-two_k4_absmin);
            const sl2cfoam_dvector w6j3 = sjC->d + TENSOR_INDEX(sjC, 6, 0, ik4, ik3, il24, il34, il35);

        for (dspin two_k5 = two_k5_min; two_k5 <= two_k5_max; two_k5 += 2) {

            int ik5 = DIV2(two_k5-two_k5_absmin);
            const sl2cfoam_dvector w6j2 = sjB->d + TENSOR_INDEX(sjB, 5, 0, ik5, ik4, il35, il45);

        for (dspin two_i1 = two_i1_min; two_i1 <= two_i1_max; two_i1 += 2) {

            int ii1 = DIV2(two_i1-two_i1_absmin);
            const sl2cfoam_dvector w6j1 = sjA->d + TENSOR_INDEX(sjA, 4, 0, ii1, ik5, il25);
            const sl2cfoam_dvector w6j5 = sjE->d + TENSOR_INDEX(sjE, 4, 0, ii1, ik2, il25);

            double w15j = 0.0;
            for (int i = ix0; i < ix0 + x_size_l; i++) {
                w15j += xdims[i] * w6j1[i] * w6j2[i] * w6j3[i] * w6j4[i] * w6j5[i];
            }

            if (w15j != 0) {
        
                // add sign
                w15j *= real_negpow(jl_sum + two_i1 + two_k2 + two_k3 + two_k4 + two_k5);

                // add dimensions
                w15j *= sqrt(DIM(two_i1) * DIM(two_k2) * DIM(two_k3) * DIM(two_k4) * DIM(two_k5));

                TENSOR_SET(w15j, w15j_tensor, 5, DIV2(two_i1-two_i1_min), DIV2(two_k5-two_k5_min), 
                                                 DIV2(two_k4-two_k4_min), DIV2(two_k3-two_k3_min), 
                                                 DIV2(two_k2-two_k2_min));

            }

        } // k2
        } // k3
        } // k4
        } // k5
        } // i1

        // get booster submatrices
        const sl2cfoam_dmatrix booster_2_subm = booster_2->d + TENSOR_INDEX(booster_2, 6, DIV2(two_i2_min-b2_two_i_min), 0, il23, il24, il25,    0);
        const sl2cfoam_dmatrix booster_3_subm = booster_3->d + TENSOR_INDEX(booster_3, 6, DIV2(two_i3_min-b3_two_i_min), 0, il34, il35,    0, il23);
        const sl2cfoam_dmatrix booster_4_subm = booster_4->d + TENSOR_INDEX(booster_4, 6, DIV2(two_i4_min-b4_two_i_min), 0, il45,    0, il24, il34);
        const sl2cfoam_dmatrix booster_5_subm = booster_5->d + TENSOR_INDEX(booster_5, 6, DIV2(two_i5_min-b5_two_i_min), 0,    0, il25, il35, il45);

        //////////////////////////////////////////////////////////////////
        // multiply using BLAS
        //////////////////////////////////////////////////////////////////

        // for contracting k indices
        size_t pampl_r1_cols = k3_size_l * k4_size_l * k5_size_l * i1_size;
        size_t pampl_r2_cols = k4_size_l * k5_size_l * i1_size * i2_size;
        size_t pampl_r3_cols = k5_size_l * i1_size * i2_size * i3_size;

        sl2cfoam_dmatrix pampl_tmp = dmatrix_alloc(i2_size, pampl_r1_cols);

        long m, n, k;

        // contract over k2
        // B2_(i2 k2) . ( W_([i1 k5 k4 k3] k2) )' -> T_(i2 [i1 k5 k4 k3])
        m = i2_size;
        n = pampl_r1_cols;
        k = k2_size_l;
        BLASW_DGEMM_NT(1.0, 0.0, m, n, k, booster_2_subm, booster_2->dims[0], w15j_tensor->d, n, pampl_tmp, m);

        // reuse 15j tensor memory
        w15j_tensor->d = dmatrix_realloc(w15j_tensor->d, i3_size, pampl_r2_cols);

        // contract over k3
        // B3_(i3 k3) . ( T_([i2 i1 k5 k4] k3) )' -> T_(i3 [i2 i1 k5 k4])
        m = i3_size;
        n = pampl_r2_cols;
        k = k3_size_l;
        BLASW_DGEMM_NT(1.0, 0.0, m, n, k, booster_3_subm, booster_3->dims[0], pampl_tmp, n, w15j_tensor->d, m);

        pampl_tmp = dmatrix_realloc(pampl_tmp, i4_size, pampl_r3_cols);

        // contract over k4
        // B4_(i4 k4) . ( T_([i3 i2 i1 k5] k4) )' -> T_(i4 [i3 i2 i1 k5])
        m = i4_size;
        n = pampl_r3_cols;
        k = k4_size_l;
        BLASW_DGEMM_NT(1.0, 0.0, m, n, k, booster_4_subm, booster_4->dims[0], w15j_tensor->d, n, pampl_tmp, m);

        // contract over k5
        // B5_(i5 k5) . ( T_([i4 i3 i2 i1] k5) )' -> T_(i5 [i4 i3 i2 i1])
        m = i5_size;
        n = pampl_cols;
        k = k5_size_l;
        BLASW_DGEMM_NT(1.0, 1.0, m, n, k, booster_5_subm, booster_5->dims[0], pampl_tmp, n, pampl, m);

        //////////////////////////////////////////////////////////////////

        matrix_free(pampl_tmp);
        TENSOR_FREE(w15j_tensor);

    } // lind

    // reduce over all threads

    #ifdef USE_OMP
    #pragma omp critical (ls_reduction)
    #endif
    for (size_t s = 0; s < tot_size; s++) {
        t_range->d[s] += pampl[s];
    }

    vector_free(pampl);

    #ifdef USE_OMP
    } // omp parallel
    #endif

    TENSOR_FREE(sjA);
    TENSOR_FREE(sjB);
    TENSOR_FREE(sjC);
    TENSOR_FREE(sjD);
    TENSOR_FREE(sjE);
    TENSOR_FREE(booster_2);
    TENSOR_FREE(booster_3);
    TENSOR_FREE(booster_4);
    TENSOR_FREE(booster_5);
    vector_free(xdims);
    free(ls_todo);

    #ifdef USE_MPI

    // reduce tensors over all nodes to master
    if (mpi_rank == MPI_MASTER) {
        MPI_Reduce(MPI_IN_PLACE, t_range->d, tot_size, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(t_range->d, t_range->d, tot_size, MPI_DOUBLE, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    }

    #endif

MPI_MASTERONLY_START

    // accumulate with lower Dl tensor if it was found
    if (found) {

        if (t_range->dim != t_found->dim)
            error("found tensor and computed tensor have different sizes");

        for (size_t s = 0; s < tot_size; s++) {
            t_range->d[s] += t_found->d[s];
        } 

    }

    TOC();
    verb(SL2CFOAM_VERBOSE_HIGH, "... done in %.3f seconds.\n", ELAPSED());

    // add tag
    if (t_range != NULL) {

        sl2cfoam_vertex_tag tag;

        tag.Immirzi = IMMIRZI;
        memcpy(tag.two_js, two_js, 10 * sizeof(dspin));
        tag.two_i1_range[0] = two_i1_min;
        tag.two_i1_range[1] = two_i1_max;
        tag.two_i2_range[0] = two_i2_min;
        tag.two_i2_range[1] = two_i2_max;
        tag.two_i3_range[0] = two_i3_min;
        tag.two_i3_range[1] = two_i3_max;
        tag.two_i4_range[0] = two_i4_min;
        tag.two_i4_range[1] = two_i4_max;
        tag.two_i5_range[0] = two_i5_min;
        tag.two_i5_range[1] = two_i5_max;
        tag.Dl = Dl;

        TENSOR_TAG(t_range, &tag, sizeof(tag));

    }

#ifndef NO_IO

    // store if configured
    if (bitmask_has(tresult, TENSOR_RESULT_STORE)) {
        
        char path[strlen(DIR_AMPLS) + 512];
        build_vertex_path(path, two_js, 
                                two_i1_min, two_i1_max, 
                                two_i2_min, two_i2_max, 
                                two_i3_min, two_i3_max, 
                                two_i4_min, two_i4_max, 
                                two_i5_min, two_i5_max,
                                Dl);

        TENSOR_STORE(t_range, path);

    }

#endif

MPI_MASTERONLY_END

tensor_range_return:

MPI_MASTERONLY_START

    // return the result if configured (and master)
    if (bitmask_has(tresult, TENSOR_RESULT_RETURN))
        return t_range;

MPI_MASTERONLY_END
        
    if (t_range != NULL) TENSOR_FREE(t_range);
    return NULL;

}

sl2cfoam_tensor_vertex* sl2cfoam_vertex_range_load(sl2cfoam_dspin two_js[10],
                                                   sl2cfoam_dspin two_i1_min, sl2cfoam_dspin two_i1_max, 
                                                   sl2cfoam_dspin two_i2_min, sl2cfoam_dspin two_i2_max, 
                                                   sl2cfoam_dspin two_i3_min, sl2cfoam_dspin two_i3_max, 
                                                   sl2cfoam_dspin two_i4_min, sl2cfoam_dspin two_i4_max, 
                                                   sl2cfoam_dspin two_i5_min, sl2cfoam_dspin two_i5_max, 
                                                   int Dl) {

    
    char path[strlen(DIR_AMPLS) + 512];
    build_vertex_path(path, two_js, 
                            two_i1_min, two_i1_max, 
                            two_i2_min, two_i2_max, 
                            two_i3_min, two_i3_max, 
                            two_i4_min, two_i4_max, 
                            two_i5_min, two_i5_max,
                            Dl);

    return sl2cfoam_vertex_load(path);

}

sl2cfoam_tensor_vertex* sl2cfoam_vertex_fullrange_load(sl2cfoam_dspin two_js[10], int Dl) {

    char path[strlen(DIR_AMPLS) + 512];
    build_vertex_path_fullrange(path, two_js, Dl);

    return sl2cfoam_vertex_load(path);

}

sl2cfoam_tensor_vertex* sl2cfoam_vertex_load(char* path) {

    sl2cfoam_tensor_vertex* t;
    TENSOR_LOAD(vertex, t, 5, path);
    return t;

}

void sl2cfoam_vertex_free(sl2cfoam_tensor_vertex* t) {
    TENSOR_FREE(t);
}