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
#include "jsymbols.h"
#include "verb.h"
#include "error.h"
#include "timing.h"

// computes the approximate memory used in a single thread
// returns in MBytes
static size_t max_mem_needed_per_thread_BF(dspin two_js[10], 
                                           dspin two_i1_min, dspin two_i1_max, 
                                           dspin two_i2_min, dspin two_i2_max, 
                                           dspin two_i3_min, dspin two_i3_max, 
                                           dspin two_i4_min, dspin two_i4_max, 
                                           dspin two_i5_min, dspin two_i5_max) {

    size_t i1_size, i2_size, i3_size, i4_size, i5_size;
    i1_size = DIV2(two_i1_max-two_i1_min) + 1;
    i2_size = DIV2(two_i2_max-two_i2_min) + 1;
    i3_size = DIV2(two_i3_max-two_i3_min) + 1;
    i4_size = DIV2(two_i4_max-two_i4_min) + 1;
    i5_size = DIV2(two_i5_max-two_i5_min) + 1;

    size_t tot_bytes = sizeof(double) * i1_size * i2_size * i3_size * i4_size * i5_size;

    size_t tot_MB = tot_bytes / (1 << 20); // integer truncated division
    return tot_MB;

}

sl2cfoam_tensor_vertex_BF* sl2cfoam_vertex_BF_fullrange(dspin two_js[10]) {

    MPI_FUNC_INIT();

    ///////////////////////////////////////////////////
    // NO MPI for BF vertex
    // master does the job, others return NULL
    ///////////////////////////////////////////////////

    #ifdef USE_MPI

    if (mpi_rank != MPI_MASTER) return NULL;

    #endif

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
    tensor_ptr(vertex_BF) t_full = NULL;

    verb(SL2CFOAM_VERBOSE_LOW, "Computing full BF tensor...\n");
    TIC();  

    // check if the max memory is capped
    // if yes the computation is batched
    bool batch_computation = false;

    size_t needed_MB;
    size_t max_MB = CONFIG.max_MB_mem_per_thread;
    if (max_MB > 0) {

        needed_MB = max_mem_needed_per_thread_BF(two_js,
                                                 two_i1_min_allowed, two_i1_max_allowed, 
                                                 two_i2_min_allowed, two_i2_max_allowed, 
                                                 two_i3_min_allowed, two_i3_max_allowed, 
                                                 two_i4_min_allowed, two_i4_max_allowed, 
                                                 two_i5_min_allowed, two_i5_max_allowed);

        // batch if memory is not enough (required memory is more than 90% available)
        if ((double)needed_MB > (0.9 * max_MB)) {

            batch_computation = true;

        }

    }

    if (batch_computation) {

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
        verb(SL2CFOAM_VERBOSE_LOW, "[ batching full tensor into %d small tensors ]\n", num_batches);

        dspin two_i1_min_batch, two_i1_max_batch;
        two_i1_min_batch = two_i1_min_allowed;

        TENSOR_CREATE(vertex_BF, t_full, 5, i5_size, i4_size, i3_size, i2_size, i1_size);

        tensor_ptr(vertex_BF) tbatch;
        size_t skip = 0;
        while (true) {

            two_i1_max_batch = min(two_i1_min_batch + 2*(i1_size_per_thread - 1), two_i1_max_allowed);
            
            verb(SL2CFOAM_VERBOSE_HIGH, "[ batching: computing two_i1 from %d to %d...]\n", two_i1_min_batch, two_i1_max_batch);
            tbatch = sl2cfoam_vertex_BF_range(two_js,
                                              two_i1_min_batch, two_i1_max_batch, 
                                              two_i2_min_allowed, two_i2_max_allowed, 
                                              two_i3_min_allowed, two_i3_max_allowed, 
                                              two_i4_min_allowed, two_i4_max_allowed, 
                                              two_i5_min_allowed, two_i5_max_allowed);

            // accumulate 
            memcpy(t_full->d + skip, tbatch->d, tbatch->dim * sizeof(double));
            skip += tbatch->dim;
            TENSOR_FREE(tbatch);      

            if (two_i1_max_batch == two_i1_max_allowed) break;

            two_i1_min_batch = two_i1_max_batch + 2;

        }

    } else {

        // do not batch
        t_full = sl2cfoam_vertex_BF_range(two_js,
                                          two_i1_min_allowed, two_i1_max_allowed, 
                                          two_i2_min_allowed, two_i2_max_allowed, 
                                          two_i3_min_allowed, two_i3_max_allowed, 
                                          two_i4_min_allowed, two_i4_max_allowed, 
                                          two_i5_min_allowed, two_i5_max_allowed);

    }

    TOC();
    verb(SL2CFOAM_VERBOSE_LOW, "... done in %.3f seconds.\n", ELAPSED());

    return t_full;

}

sl2cfoam_tensor_vertex_BF* sl2cfoam_vertex_BF_range(dspin two_js[10],
                                                    dspin two_i1_min, dspin two_i1_max, 
                                                    dspin two_i2_min, dspin two_i2_max, 
                                                    dspin two_i3_min, dspin two_i3_max, 
                                                    dspin two_i4_min, dspin two_i4_max, 
                                                    dspin two_i5_min, dspin two_i5_max) {

    MPI_FUNC_INIT();

    ///////////////////////////////////////////////////
    // NO MPI for BF vertex
    // master does the job, others return NULL
    ///////////////////////////////////////////////////

    #ifdef USE_MPI

    if (mpi_rank != MPI_MASTER) return NULL;

    #endif

    dspin two_j12, two_j13, two_j14, two_j15, two_j23,
          two_j24, two_j25, two_j34, two_j35, two_j45;

    ASSIGN_SPINS(two_js);

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
    tensor_ptr(vertex_BF) t_range = NULL;

    // compute dimensions
    size_t i1_size, i2_size, i3_size, i4_size, i5_size;
    i1_size = DIV2(two_i1_max-two_i1_min) + 1;
    i2_size = DIV2(two_i2_max-two_i2_min) + 1;
    i3_size = DIV2(two_i3_max-two_i3_min) + 1;
    i4_size = DIV2(two_i4_max-two_i4_min) + 1;
    i5_size = DIV2(two_i5_max-two_i5_min) + 1;

    // create result tensor
    TENSOR_CREATE(vertex_BF, t_range, 5, i5_size, i4_size, i3_size, i2_size, i1_size);

    ///////////////////////////////////////////////////////////////
    // compute 6js tensors
    ///////////////////////////////////////////////////////////////

    verb(SL2CFOAM_VERBOSE_HIGH, "Computing 6js tensors...\n");
    TIC();

    tensor_ptr(6j_4) sjA;
    tensor_ptr(6j_5) sjB;
    tensor_ptr(6j_6) sjC;
    tensor_ptr(6j_5) sjD;
    tensor_ptr(6j_4) sjE;

    dspin bounds[10];

    sl2cfoam_6j_tensors_vertex(two_js, two_i1_min, two_i1_max, 0, &sjA, &sjB, &sjC, &sjD, &sjE, bounds);

    TOC();
    verb(SL2CFOAM_VERBOSE_HIGH, "... done in %.3f seconds.\n", ELAPSED());

    dspin two_x_absmin = bounds[8];
    dspin two_x_absmax = bounds[9];

    // sometimes triangular inequalities leave nothing to compute
    // (e.g. in batching --- I am not sure this is relevant for BF ...)
    if (two_x_absmax < two_x_absmin) return t_range; // all zeros

    size_t x_size = DIV2(two_x_absmax - two_x_absmin) + 1;

    sl2cfoam_dvector xdims = dvector_alloc(x_size);
    for (dspin two_x = two_x_absmin; two_x <= two_x_absmax; two_x += 2) {
        xdims[DIV2(two_x-two_x_absmin)] = (double)DIM(two_x);
    }

    ///////////////////////////////////////////////////////////////
    // assemble 6js (to 15j)
    // start of parallel OMP code
    ///////////////////////////////////////////////////////////////

    verb(SL2CFOAM_VERBOSE_HIGH, "Assembling BF amplitudes...\n");
    TIC();
 
    dspin jj_sum = two_j12 + two_j13 + two_j14 + two_j15 + two_j23 +  
                   two_j24 + two_j25 + two_j34 + two_j35 + two_j45;

    // build 15j
    #ifdef USE_OMP
    #pragma omp parallel for collapse(2) if(OMP_PARALLELIZE)
    #endif
    for (dspin two_i2 = two_i2_min; two_i2 <= two_i2_max; two_i2 += 2) {
    for (dspin two_i3 = two_i3_min; two_i3 <= two_i3_max; two_i3 += 2) {

        int ii2 = DIV2(two_i2-two_i2_min);
        int ii3 = DIV2(two_i3-two_i3_min);
        const sl2cfoam_dvector w6j4 = sjD->d + TENSOR_INDEX(sjD, 5, 0, ii3, ii2, 0, 0);

    for (dspin two_i4 = two_i4_min; two_i4 <= two_i4_max; two_i4 += 2) {

        int ii4 = DIV2(two_i4-two_i4_min);
        const sl2cfoam_dvector w6j3 = sjC->d + TENSOR_INDEX(sjC, 6, 0, ii4, ii3, 0, 0, 0);

    for (dspin two_i5 = two_i5_min; two_i5 <= two_i5_max; two_i5 += 2) {

        int ii5 = DIV2(two_i5-two_i5_min);
        const sl2cfoam_dvector w6j2 = sjB->d + TENSOR_INDEX(sjB, 5, 0, ii5, ii4, 0, 0);

    for (dspin two_i1 = two_i1_min; two_i1 <= two_i1_max; two_i1 += 2) {

        int ii1 = DIV2(two_i1-two_i1_min);
        const sl2cfoam_dvector w6j1 = sjA->d + TENSOR_INDEX(sjA, 4, 0, ii1, ii5, 0);
        const sl2cfoam_dvector w6j5 = sjE->d + TENSOR_INDEX(sjE, 4, 0, ii1, ii2, 0);

        double w15j = 0.0;
        for (int i = 0; i < x_size; i++) {
            w15j += xdims[i] * w6j1[i] * w6j2[i] * w6j3[i] * w6j4[i] * w6j5[i];
        }

        if (w15j != 0) {
    
            // add sign
            w15j *= real_negpow(jj_sum + two_i1 + two_i2 + two_i3 + two_i4 + two_i5);

            // add dimensions
            w15j *= sqrt(DIM(two_i1) * DIM(two_i2) * DIM(two_i3) * DIM(two_i4) * DIM(two_i5));

            TENSOR_SET(w15j, t_range, 5, DIV2(two_i5-two_i5_min), 
                                         DIV2(two_i4-two_i4_min), DIV2(two_i3-two_i3_min), 
                                         DIV2(two_i2-two_i2_min), DIV2(two_i1-two_i1_min));

        }

    } // i2
    } // i3
    } // i4
    } // i5
    } // i1

    TENSOR_FREE(sjA);
    TENSOR_FREE(sjB);
    TENSOR_FREE(sjC);
    TENSOR_FREE(sjD);
    TENSOR_FREE(sjE);
    vector_free(xdims);

    TOC();
    verb(SL2CFOAM_VERBOSE_HIGH, "... done in %.3f seconds.\n", ELAPSED());

    return t_range;

}

void sl2cfoam_vertex_BF_free(sl2cfoam_tensor_vertex_BF* t) {
    TENSOR_FREE(t);
}
