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

#ifndef __SL2CFOAM_BLAS_WRAPPER_H__
#define __SL2CFOAM_BLAS_WRAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdint.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////
// Wrappers to BLAS functions of different libraries.
////////////////////////////////////////////////////////////////////

// BLAS libraries
#ifdef USE_MKL
#include <mkl.h>
#elif USE_BLASFEO
#include <cblas.h>
#include <blasfeo.h>
#else
#include <cblas.h>
#endif

// DGEMM routines
// TODO: implement BLASFEO wrapper
#ifdef USE_BLASFEO
#define BLASW_DGEMM_NN(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    blas_dgemm("N", "N", &(int){m}, &(int){n}, &(int){k}, &(double){alpha}, A, &(int){lda}, B, &(int){ldb}, &(double){beta}, C, &(int){ldc});
#define BLASW_DGEMM_NT(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    blas_dgemm("N", "T", &(int){m}, &(int){n}, &(int){k}, &(double){alpha}, A, &(int){lda}, B, &(int){ldb}, &(double){beta}, C, &(int){ldc});
#define BLASW_DGEMM_TN(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    blas_dgemm("T", "N", &(int){m}, &(int){n}, &(int){k}, &(double){alpha}, A, &(int){lda}, B, &(int){ldb}, &(double){beta}, C, &(int){ldc});
#define BLASW_DGEMM_TT(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    blas_dgemm("T", "T", &(int){m}, &(int){n}, &(int){k}, &(double){alpha}, A, &(int){lda}, B, &(int){ldb}, &(double){beta}, C, &(int){ldc});
#else
#define BLASW_DGEMM_NN(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#define BLASW_DGEMM_NT(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#define BLASW_DGEMM_TN(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#define BLASW_DGEMM_TT(alpha, beta, m, n, k, A, lda, B, ldb, C, ldc) \
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                                                                               
#endif

// DGEMV routines
#ifdef USE_BLASFEO
// fallback to CBLAS
#define BLASW_DGEMV(m, n, a, x, y) cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, a, m, x, 1, 0.0, y, 1);
#else
#define BLASW_DGEMV(m, n, a, x, y) cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0, a, m, x, 1, 0.0, y, 1);
#endif

// DDOT routines
#ifdef USE_BLASFEO
#define BLASW_DDOT(r, d, v1, v2) r = blas_ddot(d, v1, &(int){1}, v2, &(int){1})
#else
#define BLASW_DDOT(r, d, v1, v2) r = cblas_ddot(d, v1, 1, v2, 1)
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_BLAS_WRAPPER_H__*/