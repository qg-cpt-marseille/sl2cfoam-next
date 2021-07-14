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

#ifndef __SL2CFOAM_TENSORS_H__
#define __SL2CFOAM_TENSORS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <sys/file.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include <complex.h>
#include <quadmath.h>
#include <string.h>
#include <errno.h>

////////////////////////////////////////////////////////////////////
// Routines for aligned memory allocation.
////////////////////////////////////////////////////////////////////

// align memory addresses to 64 bits
#define MEMORY_ALIGNMENT 64

// aligned memory allocation
// input is the number of elements
#ifdef USE_MKL
#include <mkl.h>
#define sl2cfoam_aligned_alloc(nel) mkl_malloc(nel * sizeof(double), MEMORY_ALIGNMENT)
#define sl2cfoam_aligned_alloc2(nbytes) mkl_malloc(nbytes, MEMORY_ALIGNMENT)
#define sl2cfoam_aligned_realloc(ptr, nel) mkl_realloc(ptr, nel * sizeof(double))
#define sl2cfoam_aligned_realloc2(ptr, nbytes) mkl_realloc(ptr, nbytes)
#define sl2cfoam_aligned_calloc(nel) mkl_calloc(nel, sizeof(double), MEMORY_ALIGNMENT)
#define sl2cfoam_aligned_calloc2(nbytes) mkl_calloc(nbytes, sizeof(char), MEMORY_ALIGNMENT)
#define sl2cfoam_aligned_free(buf) mkl_free(buf)
#else
// TODO: I get wrong results if I allocate memory with aligned_alloc
//        and ALIGNMENT > 16 bytes on my laptop... this is weird
//#define sl2cfoam_aligned_alloc(nel) aligned_alloc(MEMORY_ALIGNMENT, nel * sizeof(double))
#define sl2cfoam_aligned_alloc(nel) malloc(nel * sizeof(double))
#define sl2cfoam_aligned_alloc2(nbytes) malloc(nbytes)
#define sl2cfoam_aligned_realloc(ptr, nel) realloc(ptr, nel * sizeof(double))
#define sl2cfoam_aligned_realloc2(ptr, nbytes) realloc(ptr, nbytes)
#define sl2cfoam_aligned_calloc(nel) calloc(nel, sizeof(double))
#define sl2cfoam_aligned_calloc2(nbytes) calloc(nbytes, sizeof(char))
#define sl2cfoam_aligned_free(buf) free(buf)
#endif


////////////////////////////////////////////////////////////////////
// A set of macros to deal with multidimensional arrays of double-precision values.
// We call a multi-dimensional array a 'tensor'.
// Indices are integers starting from 0.
// Entries are stored in COLUMN-MAJOR ordering.
////////////////////////////////////////////////////////////////////

// Macros to set elements of an array of fixed length.
// nkeys > 10 not implemented
#define  TSET_1(a, v0)                                     a[0] = v0;
#define  TSET_2(a, v0, v1)                                 TSET_1(a, (v0)) a[1] = v1;
#define  TSET_3(a, v0, v1, v2)                             TSET_2(a, (v0), (v1)) a[2] = v2;
#define  TSET_4(a, v0, v1, v2, v3)                         TSET_3(a, (v0), (v1), (v2)) a[3] = v3;
#define  TSET_5(a, v0, v1, v2, v3, v4)                     TSET_4(a, (v0), (v1), (v2), (v3)) a[4] = v4;
#define  TSET_6(a, v0, v1, v2, v3, v4, v5)                 TSET_5(a, (v0), (v1), (v2), (v3), (v4)) a[5] = v5;
#define  TSET_7(a, v0, v1, v2, v3, v4, v5, v6)             TSET_6(a, (v0), (v1), (v2), (v3), (v4), (v5)) a[6] = v6;
#define  TSET_8(a, v0, v1, v2, v3, v4, v5, v6, v7)         TSET_7(a, (v0), (v1), (v2), (v3), (v4), (v5), (v6)) a[7] = v7;
#define  TSET_9(a, v0, v1, v2, v3, v4, v5, v6, v7, v8)     TSET_8(a, (v0), (v1), (v2), (v3), (v4), (v5), (v6), (v7)) a[8] = v8;
#define TSET_10(a, v0, v1, v2, v3, v4, v5, v6, v7, v8, v9) TSET_9(a, (v0), (v1), (v2), (v3), (v4), (v5), (v6), (v7), (v8)) a[9] = v9;

// Computes an index in COLUMN-MAJOR ORDERING.
// nkeys > 10 not implemented
#define  TIND_1(t, i0)                                     ( (i0) )
#define  TIND_2(t, i0, i1)                                 ( TIND_1(t, (i0)) + t->strides[1] * (i1) )
#define  TIND_3(t, i0, i1, i2)                             ( TIND_2(t, (i0), (i1)) + t->strides[2] * (i2) )
#define  TIND_4(t, i0, i1, i2, i3)                         ( TIND_3(t, (i0), (i1), (i2)) + t->strides[3] * (i3) )
#define  TIND_5(t, i0, i1, i2, i3, i4)                     ( TIND_4(t, (i0), (i1), (i2), (i3)) + t->strides[4] * (i4) )
#define  TIND_6(t, i0, i1, i2, i3, i4, i5)                 ( TIND_5(t, (i0), (i1), (i2), (i3), (i4)) + t->strides[5] * (i5) )
#define  TIND_7(t, i0, i1, i2, i3, i4, i5, i6)             ( TIND_6(t, (i0), (i1), (i2), (i3), (i4), (i5)) + t->strides[6] * (i6) )
#define  TIND_8(t, i0, i1, i2, i3, i4, i5, i6, i7)         ( TIND_7(t, (i0), (i1), (i2), (i3), (i4), (i5), (i6)) + t->strides[7] * (i7) )
#define  TIND_9(t, i0, i1, i2, i3, i4, i5, i6, i7, i8)     ( TIND_8(t, (i0), (i1), (i2), (i3), (i4), (i5), (i6), (i7)) + t->strides[8] * (i8) )
#define TIND_10(t, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9) ( TIND_9(t, (i0), (i1), (i2), (i3), (i4), (i5), (i6), (i7), (i8) + t->strides[9] * (i9) )

// C type of tensor
#define tensor_t(name) sl2cfoam_tensor_##name
#define tensor_ptr(name) sl2cfoam_tensor_##name *

#define __NUM_KEYS_MAX 128
#define __TAG_BYTES 1024

// Inits the tensor type given name and number of keys (indices).
//
// the fields are:
// num_keys: number of dimensions
// dims: list of dimensions
// dim: the (maximum) total number of elements in the array
// d: contains the data
// tag: contains optional infos (MAX 1024 bytes in size)
#define TENSOR_INIT(name, nkeys)          \
typedef struct sl2cfoam_tensor_##name {   \
    uint8_t num_keys;                     \
    size_t dims[nkeys];                   \
    ptrdiff_t strides[nkeys];             \
    size_t dim;                           \
    double* d;                            \
    void* tag;                            \
} sl2cfoam_tensor_##name;

// Creates a tensor given a pointer to an empty tensor,
// the number of keys and a list of dimensions.
// The tensor is initialized to 0.
#define TENSOR_CREATE(name, t, nkeys, ...)                   \
    {                                                        \
    if (nkeys > __NUM_KEYS_MAX)                              \
        error("too many indices")                            \
    t = malloc(sizeof(sl2cfoam_tensor_##name));              \
    t->num_keys = nkeys;                                     \
    TSET_##nkeys(t->dims, __VA_ARGS__ );                     \
    size_t dim = 1;                                          \
    for (int i = 0; i < nkeys; i++) {                        \
        t->strides[i] = (ptrdiff_t) dim;                     \
        dim *= (size_t) t->dims[i];                          \
    }                                                        \
    t->dim = dim;                                            \
    t->d = sl2cfoam_aligned_calloc(dim);                     \
    t->tag = calloc(__TAG_BYTES, sizeof(uint8_t));           \
    }

// Frees the memory of a tensor.
#define TENSOR_FREE(t)            \
    {                             \
    sl2cfoam_aligned_free(t->d);  \
    free(t->tag);                 \
    free(t);                      \
    }

// Fills a tensor with data from an array d.
#define TENSOR_FILL(t, arr) \
    memcpy(t->d, arr, t->dim * sizeof(double))

// Fills a tensor with zeros.
#define TENSOR_ZERO(t) \
    memset(t->d, 0, t->dim * sizeof(double))

// Computes the multidimensional index given a list of indices.
#define TENSOR_INDEX(t, nkeys, ...) \
    TIND_##nkeys(t, __VA_ARGS__)

// Gets a value from a tensor given a list of indices.
#define TENSOR_GET(t, nkeys, ...) \
    t->d[TENSOR_INDEX(t, nkeys, __VA_ARGS__)]

// Sets a value from a tensor given a list of indices.
#define TENSOR_SET(v, t, nkeys, ...) \
    t->d[TENSOR_INDEX(t, nkeys, __VA_ARGS__)] = v

// Writes the tag field with data from dat
// (previous tag if any is deleted).)
#define TENSOR_TAG(t, dat, nbytes)   \
    {                                \
    memset(t->tag, 0, __TAG_BYTES);  \
    memcpy(t->tag, dat, nbytes);     \
    }
    

///////////////////////////////////////////////////////////////
// File system macros.
///////////////////////////////////////////////////////////////

// Stores a tensor to disk (dimensions, data and tag).
// The header contain the dimensions and the tag.
#define TENSOR_STORE(t, path)                                                 \
    {                                                                         \
    size_t __dim_bytes = sizeof(size_t) * __NUM_KEYS_MAX;                     \
    size_t __header_bytes = __dim_bytes + __TAG_BYTES;                        \
    uint8_t* __tensor_header = calloc(__header_bytes, sizeof(uint8_t));       \
    memcpy(__tensor_header, t->dims, sizeof(size_t) * t->num_keys);           \
    memcpy(__tensor_header + __dim_bytes, t->tag, __TAG_BYTES);               \
    FILE *ptr;                                                                \
    ptr = fopen(path, "wbx");                                                 \
    if (ptr == NULL) {                                                        \
        if (errno != EEXIST)                                                  \
            warning("error writing to file %s: %s", path, strerror(errno));   \
        goto store_##t##_end;                                                 \
    }                                                                         \
    if (flock(fileno(ptr), LOCK_EX) != 0) {                                   \
        warning("error locking file for writing: %s", strerror(errno));       \
        goto store_##t##_end;                                                 \
    }                                                                         \
    size_t ret;                                                               \
    ret = fwrite(__tensor_header, sizeof(uint8_t), __header_bytes, ptr);      \
    if (ret != __header_bytes) {                                              \
        warning("error storing tensor header, err: %s", strerror(errno));     \
        goto store_##t##_end;                                                 \
    }                                                                         \
    ret = fwrite(t->d, sizeof(double), t->dim, ptr);                          \
    if (ret != t->dim) {                                                      \
        warning("error storing tensor data, err: %s", strerror(errno));       \
        goto store_##t##_end;                                                 \
    }                                                                         \
    fflush(ptr);                                                              \
    if (flock(fileno(ptr), LOCK_UN) != 0) {                                   \
        warning("error unlocking file: %s", strerror(errno));                 \
    }                                                                         \
store_##t##_end:                                                              \
    if (ptr != NULL) fclose(ptr);                                             \
    free(__tensor_header);                                                    \
    }

// Reads a tensor from disk.
// It must be later deallocated with TENSOR_FREE(t).
#define TENSOR_LOAD(name, t, nkeys, path)                                     \
    {                                                                         \
    FILE *ptr;                                                                \
    ptr = fopen(path,"rb");                                                   \
    if (ptr == NULL) {                                                        \
        warning("error opening file %s: %s", path, strerror(errno));          \
        t = NULL; goto load_##t##_end;                                        \
    }                                                                         \
    if (flock(fileno(ptr), LOCK_SH) != 0) {                                   \
        warning("error locking file for reading: %s", strerror(errno));       \
        goto load_##t##_end;                                                  \
    }                                                                         \
    t = malloc(sizeof(sl2cfoam_tensor_##name));                               \
    t->num_keys = nkeys;                                                      \
    t->tag = calloc(__TAG_BYTES, sizeof(uint8_t));                            \
    size_t ret;                                                               \
    ret = fread(t->dims, sizeof(size_t), nkeys, ptr);                         \
    fseek(ptr, sizeof(size_t) * __NUM_KEYS_MAX, SEEK_SET);                    \
    ret += fread(t->tag, sizeof(uint8_t), __TAG_BYTES, ptr);                  \
    if (ret != t->num_keys + __TAG_BYTES) {                                   \
        warning("error reading tensor header, err: %s", strerror(errno));     \
        free(t);                                                              \
        t = NULL; goto load_##t##_end;                                        \
    }                                                                         \
    size_t dim = 1;                                                           \
    for (int i = 0; i < nkeys; i++) {                                         \
        t->strides[i] = (ptrdiff_t) dim;                                      \
        dim *= (size_t) t->dims[i];                                           \
    }                                                                         \
    t->dim = (size_t) dim;                                                    \
    t->d = sl2cfoam_aligned_calloc(dim);                                      \
    ret = fread(t->d, sizeof(double), dim, ptr);                              \
    if (ret != t->dim) {                                                      \
        warning("error reading tensor data, err: %s", strerror(errno));       \
        sl2cfoam_aligned_free(t->d); free(t);                                 \
        t = NULL; goto load_##t##_end;                                        \
    }                                                                         \
    if (flock(fileno(ptr), LOCK_UN) != 0) {                                   \
        warning("error unlocking file: %s", strerror(errno));                 \
    }                                                                         \
load_##t##_end:                                                               \
    if (ptr != NULL) fclose(ptr);                                             \
    }


///////////////////////////////////////////////////////////////
// Simple custom types for matrices and vectors.
// Everything is initialized to 0 when allocated.
// NB: these matrices are NOT the same as 2d C-arrays
//     because of column-major ordering
///////////////////////////////////////////////////////////////

// 2d matrix in COLUMN-MAJOR format
typedef void*           sl2cfoam_matrix;
typedef double*         sl2cfoam_dmatrix;
typedef double complex* sl2cfoam_cmatrix;
typedef __float128*     sl2cfoam_qmatrix;
typedef __complex128*   sl2cfoam_zmatrix;

#define dmatrix_alloc(d1, d2) sl2cfoam_aligned_calloc(d1*d2)
#define dmatrix_realloc(m, d1, d2) sl2cfoam_aligned_realloc(m, d1*d2)

#define cmatrix_alloc(d1, d2) sl2cfoam_aligned_calloc2(d1*d2*sizeof(double complex))
#define qmatrix_alloc(d1, d2) sl2cfoam_aligned_calloc2(d1*d2*sizeof(__float128))
#define zmatrix_alloc(d1, d2) sl2cfoam_aligned_calloc2(d1*d2*sizeof(__complex128))

#define dmatrix_zero(m, d1, d2) memset(m, 0, d1 * d2 * sizeof(double))
#define cmatrix_zero(m, d1, d2) memset(m, 0, d1 * d2 * sizeof(double complex))
#define qmatrix_zero(m, d1, d2) memset(m, 0, d1 * d2 * sizeof(__float128))
#define zmatrix_zero(m, d1, d2) memset(m, 0, d1 * d2 * sizeof(__complex128))

#define matrix_get(m, d1, i, j) m[i + d1*j]
#define matrix_column(m, d1, j) (m + d1*j)
#define matrix_set(a, m, d1, i, j) (m[i + d1*j] = a)
#define matrix_free(m) sl2cfoam_aligned_free(m)

// vectors
typedef void*           sl2cfoam_vector;
typedef double*         sl2cfoam_dvector;
typedef double complex* sl2cfoam_cvector;
typedef __float128*     sl2cfoam_qvector;
typedef __complex128*   sl2cfoam_zvector;

#define dvector_alloc(d) sl2cfoam_aligned_calloc(d)
#define cvector_alloc(d) sl2cfoam_aligned_calloc2(d*sizeof(double complex))
#define qvector_alloc(d) sl2cfoam_aligned_calloc2(d*sizeof(__float128))
#define zvector_alloc(d) sl2cfoam_aligned_calloc2(d*sizeof(__complex128))

#define dvector_zero(m, d) memset(m, 0, d * sizeof(double))
#define cvector_zero(m, d) memset(m, 0, d * sizeof(double complex))
#define qvector_zero(m, d) memset(m, 0, d * sizeof(__float128))
#define zvector_zero(m, d) memset(m, 0, d * sizeof(__complex128))

#define vector_get(v, i) v[i]
#define vector_set(a, v, i) (v[i] = a)
#define vector_free(v) sl2cfoam_aligned_free(v)

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_TENSORS_H__*/