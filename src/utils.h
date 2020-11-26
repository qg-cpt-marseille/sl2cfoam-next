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


#ifndef __SL2CFOAM_UTILS_H__
#define __SL2CFOAM_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <omp.h>

#include "common.h"
#include "error.h"

/////////////////////////////////////////////////////////////////////////
// Bitmask utilities.
/////////////////////////////////////////////////////////////////////////

#define bitmask_set(bm, v) (bm |= v)
#define bitmask_unset(bm, v) (bm &= ~v)
#define bitmask_has(bm, v) ((bm & v) != 0)

/////////////////////////////////////////////////////////////////////////
// Various utilities.
/////////////////////////////////////////////////////////////////////////

// Checks if files exists.
static inline int file_exist(char* filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

// Returns true if obj is in list.
static inline bool in_list(void* list, size_t list_size, 
                           void* obj, size_t obj_size) {

    void* list_pt;
    for (int i = 0; i < list_size; i++) {
        list_pt = list + i * obj_size;
        if (memcmp(list_pt, obj, obj_size) == 0) {
            return true;
        }
    }

    return false;

}

// Warns the user if calling a non thread-safe function 
// inside an active parallel region.
static inline void not_thread_safe() {
    if (omp_in_parallel()) {
        warning("calling this function from a parallel region may result in unexpected behavior")
    }
}

/////////////////////////////////////////////////////////////////////////
// Min/max utilities.
/////////////////////////////////////////////////////////////////////////

// Min/max macros.
#define max(two_j1, two_j2)                  sl2cfoam_imax(two_j1, two_j2)
#define max3(two_j1, two_j2, two_j3)         max(max(two_j1, two_j2), two_j3)
#define max4(two_j1, two_j2, two_j3, two_j4) max(max3(two_j1, two_j2, two_j3), two_j4)

#define min(two_j1, two_j2)                  sl2cfoam_imin(two_j1, two_j2)
#define min3(two_j1, two_j2, two_j3)         min(min(two_j1, two_j2), two_j3)
#define min4(two_j1, two_j2, two_j3, two_j4) min(min3(two_j1, two_j2, two_j3), two_j4)

// Returns the maximum of two integers.
static inline int sl2cfoam_imax(int n1, int n2) {

	if (n1 > n2) {
		return n1;
	}
	return n2;

}

// Returns the minimum of two integers.
static inline int sl2cfoam_imin(int n1, int n2) {

	if (n1 < n2) {
		return n1;
	}
	return n2;
	
}


/////////////////////////////////////////////////////////////////////////
// Spin utilities.
/////////////////////////////////////////////////////////////////////////

// Returns true if spin is integer.
#define is_integer(two_j) (((two_j) % 2) == 0)

// Returns true if spin is semi-integer.
#define is_semi_integer(two_j) (((two_j) % 2) == 1)

// Checks that the dspin value corresponds to an integer spin.
#define ensure_integer_spin(two_j) \
    { if ((two_j) % 2 != 0) { error("integer check failed"); } }

// Absolute bounds for intertwiner dependent on l indices.
// gf is the gauge-fixed index (1 to 4) 
static inline void find_k_absolute_bounds(size_t* k_size, dspin* two_k_absmin, dspin* two_k_absmax, 
                                          dspin two_ja, dspin two_jb, dspin two_jc, dspin two_jd,
                                          dspin two_Dl, int gf) {

    dspin two_zero = is_integer(two_ja + two_jb) ? 0 : 1;

    switch (gf)
    {

    case 1:
        *two_k_absmin = max4(two_zero, two_ja-two_jb-two_Dl, two_jb-two_ja, abs(two_jc-two_jd)-two_Dl);
        *two_k_absmax = min(two_ja + two_jb + two_Dl, two_jc + two_jd + 2*two_Dl);
        break;

    case 2:
        *two_k_absmin = max4(two_zero, two_ja-two_jb, two_jb-two_ja-two_Dl, abs(two_jc-two_jd)-two_Dl);
        *two_k_absmax = min(two_ja + two_jb + two_Dl, two_jc + two_jd + 2*two_Dl);
        break;

    case 3:
        *two_k_absmin = max4(two_zero, abs(two_ja-two_jb)-two_Dl, two_jc-two_jd-two_Dl, two_jd-two_jc);
        *two_k_absmax = min(two_ja + two_jb + 2*two_Dl, two_jc + two_jd + two_Dl);
        break;

    case 4:
        *two_k_absmin = max4(two_zero, abs(two_ja-two_jb)-two_Dl, two_jc-two_jd, two_jd-two_jc-two_Dl);
        *two_k_absmax = min(two_ja + two_jb + 2*two_Dl, two_jc + two_jd + two_Dl);
        break;
    
    default:
        error("wrong gf index");

    }

    if (*two_k_absmax < *two_k_absmin) {
        *k_size = 0;
    } else {
        *k_size = DIV2(*two_k_absmax - *two_k_absmin) + 1;
    }

}


/////////////////////////////////////////////////////////////////////////
// Math utilities.
/////////////////////////////////////////////////////////////////////////

// Compensated summation.
#define comp_sum(inp, sum, c, y, t) \
	{ y = inp - c; t = sum + y; c = (t - sum) - y; sum = t; }

// Squaring macro.
#define SQ(d) ((d)*(d))

// Cube macro.
#define CUBE(d) ((d)*(d)*(d))

// Computes (-1)^j for generic spin j.
static inline double complex sl2cfoam_complex_negpow(dspin two_j) {

    int k = two_j % 2; // i factor
    int j = two_j / 2;

    if (k == 1) {
        if (j % 2 == 0) {
            return I;
        }
        return -I;
    }

    if (k == -1) {
        if (j % 2 == 0) {
            return -I;
        }
        return I;
    }

    if (j % 2 == 0) {
        return 1;
    }
    return -1;

}

// Computes (-1)^j for INTEGER spin j.
// Throws an error if j is not integer.
static inline int sl2cfoam_real_negpow(dspin two_j) {

    ensure_integer_spin(two_j);

    if (two_j % 4 == 0) {
        return 1;
    }
    return -1;

}


#ifdef SL2CFOAM_SHORT_NAMES
#define real_negpow(...) sl2cfoam_real_negpow(__VA_ARGS__)
#define complex_negpow(...) sl2cfoam_complex_negpow(__VA_ARGS__)
#define imin(...) sl2cfoam_imin(__VA_ARGS__)
#define imax(...) sl2cfoam_imax(__VA_ARGS__)
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_UTILS_H__*/
