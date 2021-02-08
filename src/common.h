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

#ifndef __SL2CFOAM_COMMON_H__
#define __SL2CFOAM_COMMON_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

///////////////////////////////////////////////////////////
// Common includes, macros and defines for the library.
///////////////////////////////////////////////////////////

#include <stddef.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <complex.h>

#include "sl2cfoam.h"
#include "config.h"

// This enables shortened names for non-exposed functions.
#define SL2CFOAM_SHORT_NAMES

// Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////
// Root folder. Set at library initialization.
///////////////////////////////////////////////////////////////

extern char* DATA_ROOT;

///////////////////////////////////////////////////////////////
// Barbero-Immirzi parameter and related folders.
///////////////////////////////////////////////////////////////

extern double IMMIRZI;
extern char* DIR_BOOSTERS;
extern char* DIR_AMPLS;

///////////////////////////////////////////////////////////////
// Accuracy parameter.
// This define also the default bits of precision.
///////////////////////////////////////////////////////////////

extern int ACCURACY;
extern int MPBITS;

///////////////////////////////////////////////////////////////
// Controls OMP parallelization. 
// Parallelization is ENABLED at library initialization.
///////////////////////////////////////////////////////////////

extern bool OMP_PARALLELIZE;

///////////////////////////////////////////////////////////////
// Global configuration object. Set at library initialization.
///////////////////////////////////////////////////////////////

extern struct sl2cfoam_config CONFIG;

///////////////////////////////////////////////////////////////
// Spin types (short names). 
///////////////////////////////////////////////////////////////

#ifdef SL2CFOAM_SHORT_NAMES
typedef sl2cfoam_dspin dspin;
typedef sl2cfoam_spin spin;
#endif

// Dimension of the spin-j representation.
#define DIM(two_j) ((long)(two_j) + 1)

// Converts a dspin to corresponding spin.
#define SPIN(two_j) ((spin)(two_j) * 0.5)

// Divide an INTEGER dspin by 2 (exact integer divison).
#define DIV2(two_j) ((int)((two_j) >> 1))

///////////////////////////////////////////////////////////////
// Spin labeling and recoupling
///////////////////////////////////////////////////////////////

#define ASSIGN_SPINS(two_js) \
    two_j12 = two_js[0]; \
    two_j13 = two_js[1]; \
    two_j14 = two_js[2]; \
    two_j15 = two_js[3]; \
    two_j23 = two_js[4]; \
    two_j24 = two_js[5]; \
    two_j25 = two_js[6]; \
    two_j34 = two_js[7]; \
    two_j35 = two_js[8]; \
    two_j45 = two_js[9];

#define ASSIGN_INTW_RANGES(suffix) \
    two_i1_min_##suffix  = (dspin) max(abs(two_j12-two_j13), abs(two_j14-two_j15)); \
    two_i1_max_##suffix  = (dspin) min(two_j12+two_j13, two_j14+two_j15);           \
    two_i2_min_##suffix  = (dspin) max(abs(two_j23-two_j24), abs(two_j25-two_j12)); \
    two_i2_max_##suffix  = (dspin) min(two_j23+two_j24, two_j25+two_j12);           \
    two_i3_min_##suffix  = (dspin) max(abs(two_j34-two_j35), abs(two_j13-two_j23)); \
    two_i3_max_##suffix  = (dspin) min(two_j34+two_j35, two_j13+two_j23);           \
    two_i4_min_##suffix  = (dspin) max(abs(two_j45-two_j14), abs(two_j24-two_j34)); \
    two_i4_max_##suffix  = (dspin) min(two_j45+two_j14, two_j24+two_j34);           \
    two_i5_min_##suffix  = (dspin) max(abs(two_j15-two_j25), abs(two_j35-two_j45)); \
    two_i5_max_##suffix  = (dspin) min(two_j15+two_j25, two_j35+two_j45);

#define CHECK_NULL_INTW(intw, err) \
    if (two_##intw##_min_allowed > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " has empty range"); \
        err = true; \
    } \

#define CHECK_ALLOWED_INTW(intw, err) \
if (two_##intw < two_##intw##_min_allowed || two_##intw > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " must be in [%d %d]", DIV2(two_##intw##_min_allowed), DIV2(two_##intw##_max_allowed)); \
        err = true; \
    } \

#define CHECK_ALLOWED_INTW_RANGE(intw, err) \
    if (two_##intw##_min < two_##intw##_min_allowed || two_##intw##_max > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " range must be in [%d %d]", DIV2(two_##intw##_min_allowed), DIV2(two_##intw##_max_allowed)); \
        err = true; \
    } \


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_COMMON_H__*/
