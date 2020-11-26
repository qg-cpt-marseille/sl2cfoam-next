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

/**********************************************************************/

#ifndef __SL2CFOAM_CGAMMA_H__
#define __SL2CFOAM_CGAMMA_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <mpfr.h>
#include <mpc.h>

#include "common.h"

////////////////////////////////////////////////////////////////////////
// Computes the logarithm of complex Gamma function
// in arbitrary precision (ab) using Lanczos methos.
////////////////////////////////////////////////////////////////////////

// Holds the coefficients for the Lanczos expansion.
typedef struct cgamma_lanczos {
    int n;
    double g;
    mpfr_ptr* c;
} cgamma_lanczos;

// Computes the (logarithm of) complex Gamma function of op and stores
// result in rop. Last argument must be precomputed Lanczos coefficients.
int sl2cfoam_complex_lngamma(mpc_t rop, mpc_t op, cgamma_lanczos* lanczos);
int sl2cfoam_complex_gamma(mpc_t rop, mpc_t op, cgamma_lanczos* lanczos);

// Computes the coefficients for the Lanczos expansion.
// Always call this before using complex gamma.
void sl2cfoam_cgamma_lanczos_fill(cgamma_lanczos* lanczos);

// Release memory of the Lanczos coefficients.
void sl2cfoam_cgamma_lanczos_free(cgamma_lanczos* lanczos);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_CGAMMA_H__*/
