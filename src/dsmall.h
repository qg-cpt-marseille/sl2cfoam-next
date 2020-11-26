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

#ifndef __SL2CFOAM_DSMALL_H__
#define __SL2CFOAM_DSMALL_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <quadmath.h>
#include <mpc.h>

#include "common.h"

////////////////////////////////////////////////////////////////
// Functions for computing the dsmall coefficients 
// using Collet's formula.
//
// TODO: only for rho != 0. The 0 case is left unimplemented.
////////////////////////////////////////////////////////////////

// Precomputes the Y coefficients.
void sl2cfoam_dsmall_Yc(mpc_ptr Ys[], int prec, double rho, 
                        dspin two_k, dspin two_j, dspin two_l, dspin two_p);

// Computes the dsmall coefficients.
// Input is a sequence (x) of N real numbers in (0 1),
// an array where to put the result, the precision and the
// Y coefficients for the sums over indices m and n.
// An optional array prefactor (can be NULL) may be given
// with precomputed prefactor independent of p for each r.
// 
// For each x = exp(-r) computes d^(rho, k)_jlp (r).
void sl2cfoam_dsmall(__complex128 ds[], __float128 xs[], size_t N,
                     int prec, mpc_ptr Ym[], mpc_ptr Yn[], mpc_ptr prefactor[],
                     double rho, dspin two_k, dspin two_j, dspin two_l, dspin two_p);


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_DSMALL_H__*/
