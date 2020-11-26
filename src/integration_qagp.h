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

#ifndef __SL2CFOAM_INTEGRATION_QAGP_H__
#define __SL2CFOAM_INTEGRATION_QAGP_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <complex.h>
#include <quadmath.h>

#include "common.h"

// Number of Gauss-Kronrod points, default rule: 30-61
#define GK_POINTS 61

////////////////////////////////////////////////////////////////
// Adaptive quadrature using 61-point gauss-konrod method.
// Inspired by GSL (QUADPACK) but I need also:
//    - evaluation in all points at once to exploit parallelism
//    - complex-integration with quad precision
//    - more robust (see Numerical Recipes sec 4.7)
////////////////////////////////////////////////////////////////

// Function to be integrated over (0 1) with parameters.
// There is an optional measure function (set to NULL if not defined).
typedef struct sl2cfoam_integration_function 
{
  void (*function) (__complex128 y[], __float128 x[], size_t N, void* params);
  void  (*measure) (__complex128 y[], __float128 x[], size_t N, void* params);
  void* params;
} sl2cfoam_integration_function;

// Return codes for adaptive integration.
typedef enum sl2cfoam_qagp_retcode
{
    QAGP_SUCCESS   = 0,      // integration ok
    QAGP_MAX_EVAL  = 1 << 0, // maximum number of evaluation reached
    QAGP_MAX_DEPTH = 1 << 1  // maximum depth in subdivision reached
} sl2cfoam_qagp_retcode;

// Integrates a complex function over (0 1) within estimated epsrel precision,
// stores result, estimated absolute error and number of evaluation.
// The function terminates and best approximation is returned if dx becomes less
// than dx_min or number of evaluations overcomes max_eval.
sl2cfoam_qagp_retcode sl2cfoam_qagp(sl2cfoam_integration_function* f, double epsrel, __complex128* result, __float128* abserr,
                                    size_t* neval, __float128 dx_min, size_t max_eval);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_INTEGRATION_QAGP_H__*/
