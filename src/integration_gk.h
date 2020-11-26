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

#ifndef __SL2CFOAM_INTEGRATION_GK_H__
#define __SL2CFOAM_INTEGRATION_GK_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <quadmath.h>

#include "common.h"


////////////////////////////////////////////////////////////////
// Fixed quadrature using 61-point gauss-konrod method on 
// subintervals of (0 1) defined by a grid.
////////////////////////////////////////////////////////////////

// Number of Gauss-Kronrod points, default rule: 30-61
#define GK_POINTS 61

// Computes a grid in (0 1) with "harmonic progression"
// that is it divides (0 1) in N intervals with 2nd interval
// two times the first, third interval three times the first 
// and so on...
// Useful here since the dsmall functions oscillate faster and faster
// approaching 0 as rho increases (but only for high magnetic index p).
__float128* sl2cfoam_grid_harmonic(int intervals);

// Computes a uniform grid in (0 1).
__float128* sl2cfoam_grid_uniform(int intervals);

// Computes all the Gauss-Kronrod abscissae for the subintervals
// of (0 1) defined by the given grid.
// The number of points is GK_POINTS * intervals
// Default rule uses 61 GK points.
__float128* sl2cfoam_gk_grid_abscissae(int intervals, __float128* grid);

// Computes all the Gauss-Kronrod weights for the subintervals
// of (0 1) defined by the given grid.
double* sl2cfoam_gk_grid_weights(int intervals, __float128* grid);

// Computes all the (plain) Gauss weights  for the subintervals
// of (0 1) defined by the given grid.
double* sl2cfoam_gk_grid_weights_gauss(int intervals, __float128* grid);

// Computes the GK sum on a given grid.
// Inputs are function evaluations, measures and weigths.
// Returns the result of Kronrod extension and optionally
// an absolute error given by the difference with the Gauss sum
// (for this the wgs parameter must be not NULL);
double sl2cfoam_gk_grid(int intervals, double* ys, double* ms, 
                        double* wgks, double* wgs, double* abserr);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_INTEGRATION_GK_H__*/
