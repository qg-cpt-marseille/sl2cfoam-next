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

#ifndef __SL2CFOAM_BOOSTERS_H__
#define __SL2CFOAM_BOOSTERS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include "common.h"
#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"

////////////////////////////////////////////////////////////////
// Constructs the tensors with booster coefficients.
////////////////////////////////////////////////////////////////

// Precompute all the needed booster tensors for a vertex computation
// for given boundary spins.
void sl2cfoam_boosters_tensors_vertex(dspin two_js[10], int Dl,
                                      tensor_ptr(boosters)* b2, tensor_ptr(boosters)* b3,
                                      tensor_ptr(boosters)* b4, tensor_ptr(boosters)* b5,
                                      dspin b_two_i_mins[4]);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_BOOSTERS_H__*/
