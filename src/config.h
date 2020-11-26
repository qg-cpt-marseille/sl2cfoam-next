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

#ifndef __SL2CFOAM_CONFIG_H__
#define __SL2CFOAM_CONFIG_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

// Compile-time configuration.

// Y-map irrep reduction
#ifdef RHO_GJ
#define RHO(j) (IMMIRZI * (j))         // rho = gamma * j
#elif RHO_GJP1
#define RHO(j) (IMMIRZI * ((j) + 1.0)) // rho = gamma * ( j + 1 )
#else
#define RHO_GJP1 
#define RHO(j) (IMMIRZI * ((j) + 1.0)) // default [ gamma * ( j + 1 ) ]
#endif


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_CONFIG_H__*/