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

#ifndef __SL2CFOAM_TIMING_H__
#define __SL2CFOAM_TIMING_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <time.h>

////////////////////////////////////////////////////////////////////////
// Timing utils.
////////////////////////////////////////////////////////////////////////

static __thread int __tri;
static __thread struct timespec __tstart[10];
static __thread struct timespec __tstop[10];
#define TIC() { clock_gettime(CLOCK_MONOTONIC, &__tstart[__tri++]);}
#define TOC() { clock_gettime(CLOCK_MONOTONIC, &__tstop[--__tri]); }
#define ELAPSED() ((double)((__tstop[__tri].tv_sec - __tstart[__tri].tv_sec) + (__tstop[__tri].tv_nsec - __tstart[__tri].tv_nsec) / 1000000000.0 ))

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_TIMING_H__*/