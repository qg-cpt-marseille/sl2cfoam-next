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

#ifndef __SL2CFOAM_VERB_H__
#define __SL2CFOAM_VERB_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdarg.h>
#include <stdio.h>

#include "mpi_utils.h"

////////////////////////////////////////////////////////////////////////
// Functions for verbose outputs.
////////////////////////////////////////////////////////////////////////

extern int VERBOSITY;

// Verbose output. First param is verbosity value for
// the message. If it is higher than verbosity level
// then nothing is printed.
static inline void sl2cfoam_verb(int v, const char* restrict format, ...) {

    if (v > VERBOSITY) {
        return;
    }

    #ifdef USE_MPI

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    printf("sl2cfoam (node %d): ", mpi_rank);

    #else 

    printf("sl2cfoam: ");

    #endif

    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);

}

#ifdef SL2CFOAM_SHORT_NAMES
#define verb(...) sl2cfoam_verb(__VA_ARGS__)
#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_VERB_H__*/