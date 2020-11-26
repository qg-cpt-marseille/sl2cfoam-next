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

#ifndef __SL2CFOAM_ERROR_H__
#define __SL2CFOAM_ERROR_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "mpi_utils.h"

////////////////////////////////////////////////////////////////////////
// Functions to handle exceptions.
////////////////////////////////////////////////////////////////////////

// Prints the file and line of the error location, an optional message
// and then interrupts program.
#ifdef USE_MPI
#define sl2cfoam_error(format, ...) {\
                int rank;  \
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                     \
                fprintf(stderr, "sl2cfoam (node %d): ERROR at file %s, line %d. Terminating all nodes...\n", rank, __FILE__, __LINE__);\
                fprintf(stderr, "sl2cfoam (node %d): ", rank);\
                fprintf(stderr, format, ##__VA_ARGS__);\
                fprintf(stderr, "\n");\
                MPI_Abort(MPI_COMM_WORLD, 1);\
                }
#else
#define sl2cfoam_error(format, ...) {\
                fprintf(stderr, "sl2cfoam: ERROR at file %s, line %d. Terminating...\n", __FILE__, __LINE__);\
                fprintf(stderr, "sl2cfoam: ");\
                fprintf(stderr, format, ##__VA_ARGS__);\
                fprintf(stderr, "\n");\
                exit(EXIT_FAILURE);\
                }
#endif

// Prints the file and line of the error location and an optional message.
#ifdef USE_MPI
#define sl2cfoam_warning(format, ...) {\
                int rank;  \
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                     \
                fprintf(stderr, "sl2cfoam (node %d): WARNING at file %s, line %d. Continuing...\n", rank, __FILE__, __LINE__);\
                fprintf(stderr, "sl2cfoam (node %d): ", rank);\
                fprintf(stderr, format, ##__VA_ARGS__);\
                fprintf(stderr, "\n");\
                }
#else
#define sl2cfoam_warning(format, ...) {\
                fprintf(stderr, "sl2cfoam: WARNING at file %s, line %d. Continuing...\n", __FILE__, __LINE__);\
                fprintf(stderr, "sl2cfoam: ");\
                fprintf(stderr, format, ##__VA_ARGS__);\
                fprintf(stderr, "\n");\
                }
#endif

// Raises an error when a GSL command failed.
#define sl2cfoam_check_gsl(status) \
        if (status) { sl2cfoam_error("GSL error code %d", status); }

// Fails if condition is not satisfied.
#define sl2cfoam_ensure(cond) \
        if (!(cond)) { sl2cfoam_error("failed check"); }

#ifdef SL2CFOAM_SHORT_NAMES
#define error(...) sl2cfoam_error(__VA_ARGS__)
#define warning(...) sl2cfoam_warning(__VA_ARGS__)
#define check_gsl(...) sl2cfoam_check_gsl(__VA_ARGS__)
#define ensure(...) sl2cfoam_ensure(__VA_ARGS__)
#endif


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__SL2CFOAM_ERROR_H__*/