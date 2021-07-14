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


#ifndef __SL2CFOAM_MPI_UTILS_H__
#define __SL2CFOAM_MPI_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <unistd.h>
#include <sys/stat.h>

/////////////////////////////////////////////////////////////////////////
// MPI utilities.
/////////////////////////////////////////////////////////////////////////

#define MPI_MASTER 0

#ifdef USE_MPI
#include <mpi.h>

// Checks if MPI is initialized and assigns the local variables
// mpi_size and mpi_rank.
#define MPI_FUNC_INIT()                                                           \
    int mpi_initialized;                                                          \
    MPI_Initialized(&mpi_initialized);                                            \
    if (!mpi_initialized) {                                                       \
        fprintf(stderr, "MPI must be initialized before using the library.");     \
        MPI_Abort(MPI_COMM_WORLD, 1);                                             \
    }                                                                             \
    int mpi_rank, mpi_size;                                                       \
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);                                     \
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);                                     

#define MPI_MASTERONLY_START if (mpi_rank == MPI_MASTER) {
#define MPI_MASTERONLY_END }
#define MPI_MASTERONLY_DO if (mpi_rank == MPI_MASTER)

#define TENSOR_BCAST(name, t, nkeys, ...)                                         \
    {                                                                             \
    if (mpi_rank != MPI_MASTER && t == NULL) {                                    \
        TENSOR_CREATE(name, t, nkeys, __VA_ARGS__);                               \
    }                                                                             \
    MPI_Bcast(t->d, t->dim, MPI_DOUBLE, MPI_MASTER, MPI_COMM_WORLD);              \
    }

#else

#define MPI_FUNC_INIT() {}
#define MPI_MASTERONLY_START
#define MPI_MASTERONLY_END
#define MPI_MASTERONLY_DO

#define TENSOR_BCAST(name, t, nkeys, ...) {}

#endif

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_MPI_UTILS_H__*/
