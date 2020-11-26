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

///////////////////////////////////////////////////////////////
// Library interface and global variables.
///////////////////////////////////////////////////////////////

#include <mpfr.h>
#include <omp.h>

#include "common.h"
#include "sl2cfoam.h"
#include "sl2cfoam_tensors.h"
#include "utils.h"
#include "error.h"

// root folder
char* DATA_ROOT;

// global full configuration object
struct sl2cfoam_config CONFIG;

// global Immirzi parameter
double IMMIRZI;

// global verbosity level
int VERBOSITY;

// global accuracy level and default precision
int ACCURACY;
int MPBITS;

// flag to enable or disable OpenMP parallelization
bool OMP_PARALLELIZE;

void sl2cfoam_set_verbosity(int verbosity) {

    not_thread_safe();
    VERBOSITY = verbosity;
    CONFIG.verbosity = verbosity;

}

void sl2cfoam_set_accuracy(int accuracy) {

    not_thread_safe();
    ACCURACY = accuracy;
    CONFIG.accuracy = accuracy;

    switch (accuracy)
    {
    case SL2CFOAM_ACCURACY_NORMAL:
        MPBITS = 128;
        break;

    case SL2CFOAM_ACCURACY_HIGH:
        MPBITS = 256;
        break;

    case SL2CFOAM_ACCURACY_VERYHIGH:
        MPBITS = 512;
        break;
    
    default:
        error("wrong accuracy value");
    }

    mpfr_set_default_prec(MPBITS);
    mpfr_set_default_rounding_mode(MPFR_RNDN);

    // check if accuracy might be low;
    // normal accuracy looks OK for spins as high as 50 (and probably more);
    // for avg boundary spins of order j there can spins up
    // to 3*j in the 6j symbols
    // hence ...
    if (CONFIG.max_two_spin > 3 * (2 * 50) && accuracy < SL2CFOAM_ACCURACY_HIGH) {
        warning("accuracy might be too low for given maximum spin");
    }

}

void sl2cfoam_set_Immirzi(double Immirzi) {

    not_thread_safe();
    IMMIRZI = Immirzi;

}

bool sl2cfoam_is_MPI() {

#ifdef USE_MPI
    return true;
#else
    return false;
#endif

}

const char* sl2cfoam_BLAS_vendor() {

#ifdef USE_MKL
    return "MKL";
#elif  USE_BLASFEO
    return "BLASFEO";
#else
    return "SYSTEM";
#endif

}

void sl2cfoam_set_omp(bool enable) {
    OMP_PARALLELIZE = enable;
}

void sl2cfoam_vector_free(sl2cfoam_vector v) {
    vector_free(v);
}

void sl2cfoam_matrix_free(sl2cfoam_matrix m) {
    matrix_free(m);
}

