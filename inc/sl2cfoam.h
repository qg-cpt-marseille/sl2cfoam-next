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

#ifndef __SL2CFOAM_LIB_H__
#define __SL2CFOAM_LIB_H__

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/

#include <complex.h>
#include <stdbool.h>

#include "sl2cfoam_tensors.h"

///////////////////////////////////////////////////////////
// Spin machine types.
//
// By convention, all spin arguments to functions are
// (positive) integers of the form 2 * spin.
// Nomenclature is as follows: a dspin variable that starts
// with 'two_' means that the variable holds 2 * spin value.
// Example: j == 0.5, two_j == 1.
///////////////////////////////////////////////////////////    

// Machine type for double-spins.
typedef int sl2cfoam_dspin;

// Machine type for half-integer spins.
typedef double sl2cfoam_spin;


////////////////////////////////////////////////////////////////////////
// Configuration options
////////////////////////////////////////////////////////////////////////

// Verbosity levels.
#define SL2CFOAM_VERBOSE_OFF   0
#define SL2CFOAM_VERBOSE_LOW   1
#define SL2CFOAM_VERBOSE_HIGH  2

// Accuracy levels.
#define SL2CFOAM_ACCURACY_NORMAL    0
#define SL2CFOAM_ACCURACY_HIGH      1
#define SL2CFOAM_ACCURACY_VERYHIGH  2

// Contains general parameters for setup of the library.
struct sl2cfoam_config {
    int             verbosity;               // verbosity level
    int             accuracy;                // accuracy level
    sl2cfoam_dspin  max_two_spin;            // maximum spin that library needs to manage
    size_t          max_MB_mem_per_thread;   // maximum amount of MBs used by a single program thread (0 = unlimited)
};

// Specifies how to handle a computed tensors (store or not, return or not...).
typedef enum sl2cfoam_tensor_result {
    TENSOR_RESULT_RETURN        = 1 << 0, // return the computed tensor
    TENSOR_RESULT_STORE         = 1 << 1, // store the computed tensor
    TENSOR_RESULT_STORE_BATCHES = 1 << 2  // store the computed batches (if any)
} sl2cfoam_tensor_result;

////////////////////////////////////////////////////////////////////////
// Global setup functions.
////////////////////////////////////////////////////////////////////////

// Call this function to initialize the library at the start 
// of the program with the chosen Barbero-Immirzi constant.
void sl2cfoam_init_conf(char* root_folder, double Immirzi, struct sl2cfoam_config* conf);

// Initalize with default configuration.
void sl2cfoam_init(char* root_folder, double Immirzi);

// Sets the verbosity at runtime.
void sl2cfoam_set_verbosity(int verbosity);

// Sets the accuracy at runtime.
void sl2cfoam_set_accuracy(int accuracy);

// Sets the Immirzi-Barbero parameter at runtime.
// WARNING: this function is not thread-safe.
void sl2cfoam_set_Immirzi(double Immirzi);

// Call this function when finished using the library
// at the end of the program.
void sl2cfoam_free();

// Returns true if the library was compiled 
// for MPI (flag MPI=...).
bool sl2cfoam_is_MPI();

// Returns the BLAS vendor set at compile time
// (flag BLAS=...).
const char* sl2cfoam_BLAS_vendor();

// Enables or disables internal parallelization with OpenMP at runtime.
// (Parallelization is ON by default if not disabled at compile time 
// with OMP=0; in this latter case this function has no effect)
void sl2cfoam_set_OMP(bool enable);

// Returns if internal parallelization with OpenMP is enabled at runtime.
bool sl2cfoam_get_OMP();

///////////////////////////////////////////////////////////////////////////
// Vertex amplitude functions.
///////////////////////////////////////////////////////////////////////////

// 4D vertex (4-simplex) tensor in the 5 intertwiner indices.
// NB: indices order is reversed -> (i5 i4 i3 i2 i1)
//     and memory layout is COLUMN-MAJOR.
// Macro defines type sl2cfoam_tensor_vertex.
TENSOR_INIT(vertex, 5);

// Additional (optional) infos added to vertex tensor object
// specifying the vertex properties.
typedef struct sl2cfoam_vertex_tag {
    double Immirzi;
    sl2cfoam_dspin two_js[10];
    sl2cfoam_dspin two_i1_range[2];
    sl2cfoam_dspin two_i2_range[2];
    sl2cfoam_dspin two_i3_range[2];
    sl2cfoam_dspin two_i4_range[2];
    sl2cfoam_dspin two_i5_range[2];
    int Dl;
} sl2cfoam_vertex_tag;

// Computes a single vertex amplitude given spins js, intertwiners is
// and number of shells Dl (0, 1, 2, ...).
double sl2cfoam_vertex_amplitude(sl2cfoam_dspin two_js[10], sl2cfoam_dspin two_is[5], int Dl);

// Computes a the vertex amplitudes for all intertwiners in the given range, given spins js,
// number of shells Dl and a parameter to specify how to handle the result.
sl2cfoam_tensor_vertex* sl2cfoam_vertex_range(sl2cfoam_dspin two_js[10],
                                              sl2cfoam_dspin two_i1_min, sl2cfoam_dspin two_i1_max, 
                                              sl2cfoam_dspin two_i2_min, sl2cfoam_dspin two_i2_max, 
                                              sl2cfoam_dspin two_i3_min, sl2cfoam_dspin two_i3_max, 
                                              sl2cfoam_dspin two_i4_min, sl2cfoam_dspin two_i4_max, 
                                              sl2cfoam_dspin two_i5_min, sl2cfoam_dspin two_i5_max, 
                                              int Dl, sl2cfoam_tensor_result tresult);

// Computes the vertex tensor for all allowed intertwiners.
sl2cfoam_tensor_vertex* sl2cfoam_vertex_fullrange(sl2cfoam_dspin two_js[10], int Dl, sl2cfoam_tensor_result tresult);

// Loads vertex tensors that have been computed and stored under root folder.
// NULL pointer is returned if there is an error loading the file.
sl2cfoam_tensor_vertex* sl2cfoam_vertex_range_load(sl2cfoam_dspin two_js[10],
                                                   sl2cfoam_dspin two_i1_min, sl2cfoam_dspin two_i1_max, 
                                                   sl2cfoam_dspin two_i2_min, sl2cfoam_dspin two_i2_max, 
                                                   sl2cfoam_dspin two_i3_min, sl2cfoam_dspin two_i3_max, 
                                                   sl2cfoam_dspin two_i4_min, sl2cfoam_dspin two_i4_max, 
                                                   sl2cfoam_dspin two_i5_min, sl2cfoam_dspin two_i5_max, 
                                                   int Dl);

sl2cfoam_tensor_vertex* sl2cfoam_vertex_fullrange_load(sl2cfoam_dspin two_js[10], int Dl);

// Loads any vertex tensor given the path to the file.
// NULL is returned is there is an error loading the file.
// This function can be called without initializing the library.
sl2cfoam_tensor_vertex* sl2cfoam_vertex_load(char* path);


///////////////////////////////////////////////////////////////////////////
// BF vertex amplitude functions.
///////////////////////////////////////////////////////////////////////////

// 4D BF vertex (15j) tensor in the 5 intertwiner indices.
// NB: indices order is reversed -> (i5 i4 i3 i2 i1)
//     and memory layout is COLUMN-MAJOR.
// Macro defines type sl2cfoam_tensor_vertex_BF.
TENSOR_INIT(vertex_BF, 5);

// Computes a the BF vertex amplitudes for all intertwiners in the given range, given spins js.
sl2cfoam_tensor_vertex_BF* sl2cfoam_vertex_BF_range(sl2cfoam_dspin two_js[10],
                                                    sl2cfoam_dspin two_i1_min, sl2cfoam_dspin two_i1_max, 
                                                    sl2cfoam_dspin two_i2_min, sl2cfoam_dspin two_i2_max, 
                                                    sl2cfoam_dspin two_i3_min, sl2cfoam_dspin two_i3_max, 
                                                    sl2cfoam_dspin two_i4_min, sl2cfoam_dspin two_i4_max, 
                                                    sl2cfoam_dspin two_i5_min, sl2cfoam_dspin two_i5_max);

// Computes the BF vertex tensor for all allowed intertwiners.
sl2cfoam_tensor_vertex_BF* sl2cfoam_vertex_BF_fullrange(sl2cfoam_dspin two_js[10]);


///////////////////////////////////////////////////////////////////////////
// Booster functions.
///////////////////////////////////////////////////////////////////////////

// Tensor for 4d booster coefficients.
// Indices: (i, k, l1, l2, l3, l4)
TENSOR_INIT(boosters, 6);

// Computes a boosters tensor.
// gf parameter is the gauge-fixed index (1 to 4).
// Spins order must match the order of the symbol (anti-clockwise).
// Set store parameter to true to store the tensor after computation.
sl2cfoam_tensor_boosters* sl2cfoam_boosters(int gf,
                                            sl2cfoam_dspin two_ja, sl2cfoam_dspin two_jb,
                                            sl2cfoam_dspin two_jc, sl2cfoam_dspin two_jd, 
                                            int Dl, bool store);

// Loads a computed tensor for the boosters given gauge-fixed index,
// spins and number of shells.
sl2cfoam_tensor_boosters* sl2cfoam_boosters_load(int gf,
                                                 sl2cfoam_dspin two_ja, sl2cfoam_dspin two_jb,
                                                 sl2cfoam_dspin two_jc,  sl2cfoam_dspin two_jd, 
                                                 int Dl);

// Computes the b4^gamma(j_a, l_a; i, k) coefficients for all possible 
// intertwiner pairs (i, k).
// Result matrix is stored with indices (i, k).
sl2cfoam_dmatrix sl2cfoam_b4(sl2cfoam_dspin two_j1, sl2cfoam_dspin two_j2, sl2cfoam_dspin two_j3, sl2cfoam_dspin two_j4,
                             sl2cfoam_dspin two_l1, sl2cfoam_dspin two_l2, sl2cfoam_dspin two_l3, sl2cfoam_dspin two_l4);

// Computes the b4^gamma(j_a, l_a; i, k) coefficients using adaptive integration
// for given range of intertwiners.
// The result is more accurate but very slow. Useful for testing the fast version.
// Result matrix is stored with indices (i, k).
sl2cfoam_dmatrix sl2cfoam_b4_accurate(sl2cfoam_dspin two_j1, sl2cfoam_dspin two_j2, sl2cfoam_dspin two_j3, sl2cfoam_dspin two_j4,
                                      sl2cfoam_dspin two_l1, sl2cfoam_dspin two_l2, sl2cfoam_dspin two_l3, sl2cfoam_dspin two_l4,
                                      sl2cfoam_dspin two_i_min, sl2cfoam_dspin two_i_max,
                                      sl2cfoam_dspin two_k_min, sl2cfoam_dspin two_k_max);


///////////////////////////////////////////////////////////////////////////
// Coherent states functions.
///////////////////////////////////////////////////////////////////////////

// For given spins and range of intertwiners, computes a (normalized) coherent state vector 
// corresponding to the spherical angles (theta, phi) of the four normals.
sl2cfoam_cvector sl2cfoam_coherentstate_range(sl2cfoam_dspin two_js[4],
                                              sl2cfoam_dspin two_i_min, sl2cfoam_dspin two_i_max,
                                              double theta_phi[4][2]);

// Computes a coherent state vector for all allowed intertwiners.
sl2cfoam_cvector sl2cfoam_coherentstate_fullrange(sl2cfoam_dspin two_js[4], 
                                                  double theta_phi[4][2]);


///////////////////////////////////////////////////////////////////////////
// Cleanup utilities.
///////////////////////////////////////////////////////////////////////////

// Frees a vertex tensor.
void sl2cfoam_vertex_free(sl2cfoam_tensor_vertex* t);

// Frees a BF vertex tensor.
void sl2cfoam_vertex_BF_free(sl2cfoam_tensor_vertex_BF* t);

// Frees a boosters tensor.
void sl2cfoam_boosters_free(sl2cfoam_tensor_boosters* t);

// Functions for freeing allocated vectors and matrices.
void sl2cfoam_vector_free(sl2cfoam_vector v);
void sl2cfoam_matrix_free(sl2cfoam_matrix m);


/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif /*__SL2CFOAM_LIB_H__*/
