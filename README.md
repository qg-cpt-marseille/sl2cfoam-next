# SL2Cfoam-next: computing LQG vertex amplitudes

**SL2Cfoam-next** is a library for computing vertex amplitudes of Loop Quantum Gravity. The library is optimized for computing the EPRL vertex amplitude [Engle et al., 2008] for all possible boundary intertwiners and fixed boundary spins. The amplitudes are stored in 5-dimensional arrays of double-precision numbers and can be manipulated using the provided C bindings or a convenient Julia module.

## Features

The computation uses the EPRL splitting introduced in [Speziale, 2017]. In particular, the library splits the amplitude in the computation of: _(i)_ SU(2) recoupling symbols _(ii)_ SL(2,C) dipole amplitudes (_boosters_) and _(iii)_ shelled sums. The library exports methods to compute:

- EPRL vertex amplitudes as tensors
- a single EPRL vertex amplitude
- Livine-Speziale coherent state coefficients
- booster coefficients as tensors

The Julia interface implements all these methods as well as additional methods for contracting vertex tensors with coherent states. Two additional tools are provided as standalone programs to compute a single amplitude and to store a full vertex tensor.

## Dependencies

The library depends on:

1. GNU GMP, MPFR, MPC
2. quadmath (GCC extension)
3. _wigxjpf_ and _fastwigxj_ [Johansson et al., 2015]
4. OpenMP
5. a BLAS implementation
6. OpenMPI (optional)
7. Julia >= 1.5 (optional)

## Compilation

The library can be compiled typing `make`. This compiles the shared library, the test programs and the library tools. There are additional flags that can be provided. For example:

- type `make BLAS=mkl/blasfeo/system` to choose between different BLAS libraries [Frison et al., 2018]
- type `make DEBUG=1` to build the debug version
- type `make MPI=1` to build the MPI version
- type `make OMP=0` to disable builtin OpenMP parallelization

The previous flags can be combined. The library searches for folders `wigxjpf`, `fastwigxj` and `blasfeo` (optional) under `ext/`. The variable `MKLROOT` must be set for compiling with MKL. The library has been tested with GCC version 8.1 or greater.

You can provide flags for GCC or the linker using the option `ADD_CFLAGS=...` after `make`. This can be used for example for providing a different definition of the Y-map. The default definition sets the irrep continuous label to `Immirzi * (j+1)`. Compiling with `ADD_CFLAGS=-DRHO_GJ` sets the Y-map to `Immirzi * j`.

## Usage (C interface)

See `inc/sl2cfoam.h`.

## Usage (Julia interface)

You must init the C library calling `SL2Cfoam.cinit(...)` providing:

1. a root folder with fastwigxj tables (one file with extension `.6j` and one with `.3j`)
2. the Immirzi parameter (can be later changed at runtime)
3. a `SL2Cfoam.Config` object with additional configuration

When finished please call `SL2Cfoam.cclear()`. To compute a vertex tensor call `vertex_compute(...)` providing a list of 10 spins, the number of shells (from 0 to 100) and optionally a range of intertwiners. You can load vertex tensors computed with the C library calling `vertex_load(...)`. You can compute coherent states coefficients calling `coherentstates_compute(...)`. You can contract a vertex tensor with 1 to 5 coherent states coefficients using `contract(v, cs...)`. See `julia/SL2Cfoam.jl` for more functions and details.

## License

SL2CFOAM-NEXT is free software: you can redistribute it and/or modify 
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

SL2CFOAM-NEXT is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public License
along with SL2CFOAM-NEXT. If not, see <http://www.gnu.org/licenses/>.

If you use the library you should cite the following paper: "Francesco Gozzini, *High-performance spinfoam numerics*, in preparation." 

