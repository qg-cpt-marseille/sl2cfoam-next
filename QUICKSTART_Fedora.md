# Installation

This guide will help you install `sl2cfoam-next` on a Linux machine. All the steps have been tested on a new _Ubuntu 22.04 LTS_ installation.

## 1. Install the dependencies

We assume you have the following packages already installed on your computer `lzip, gcc, g++, m4, make`. If you miss any, you can use your package manager to install them. In the following, we will use as an example `apt` (Advanced package tool), which is the default package manager on _Debian_ based distributions like _Ubuntu 22.04 LTS_

```[bash]
sudo dnf install zip gcc g++ m4 make
```

**Step 1.** Install [libquadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) (the GCC Quad-Precision Math Library API) and [Open MP](https://www.openmp.org/) (support multi-platform shared-memory multiprocessing programming in C). Also, in this case, both should be available in your package manager.

```[bash]
sudo dnf install libquadmath libquadmath-devel libomp libomp-devel
```

**Step 2.** Install next is [GMP](https://gmplib.org) (The GNU Multiple Precision Arithmetic Library). Also in this case, we suggest using the version included in your package manager.

```[bash]
sudo dnf install gmp gmp-devel
```

**Step 3.** Install [MPFR](https://www.mpfr.org/) (a C library for multiple-precision floating-point computations with correct rounding). This library has to be downloaded and compiled from its source code. Change the version number of MPFR to match the most recent one. At the time of writing, the last stable version is 4.2.1

```[bash]
sudo dnf install mpfr mpfr-devel
```

**Step 4.** Install [MPC](http://www.multiprecision.org/mpc/) (a C high precision library for complex numbers built on GMP and MPFR). This library has to be downloaded and compiled from its source code. Change the version number of MPC to match the most recent one. At the time of writing, the last stable version is 1.3.1

```[bash]
sudo dnf install libmpc libmpc-devel
```

**Step 5.** Install [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library) and the [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) (Math Kernel Library) both are present in standard package managers therefore the installation should be simple.

```[bash]
sudo dnf install gsl gsl-devel
```
To install Intel MKL you need first to add the repo to the dnf manager. I followed the instruction on the [Intel MKL website](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2023-1/yum-dnf-zypper.html#GUID-01F72C0F-4297-49AE-ABB0-41709E2D9E2C). First create the file 
```[bash]
sudo nano /etc/yum.repos.d/oneAPI.repo
```
and paste this in it
```
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
```
then just update dnf and install the packages you need
```[bash]
sudo dnf update
sudo dnf install intel-oneapi-mkl intel-oneapi-mkl-devel
```

##  2. Install `sl2cfoam-next`

**Step 1.** First, clone the repository in a working directory

```[bash]
git clone https://github.com/qg-cpt-marseille/sl2cfoam-next.git
cd sl2cfoam-next
```

**Step 2.** Create the `ext` directory where you will download and compile [wigxjpf](http://fy.chalmers.se/subatom/wigxjpf) and [fastwigxj](http://fy.chalmers.se/subatom/fastwigxj) the libraries we use to compute SU(2) invariants. Create and enter the `ext` directory.

```[bash]
mkdir ext
cd ext
```

Second, download, untar, and compile the `wigxjpf` library. Notice how we rename the directory, removing the version number to comply with the `sl2cfoam-next` makefile.

```
wget http://fy.chalmers.se/subatom/wigxjpf/wigxjpf-1.13.tar.gz
mv wigxjpf-1.11.tar.gz wigxjpf.tar.gz
tar -xf wigxjpf.tar.gz
cd wigxjpf
make
cd ..
```

Last, download, untar, and compile the `fastwigxj` library. Notice how we rename the directory, removing the version number to comply with the `sl2cfoam-next` makefile.

```[bash]
wget http://fy.chalmers.se/subatom/fastwigxj/fastwigxj-1.4.1.tar.gz
mv fastwigxj-1.4.1.tar.gz fastwigxj.tar.gz
tar -xf fastwigxj.tar.gz
cd fastwigxj
```

If you are using a recent version of GCC there is a [problem](http://fy.chalmers.se/subatom/fastwigxj/CHANGELOG) with calculating the ${9j}$ symbols and `fastwigxj` refuses to compile. You do not worry about it, as `sl2cfoam-next` doesn't use them. To circumvent this problem, you need to edit the file `src/wigner9j_canonicalize.c` and comment the last block of code using `/* ... */`. 
Then, you can go back to the main `fastwigxj` directory and compile.

```[bash]
make
cd ..
```

**Step 3.** Create a data directory. In this directory, you will store all the symbols you create in the calculations. In the main `sl2cfoam-next` directory, create the `data_sl2cfoam` directory and enter it.

```[bash]
mkdir data_sl2cfoam
cd data_sl2cfoam
```

Precompute the $(3jm)$ and ${6j}$ symbols you will use in your calculations. Use `fastwigxj` to do it. Execute these two commands one at a time.

```[bash]
../ext/fastwigxj/bin/hash_js --max-E-3j=50 /dev/null ./table_50.3j
../ext/fastwigxj/bin/hash_js --max-E-6j=40 /dev/null ./table_40.6j
```

In the directory `data_sl2cfoam`, you should see two new files, `table_50.3j` and `table_40.6j`, respectively, containing all the $(3jm)$ symbols with spins up to 25 (50/2) and all the ${6j}$ symbols with spins up to 20 (40/2).

**Step 4.** Move back to the main directory `cd ..` and compile the library

```[bash]
make
```

If you see an error complaining that the system cannot find `mkl` please edit the `sl2cfoam-next` `MAKEFILE` and replace line 67 with the correct mkl path. On _Ubuntu 22.04 LTS_ that line should read

```[txt]
BLAS_CFLAGS = -DUSE_MKL -I/usr/include/mkl -DMKL_ILP64
```

**Step 5.** You are ready to use `sl2cfoam-next` directly in `C`. However, one of the main advantages of `sl2cfoam-next` is its `Julia` interface. You will set it up in the next section.

##  3. Install Julia

You can download and install Julia following the instructions from its [website](https://julialang.org/downloads/). Our advice is just to download the latest stable version and just untar it.

The only delicate additional step is to replace the `libmpfr` included with Julia. The version you installed at the beginning of this document. You should replace `julia-1.9.4/lib/julia/libmpfr.so.6.1.1` with `usr/local/lib/libmpfr.so.6.1.1`. Please change the version numbers of both Julia and `libmpfr` to your case. 
In our case, our system had `libmpfr.so.6.2.1` we just removed Julia's version of the library, replaced it with ours, and renamed our version to match `libmpfr.so.6.1.1`.

Next, you need to export the location of `sl2cfoam-next` to tell Julia where to find it. Just run the following two `export` commands

```[bash]
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/sl2cfoam-next/lib
export JULIA_LOAD_PATH="/path/to/sl2cfoam-next/julia/:$JULIA_LOAD_PATH"
```

Finally, you can run the Julia CLI

```[bash]
./bin/julia
```

and add the only required package.

```[julia]
import Pkg
Pkg.add("HalfIntegers")
```

If you want to use a [Jupyter notebook](https://jupyter.org/) (we strongly suggest it), install it and also add the following Julia package

```[julia]
Pkg.add("IJulia")
```

# First calculation

You are ready to use `sl2cfoam-next` to perform calculations. The simplest example is to compute the EPRL vertex amplitude with Immirzi parameter $\gamma = 1.2$ and truncation parameter $DL=5$. First, setup the library

```[julia]
using SL2Cfoam
Immirzi = 1.2
sl2c_data_folder = "/path/to/sl2cfoam-next/data_sl2cfoam"
sl2c_configuration = SL2Cfoam.Config(VerbosityOff, VeryHighAccuracy, 100, 0)
SL2Cfoam.cinit(sl2c_data_folder, Immirzi, sl2c_configuration)
```

Second, compute the vertex

```
Dl = 5
spins = j245, j125, j124, j145, j235, j234, j345, j123, j135, j134 = ones(10)
v = vertex_compute(spins, Dl);
```

You can access the values of the vertex amplitude in the variable `v.a`. Please refer to the paper [How-To compute EPRL spinfoam amplitudes](https://arxiv.org/abs/2202.04360) and references therein for a detailed guide on how to do calculations with the EPRL spinfoam model.
