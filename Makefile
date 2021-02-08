QUIET ?= @

CMD_CFLAGS = -DUSE_OMP

# OMP is enabled by default
# add 'OMP=0' to make call to disable it
OMP = 1

ifeq ($(OMP), 0)
CMD_CFLAGS = 
endif

# MPI is disabled by default
# add 'MPI=1' to make call to enable it
MPI = 0

ifeq ($(MPI), 1)
CC = mpicc
CMD_CFLAGS += -DUSE_MPI
else
CC = gcc
endif

# IO (storing tensors) is enabled by default
# add 'IO=0' to make call to disable it
IO = 1

ifeq ($(IO), 0)
CMD_CFLAGS += -DNO_IO
endif

AR = ar

# folder names 
EXTDIR = ext
SRCDIR = src
OBJDIR = obj
INCDIR = inc
LIBDIR = lib
BINDIR = bin
TESTDIR = test
TOOLSDIR = tools

# folders for external libraries
WIGDIR = $(EXTDIR)/wigxjpf
FASTWIGDIR = $(EXTDIR)/fastwigxj
BLASFEODIR = $(EXTDIR)/blasfeo

###############################################################################################

# set the BLAS library: MKL or BLASFEO or system default
# default to MKL
# call 'make target BLAS=mkl/blasfeo/system' to set
MKL = mkl
BLASFEO = blasfeo
SYSBLAS = system
BLAS = $(MKL)

ifeq ($(BLAS), $(SYSBLAS))
BLAS_CFLAGS = -DUSE_SYSBLAS
BLAS_LDLIBS = -lopenblas -lblas -lpthread
else ifeq ($(BLAS), $(BLASFEO))
BLAS_CFLAGS = -DUSE_BLASFEO -I$(BLASFEODIR)/include -I$(BLASFEODIR)/netlib/cblas/include
BLAS_LDFLAG = -L$(BLASFEODIR)/lib 
BLAS_LDLIBS = -lblasfeo -lgfortran -lblas -lpthread
else ifeq ($(BLAS), $(MKL))
BLAS_CFLAGS = -DUSE_MKL -I$(MKLROOT)/include/ -DMKL_ILP64
BLAS_LDLIBS = -Wl,--no-as-needed -lmkl_rt -lpthread -ldl
else
$(error Wrong BLAS vendor selected)
endif

###############################################################################################

# call 'make target DEBUG=1' to compile in debug mode
# set additional flags for the compiler by passing var ADD_CFLAGS=...
ifeq ($(DEBUG), 1)
CFLAGS = -std=gnu11 -fopenmp -Wall -g -Og -fPIC $(ADD_CFLAGS) $(BLAS_CFLAGS) $(CMD_CFLAGS) \
         -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -I$(EXTDIR) -Iinc -Isrc -DDEBUG_ON=1  # debug
else
CFLAGS = -std=gnu11 -fopenmp -O3 -fPIC -march=native -fno-math-errno $(ADD_CFLAGS) $(BLAS_CFLAGS) $(CMD_CFLAGS) \
         -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -I$(EXTDIR) -Iinc -Isrc -DDEBUG_OFF=1 # optimized
endif

LDFLAGS = $(CFLAGS) -L$(WIGDIR)/lib/ -L$(FASTWIGDIR)/lib/ $(BLAS_LDFLAG) -Wl,--allow-shlib-undefined
LDLIBS = $(BLAS_LDLIBS) -lmpc -lmpfr -lgmp -lquadmath -lfastwigxj -lwigxjpf -lm 
ARFLAGS = rcs

###############################################################################################

INCS = src/utils.h src/jsymbols.h src/verb.h

_OBJS = sl2cfoam.o setup.o vertex.o boosters.o coherentstates.o dsmall.o \
        integration_gk.o integration_qagp.o b4.o b4_qagp.o cgamma.o

_TESTS = lib_test tensor_test dsmall_test integration_test b4_test

_TOOLS = vertex-fulltensor vertex-amplitude

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))
TOOLS = $(patsubst %,$(BINDIR)/%,$(_TOOLS))

# library/src object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -c -o $@ $< 

# test programs
$(BINDIR)/%: $(TESTDIR)/%.c $(OBJS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -o $@ -Iinc/ $< -Llib/ $(LDFLAGS) -lsl2cfoam $(LDLIBS)


###############################################################################################

# shared library
$(LIBDIR)/libsl2cfoam.so: $(OBJS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) -shared $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS) 
	
###############################################################################################

# tools
$(BINDIR)/%: $(TOOLSDIR)/%.c $(OBJS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -o $@ -Iinc/ $< -Llib/ $(LDFLAGS) -lsl2cfoam $(LDLIBS)

###############################################################################################

.DEFAULT_GOAL := default

default: all

all: lib tests tools

lib: $(LIBDIR)/libsl2cfoam.so

tests: lib $(TESTS)

tools: lib $(TOOLS)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
	
