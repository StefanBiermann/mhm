# -*- Makefile -*-

#
# Setup file for Gnu compiler 4.7 with easybuild/goolf/1.4.10 at eve.ufz.de
#
# LICENSE
#    This file is part of the UFZ makefile project.
#
#    The UFZ makefile project is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    The UFZ makefile project is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with the UFZ makefile project. If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2013 Matthias Cuntz

# The Makefile assumes the standard bin, include and lib directories
# i.e. if SOMEDIR = /path/to/lib, the make file will define the two dirs
#   SOMEINC ?= SOMEDIR/include
#   SOMELIB ?= SOMEDIR/lib
# Define subirectories if another structure

#ToDo: this is a horrible mix of intel and gnu and does therefore not work. This is due to mpif90 being a wrapper for intel
# Paths
GNUDIR := /gpfs/software/juwels/stages/2018a/software/GCCcore/5.5.0
GNULIB := $(GNUDIR)/lib64
#GNUBIN := $(GNUDIR)/bin
GNUBIN := /gpfs/software/juwels/stages/2018a/software/ifort/2018.2.199-GCC-5.5.0/compilers_and_libraries_2018.2.199/linux/bin

# Compiling
F90 := $(GNUBIN)/intel64/ifort
FC  := $(F90)
CC  := $(GNUBIN)/gcc
CPP := /usr/bin/cpp
# GNU Fortran version >= 4.4
ifeq ($(release),debug)
    F90FLAGS += -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized
    FCFLAGS  += -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized
    CFLAGS   += -pedantic -Wall -W -O -g -Wno-maybe-uninitialized
else
    F90FLAGS += -O3
    FCFLAGS  += -O3
    CFLAGS   += -O3
endif
F90FLAGS += -cpp -ffree-form -ffixed-line-length-132
FCFLAGS  += -ffixed-form -ffixed-line-length-132
CFLAGS   +=
MODFLAG  := -I# space significant # has been J once
DEFINES  += -DGFORTRAN -DgFortran

# Linking
LIBS  += -L$(GNULIB)
RPATH += -Wl,-rpath,$(GNULIB)
iLDPATH = $(GNUDIR)/lib/gcc/x86_64-unknown-linux-gnu/4.8.1:/usr/local/cloog/0.18.0-2/lib:/usr/local/isl/0.11.1-2/lib:/usr/local/mpc/1.0.1-3/lib:/usr/local/mpfr/3.1.2-2/lib:/usr/local/gmp/5.1.2-1/lib
ifneq ($(LDPATH),)
    LDPATH += :$(iLDPATH)
else
    LDPATH := $(iLDPATH)
endif

# IMSL
IMSLDIR :=

# MKL
INTEL  := /usr/local/intel/composerxe-2011.4.191
MKLDIR := $(INTEL)/mkl
MKLINC := $(MKLDIR)/include/intel64/lp64
MKLLIB := $(MKLDIR)/lib/intel64
INTELLIB := $(INTEL)/compiler/lib/intel64
MKL95DIR :=

# NETCDF
ifeq ($(netcdf),netcdf3)
    NCDIR := 
else
    ZLIB    := /gpfs/software/juwels/stages/2018a/software/zlib/1.2.11-GCCcore-7.3.0/lib
    SZLIB   := /gpfs/software/juwels/stages/2018a/software/Szip/2.1.1-GCCcore-7.3.0/lib
   # HDF5LIB := /gpfs/software/juwels/stages/2018a/software/HDF5/1.8.20-GCC-7.3.0-serial/lib
    HDF5LIB := /gpfs/software/juwels/stages/2018a/software/HDF5/1.10.1-ipsmpi-2018a/lib
   # HDF5LIB := /gpfs/software/juwels/stages/2018a/software/HDF5/1.10.1-gpsmpi-2018a/lib
   # HDF5LIB := /gpfs/software/juwels/stages/2018a/software/HDF5/1.8.20-iimpi-2018a/lib
   # NCDIR   := /gpfs/software/juwels/stages/2018a/software/netCDF/4.6.1-GCC-7.3.0-serial
    NCDIR   := /gpfs/software/juwels/stages/2018a/software/netCDF/4.6.1-ipsmpi-2018a
   # NCDIR   := /gpfs/software/juwels/stages/2018a/software/netCDF/4.6.1-gpsmpi-2018a/lib64
   # NCDIR   := /gpfs/software/juwels/stages/2018a/software/netCDF/4.6.1-iimpi-2018a
   # NCFDIR  := /gpfs/software/juwels/stages/2018a/software/netCDF-Fortran/4.4.4-GCC-7.3.0-serial
    NCFDIR  := /gpfs/software/juwels/stages/2018a/software/netCDF-Fortran/4.4.4-ipsmpi-2018a
   # NCFDIR  := /gpfs/software/juwels/stages/2018a/software/netCDF-Fortran/4.4.4-gpsmpi-2018a/lib
   # NCFDIR  := /gpfs/software/juwels/stages/2018a/software/netCDF-Fortran/4.4.4-iimpi-2018a
endif

# PROJ
PROJ4DIR := /usr/local/proj/4.8.0-2_gcc_4.8.1_CentOS6
FPROJDIR :=

# LAPACK
LAPACKDIR   := /usr/local/lapack/3.5.0-1_gcc_4.8.1_CentOS6
GFORTRANDIR := $(GNUDIR)
GFORTRANLIB := $(GNULIB)

# MPI
#MPIDIR := /gpfs/software/juwels/stages/2018a/software/psmpi/5.2.1-1-iccifort-2018.2.199-GCC-5.5.0
MPIDIR := /gpfs/software/juwels/stages/2018a/software/psmpi/5.2.1-1-iccifort-2018.2.199-GCC-5.5.0
#MPIDIR := /usr/local/openmpi/gcc/1.8.4-2
#/usr/local/openmpi/1.8.2-1_gcc_4.8.1_CentOS6

# Documentation
DOXYGENDIR := /usr/local/doxygen/1.8.7-1_gcc_4.8.1_CentOS6/bin
DOTDIR     := /usr/bin
TEXDIR     := /usr/local/texlive/2011/bin/x86_64-linux
PERLDIR    := /usr/bin
iiLDPATH := /usr/local/flex/2.5.39-1_gcc_4.8.1_CentOS6/lib:/usr/local/bison/3.0.2-1_gcc_4.8.1_CentOS6/lib
ifneq ($(LDPATH),)
    LDPATH += :$(iiLDPATH)
else
    LDPATH := $(iiLDPATH)
endif

# Check some dependices, e.g. IMSL needs intel11 on eve
ifneq (,$(findstring $(system),eve))
    ifneq (,$(findstring $(imsl),vendor imsl))
        ifneq ($(icompiler),intel11)
            $(error Error: IMSL needs intel11.0.075, set 'compiler=intel11')
        endif
        ifeq ($(imsl),vendor)
            ifeq (,$(findstring $(mkl),mkl mkl95))
                $(error Error: IMSL vendor needs MKL, set 'mkl=mkl' or 'mkl=mkl95')
            endif
        endif
    endif
endif
