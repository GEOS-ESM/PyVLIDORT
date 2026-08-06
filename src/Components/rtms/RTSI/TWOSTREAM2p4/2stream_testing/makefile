#
# This file serves as the head Makefile for 2stream and its test programs
#

#
# Define some variables
#

UTIL_PATH = util
2S_MAIN_PATH = ../2stream_sourcecode
2S_TEST_PATH = .

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o


#
# Define default shell make will use
#

SHELL = /bin/bash


#
# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)
#

#  Intel
#FC = ifort

#  gfortran
#FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95


#
# Define FORTRAN90 compiler flags
#

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback
	FFLAGS_OPT = -O3
	FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -fbacktrace
	FFLAGS_DEBUG = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
	FFLAGS_DEBUG = -g -O2 -Wall -Wuninitialized -fbounds-check
	FFLAGS_OPT = -O3
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# For debug build, use "make DEBUG=t"
 ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
 endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif


#
# Define list of source files
# (Note: ordering is important because of dependencies)
#

# Standard sources

SOURCES =
SOURCES +=   \
   $(2S_MAIN_PATH)/2stream_Taylor.f90       \
   $(2S_MAIN_PATH)/2stream_inputs.f90       \
   $(2S_MAIN_PATH)/2stream_writemodules.f90 \
   $(2S_MAIN_PATH)/2stream_miscsetups.f90   \
   $(2S_MAIN_PATH)/2stream_solutions.f90    \
   $(2S_MAIN_PATH)/2stream_bvproblem.f90    \
   $(2S_MAIN_PATH)/2stream_thermalsup.f90   \
   $(2S_MAIN_PATH)/2stream_intensity.f90 

SOURCES_SUP = 
SOURCES_SUP +=   \
   $(2S_MAIN_PATH)/2stream_brdf_kernels.f90      \
   $(2S_MAIN_PATH)/2stream_sleave_routines.f90   \
   $(2S_MAIN_PATH)/2stream_brdf_supplement.f90   \
   $(2S_MAIN_PATH)/2stream_sleave_supplement.f90

# Linearized sources

SOURCES_LSSUP = 
SOURCES_LSSUP +=   \
   $(2S_MAIN_PATH)/2stream_ls_brdf_kernels.f90      \
   $(2S_MAIN_PATH)/2stream_ls_brdf_supplement.f90   \
   $(2S_MAIN_PATH)/2stream_lssl_sleave_supplement.f90

SOURCES_LA =
SOURCES_LA +=   \
   $(2S_MAIN_PATH)/2stream_l_inputs.f90       \
   $(2S_MAIN_PATH)/2stream_l_writemodules.f90 \
   $(2S_MAIN_PATH)/2stream_la_solutions.f90   \
   $(2S_MAIN_PATH)/2stream_l_thermalsup.f90

SOURCES_LS =
SOURCES_LS +=   \
   $(2S_MAIN_PATH)/2stream_ls_bvproblem.f90   \
   $(2S_MAIN_PATH)/2stream_ls_jacobians.f90   \
   $(2S_MAIN_PATH)/2stream_lssl_jacobians.f90

SOURCES_LC =
SOURCES_LC +=   \
   $(2S_MAIN_PATH)/2stream_lc_miscsetups.f90  \
   $(2S_MAIN_PATH)/2stream_lc_bvproblem.f90   \
   $(2S_MAIN_PATH)/2stream_lc_jacobians.f90

SOURCES_LP =
SOURCES_LP +=   \
   $(2S_MAIN_PATH)/2stream_lp_miscsetups.f90  \
   $(2S_MAIN_PATH)/2stream_lp_bvproblem.f90   \
   $(2S_MAIN_PATH)/2stream_lp_jacobians.f90

# Main tests

SOURCES_MASTER  = 
SOURCES_MASTER +=   \
   $(SOURCES) $(SOURCES_SUP)               \
   $(2S_MAIN_PATH)/2stream_master.f90      \
   $(2S_MAIN_PATH)/2stream_Planckfunc.f90  \
   $(2S_TEST_PATH)/test_2s_ONLY_2p4.f90

SOURCES_LPCSMASTER =
SOURCES_LPCSMASTER +=   \
   $(SOURCES) $(SOURCES_SUP)                 \
   $(SOURCES_LS) $(SOURCES_LSSUP)            \
   $(SOURCES_LA) $(SOURCES_LP) $(SOURCES_LC) \
   $(2S_MAIN_PATH)/2stream_lcs_master.f90    \
   $(2S_MAIN_PATH)/2stream_lps_master.f90    \
   $(2S_MAIN_PATH)/2stream_Planckfunc.f90    \
   $(2S_TEST_PATH)/test_2s_LPCS_2p4.f90

# Parallel tests using OpenMP

SOURCES_MASTER_OMP  = 
SOURCES_MASTER_OMP +=   \
   $(SOURCES) $(SOURCES_SUP)               \
   $(2S_MAIN_PATH)/2stream_master.f90      \
   $(2S_MAIN_PATH)/2stream_Planckfunc.f90  \
   $(2S_TEST_PATH)/test_2s_ONLY_2p4_OMP.f90

SOURCES_LPCSMASTER_OMP =
SOURCES_LPCSMASTER_OMP +=   \
   $(SOURCES) $(SOURCES_SUP)                 \
   $(SOURCES_LS) $(SOURCES_LSSUP)            \
   $(SOURCES_LA) $(SOURCES_LP) $(SOURCES_LC) \
   $(2S_MAIN_PATH)/2stream_lcs_master.f90    \
   $(2S_MAIN_PATH)/2stream_lps_master.f90    \
   $(2S_MAIN_PATH)/2stream_Planckfunc.f90    \
   $(2S_TEST_PATH)/test_2s_LPCS_2p4_OMP.f90

# Utilities

SOURCES_UTIL =
SOURCES_UTIL +=   \
   $(UTIL_PATH)/2s_diff.f90


#
# Define pattern rules for creating object files:
#

.SUFFIXES:

# For 2s main source files
$(OBJ_PATH)/%.o : $(2S_MAIN_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(2S_TEST_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@

# For utility source files
$(OBJ_PATH)/%.o : $(UTIL_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@


#
# Define object files
#

F90SOURCES_M := $(notdir $(filter %.f90, $(SOURCES_MASTER)))
F90OBJECTS_M := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_M)))

F90SOURCES_MLPCS := $(notdir $(filter %.f90, $(SOURCES_LPCSMASTER)))
F90OBJECTS_MLPCS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_MLPCS)))

F90SOURCES_M_OMP := $(notdir $(filter %.f90, $(SOURCES_MASTER_OMP)))
F90OBJECTS_M_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_M_OMP)))

F90SOURCES_MLPCS_OMP := $(notdir $(filter %.f90, $(SOURCES_LPCSMASTER_OMP)))
F90OBJECTS_MLPCS_OMP := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_MLPCS_OMP)))

F90SOURCES_UTIL := $(notdir $(filter %.f90, $(SOURCES_UTIL)))
F90OBJECTS_UTIL := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_UTIL)))


#
# Define desired targets
#

main: test_2s_ONLY_2p4.exe  test_2s_LPCS_2p4.exe

openmp: test_2s_ONLY_2p4_OMP.exe  test_2s_LPCS_2p4_OMP.exe

test_2s_ONLY_2p4.exe: $(F90OBJECTS_M)
	$(FC) $^ -o $@

test_2s_LPCS_2p4.exe: $(F90OBJECTS_MLPCS)
	$(FC) $^ -o $@

test_2s_ONLY_2p4_OMP.exe: $(F90OBJECTS_M_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@

test_2s_LPCS_2p4_OMP.exe: $(F90OBJECTS_MLPCS_OMP)
	$(FC) $(FFLAGS_OPENMP) $^ -o $@

#test_2s_ONLY_2p4_OMP.exe: $(F90OBJECTS_M_OMP)
#	$(FC) $^ -o $@

#test_2s_LPCS_2p4_OMP.exe: $(F90OBJECTS_MLPCS_OMP)
#	$(FC) $^ -o $@

2s_diff: $(F90OBJECTS_UTIL)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe

