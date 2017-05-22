
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common
PDAFLIB = $(CLAW)/geoclaw/src/2d/shallow/assimilation

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
PDAF_BASEDIR_INC=$(CLAW)/geoclaw/examples/tsunami/PDAF-D_V1.12/src/

#FFLAGS ?= -cpp -DCHILE_GEN_ENS
#FFLAGS ?= -cpp -DPKJ_DEBUG -DUSE_PDAF -DUSE_PDAF_CHILE -DUSE_PDAF_CHILE_PERT_RESTART -I$(PDAF_BASEDIR_INC) -L$(PDAFLIB) -I$(PDAFLIB) -lpdaf-d -lblas -llapack 
FFLAGS ?= -cpp
#FFLAGS ?= -cpp -DCHILE_GEN_ENS -g -Wall -pedantic -fbounds-check -ffpe-trap=invalid,overflow,zero

FC=mpif90

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
  $(GEOLIB)/pkj_geoclaw/mod_parallel_ens_gen.f90
  #$(GEOLIB)/pkj_geoclaw/checker_same_column.f90 \
  #$(GEOLIB)/pkj_geoclaw/common_level.f90

MODULES2 = \
  $(GEOLIB)/pkj_geoclaw/checker_same_column.f90 \
  $(GEOLIB)/pkj_geoclaw/common_level.f90

SOURCES = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \
  
SOURCES2 = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \
  $(GEOLIB)/pkj_geoclaw/normal.f90 \
  $(GEOLIB)/pkj_geoclaw/bufnst2_assim.f \
  $(GEOLIB)/pkj_geoclaw/alloc2field.f90 \
  $(GEOLIB)/pkj_geoclaw/field2alloc.f90 \
  $(GEOLIB)/pkj_geoclaw/get_dim_state.f90 \
  $(GEOLIB)/pkj_geoclaw/get_obs.f90 \
  $(GEOLIB)/pkj_geoclaw/print_num_cells.f90 \
  $(GEOLIB)/pkj_geoclaw/update_dim_state_p.f90 \
  $(GEOLIB)/pkj_geoclaw/flag2refine2_assim.f90 \
  $(GEOLIB)/pkj_geoclaw/flagregions2_assim.f90 \
  $(GEOLIB)/pkj_geoclaw/checker_ensemble_grid.f90 \
  $(GEOLIB)/pkj_geoclaw/put2zero.f90 \
  $(GEOLIB)/pkj_geoclaw/set_global_coordinate_array.f90 \
  $(GEOLIB)/pkj_geoclaw/extract_regions.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

RESTART = False
# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls
