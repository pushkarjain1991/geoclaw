# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################


# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR = /h2/pkjain/Desktop/Pushkar/clawpack/geoclaw/examples/tsunami/PDAF-D_V1.11.1
# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
#PDAF_ARCH = linux_gfortran
#PDAF_ARCH=linux_gfortran
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h


# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  mod_model.o
#USER_SUPPLIED_MODULES= alloc_field_exchange.o\
                         mapdomain.o 




# Model routines
OBJ_MODEL =	main.o\
		initialize.o \
		integrate.o

# Moduls used for PDAF
MOD_ASSIM = 	mod_parallel.o \
		parser_mpi.o \
		mod_assimilation.o 

# Model routines used with PDAF
OBJ_MODEL_PDAF =main_pdaf.o\
		initialize.o \
		integrate_pdaf.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		finalize_pdaf.o assimilate_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens.o \
		prepoststep_ens_pdaf.o \
		next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o \
		init_dim_obs_pdaf.o \
		obs_op_pdaf.o \
		init_obs_pdaf.o \

# User-supplied routines for ESTKF/SEIK/ETKF
OBJ_USER_SEIK = prodrinva_pdaf.o \
		init_obsvar_pdaf.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = add_obs_error_pdaf.o \
		init_obscovar_pdaf.o \

# User-supplied routines for localized analysis (LESTKF/LSEIK/LETKF)
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o \
		prodrinva_l_pdaf.o \
		init_obs_l_pdaf.o \
		init_dim_obs_l_pdaf.o \
		g2l_obs_pdaf.o \
		init_dim_obs_f_pdaf.o \
		init_obs_f_pdaf.o \
		obs_op_f_pdaf.o   \
		init_obsvar_l_pdaf.o \

#OBJ_USER_DOMAIN = oned_to_twod.o \
#		  twod_to_oned.o

# Full list of user-supplied routines for online modes
#OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)
#OBJ_PDAF_USER = $(OBJ_USER_DOMAIN) $(OBJ_USER_GEN) $(OBJ_USER_SEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)
OBJ_PDAF_USER = $(OBJ_USER_DOMAIN) $(OBJ_USER_GEN) $(OBJ_USER_SEIK) $(OBJ_USER_ENKF)

######################################################

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D serial model (without parallelization)";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - just the forward model";
	@echo "  model_pdaf    - assimilation model (use with -DUSE_PDAF!)";
	@echo "  (if compiled without -DUSE_PDAF, both programs behave equally)"
	@echo "  libpdaf-a     - PDAF library"
	@echo "  libonline     - PDAF library with PDAF modified geoclaw src code "
	@echo "  (to be compiled with -DUSE_PDAF)"
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo "List available settings for PDAF_ARCH with 'make listarch'";
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";

model : $(MODULES) $(OBJ_MODEL)
	$(LD)  $(OPT_LNK)  -o $@  $(MODULES) $(OBJ_MODEL)

model_pdaf : libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) \
	-L$(BASEDIR)/lib  -lpdaf-d $(LINK_LIBS)

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/src; make;

#libonline: libpdaf-d.a $(MODULES) $(OBJ_MPI) $(MOD_ASSIM) $(OBJ_PDAF_INT) $(OBJ_USER_GEN) $(OBJ_USER_ENKF) $(OBJ_USER_SEIK)
#libonline: libpdaf-d.a $(MODULES) $(OBJ_USER_DOMAIN) $(OBJ_MPI) $(MOD_ASSIM) $(OBJ_PDAF_INT) $(OBJ_USER_GEN) $(OBJ_USER_ENKF) $(OBJ_USER_SEIK) $(OBJ_USER_LOCAL)
libonline: libpdaf-d.a $(MODULES) $(OBJ_USER_DOMAIN) $(OBJ_MPI) $(MOD_ASSIM) $(OBJ_PDAF_INT) $(OBJ_USER_GEN) $(OBJ_USER_ENKF) $(OBJ_USER_SEIK)
	@echo "++++++ Generate assimilation library ++++++"
	@cp $(BASEDIR)/lib/libpdaf-d.a .
	@ar r $^
	


######################################################

%.o: %.F90
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

cleanall : cleanpdaf clean cleandata

clean :
	rm -f *.o *.mod *.a model model_pdaf

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR)/src; make clean

cleandata:
	rm -f ens*ana.txt ens*for.txt state_*ana.txt state_*for.txt

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1
debug:
#	@echo $(PDAF_ARCH)
#	@echo $(BASEDIR)/make.arch/$(PDAF_ARCH).h
	@echo $(FC)
