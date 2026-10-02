#############################################################################
# proVBFH-cs: VBF H at (N)NLO, fully differential, with a line-by-line
# projection-to-Born (see docs/DESIGN.md). No POWHEG-BOX machinery: the
# inclusive part (structure functions), the analysis and the VBFNLO matrix
# elements are taken from ../proVBFH; only a few POWHEG headers and
# utilities (the input-file reader, newunit, spinor routines) are needed by
# them.
#
# "make" builds proVBFH-cs; "make clean" removes object files.
# "make ANALYSIS=hxswg136" builds proVBFH-cs-hxswg136 with the analysis
# analysis/hxswg136_analysis.f instead of proVBFH's user_analysis.f.
#############################################################################

# external tools (hoppet, LHAPDF, FastJet), as configured for proVBFH
include ../proVBFH/Makefile.inc

LHAPDFLIBS=-Wl,-rpath,$(shell $(LHEXEC) --prefix)/lib -L$(shell $(LHEXEC) --prefix)/lib -lLHAPDF
FASTJETLIBS=$(shell $(FJEXEC) --libs --plugins)
HOPPETLIBS=-Wl,-rpath,$(shell $(HPEXEC) --prefix)/lib $(shell $(HPEXEC) --libs)
LIBS=$(LHAPDFLIBS) $(HOPPETLIBS) $(FASTJETLIBS)

FC = gfortran
FFLAGS = -O3 -fno-automatic -ffixed-line-length-132 -std=legacy
F90FLAGS = -O3 -fno-automatic -ffree-line-length-none $(shell $(HPEXEC) --fflags)
CXX ?= g++
CXXFLAGS = -O3 $(shell $(FJEXEC) --cxxflags)
CXXLIBS = -lstdc++

PV=../proVBFH
SRC=src
OBJ=obj
PWG=$(PV)/src/powheg-files
HJJJ=$(PWG)/hjjj-files

INCLUDE= -I$(SRC) -I$(PV)/src -I$(PV)/src/inclusive -I$(PV)/src/exclusive \
	-I$(PWG) -I$(PWG)/include -I$(PV)/analysis -I$(HJJJ) -I$(HJJJ)/vbfnlo-files
MODDIR= -J$(OBJ)

VPATH= ./:$(SRC)/:analysis/:$(PV)/src/inclusive/:$(PV)/src/exclusive/:$(PV)/src/exclusive/vbfnlo-files/
VPATH+=:$(PV)/analysis/:$(PWG)/:$(HJJJ)/:$(HJJJ)/vbfnlo-files/:$(OBJ)/

# inclusive part (structure functions, VBF Born phase space), from proVBFH
INCL= incl_vbfh.o phase_space.o matrix_element.o incl_parameters.o \
	integration.o io_utils.o lcl_dec.o tensor.o nonfact.o \
	nonfact_expressions.o helper.o constants.o

# VBFNLO matrix elements (option A) and their coupling set-up, from proVBFH
VBFNLO= brakets_new.o coeff.o dotproducts.o flavorsum.o \
	hjjj_amp_aux_corrected.o dilog.o HtoWW.o kinematics_aux.o loops_aux.o \
	parameters.o qqhqqj-born.o qqhqqj-virt.o qqh2q2g_me_new.o qqh4q_me_new.o \
	Born.o convert_coup.o init_couplings.o calc_als.o vbfnloinput.o \
	real_vbfnlo.o compreal_hqqqq_new.o init_processes.o virtual.o

# the few POWHEG utilities these need (input reader, unit numbers, spinor
# and boost routines, tag encoding; no integration machinery)
PWGUTIL= powheginput.o newunit.o bra_ket_subroutines.o boostrot.o cs_powheg_tags.o

# analysis, from proVBFH, or analysis/<ANALYSIS>_analysis.f
ifdef ANALYSIS
ANAL= $(ANALYSIS)_analysis.o pwhg_bookhist-multi.o fastjetfortran.o
EXE= proVBFH-cs-$(ANALYSIS)
else
ANAL= user_analysis.o pwhg_bookhist-multi.o fastjetfortran.o
EXE= proVBFH-cs
endif

# this package
CS= cs_glue.o cs_glue2.o cs_born2.o cs_kinematics.o cs_dipoles.o cs_memo.o cs_nlo2.o cs_exclusive.o cs_main.o

OBJS= $(INCL) $(VBFNLO) $(PWGUTIL) $(ANAL) $(CS)

all: $(EXE)

# unit tests of the kinematics (no matrix elements)
tests/test_kinematics: cs_kinematics.o tests/test_kinematics.f90
	$(FC) $(F90FLAGS) -I$(OBJ) tests/test_kinematics.f90 $(OBJ)/cs_kinematics.o -o $@
tests/test_four: cs_kinematics.o cs_dipoles.o tests/test_four.f90
	$(FC) $(F90FLAGS) -I$(OBJ) tests/test_four.f90 $(OBJ)/cs_kinematics.o $(OBJ)/cs_dipoles.o -o $@

# development tool: proVBFH's real flavour list and tags
tests/list_reals: $(OBJ) $(filter-out cs_main.o,$(OBJS)) tests/list_reals.f
	$(FC) $(FFLAGS) $(INCLUDE) tests/list_reals.f $(patsubst %,$(OBJ)/%,$(filter-out cs_main.o,$(OBJS))) $(LIBS) $(CXXLIBS) -o $@

$(OBJ):
	mkdir -p $(OBJ)

$(EXE): $(OBJ) $(OBJS)
	$(FC) $(F90FLAGS) $(patsubst %,$(OBJ)/%,$(OBJS)) $(LIBS) $(CXXLIBS) -o $@

%.o: %.f | $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) $(MODDIR) -c -o $(OBJ)/$@ $<
%.o: %.F | $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) $(MODDIR) -c -o $(OBJ)/$@ $<
%.o: %.f90 | $(OBJ)
	$(FC) $(F90FLAGS) $(INCLUDE) $(MODDIR) -c -o $(OBJ)/$@ $<
%.o: %.cc | $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $(OBJ)/$@ $<

# f90 module dependencies
incl_vbfh.o: matrix_element.o incl_parameters.o phase_space.o integration.o io_utils.o lcl_dec.o
incl_parameters.o: integration.o
matrix_element.o: incl_parameters.o tensor.o nonfact.o
nonfact.o: helper.o constants.o incl_parameters.o nonfact_expressions.o
constants.o: helper.o
phase_space.o: matrix_element.o incl_parameters.o
cs_kinematics.o: incl_parameters.o
cs_dipoles.o: cs_kinematics.o
cs_nlo2.o: cs_kinematics.o cs_dipoles.o cs_memo.o
cs_exclusive.o: cs_kinematics.o cs_nlo2.o cs_memo.o incl_parameters.o phase_space.o matrix_element.o
cs_main.o: cs_exclusive.o cs_nlo2.o incl_vbfh.o incl_parameters.o phase_space.o integration.o cs_memo.o

clean:
	rm -rf $(OBJ)

.PHONY: all clean
