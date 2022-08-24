
COMPILER=gfortran
FC=$(COMPILER)
ifeq ("$(FC)","gfortran")
F77 = gfortran
FFLAGS = -fno-automatic -ffixed-line-length-none -O2
endif
ifeq ("$(FC)","g77")
F77= g77 -fno-automatic -ffixed-line-length-none -O2
endif
ifeq ("$(FC)","ifort")
F77= ifort -save -132 -O3
endif


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
# object files for VBF-induced pp->Hjjj(j) to be implemented in POWHEG;
OBJECTS_PH =./dotproducts.o \
          ./parameters.o \
          ./HtoWW.o \
          ./flavorsum.o \
          ./kinematics_aux.o \
          ./brakets.o \
          ./qqhqqj-born.o  \
          ./qqhqqj-virt.o  \
          ./coeff.o \
          ./hjjj_amp_aux.o \
          ./loops_aux.o \
          ./qqhqqjj.o  \
          ./qqhqqQQ.o  \
          ./mg_qqqqggh_nc.o  \
          ./mg_qqqqggh_cc.o  \
          ./mg_qqqqqqh_nc.o  \
          ./mg_qqqqqqh_cc.o 
TARGETS_PH = libvbfnlo.a

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# rules for vbfnlo_lib:
#

 libvbfnlo.a:  $(OBJECTS_PH)
	$(AR) r libvbfnlo.a $(OBJECTS_PH) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.o.o:
	$(F77)  -c -o $*.o  $(FFLAGS) $*.o

# ----------------------------------------------------------------------------

clean:
	rm -f *.o */*.o */*/*.o */*/*/*.o core *~ *.x libvbfnlo.a 

