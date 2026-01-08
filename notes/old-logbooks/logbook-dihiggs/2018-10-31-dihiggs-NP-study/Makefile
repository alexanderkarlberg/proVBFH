
# Makefile generated automatically by ../../../Software/generator-framework/mkmk.pl ''
# run 'make make' to update it if you add new files

# get material from common/ directory
LCLINCLUDE   += -I../../../Software/generator-framework
LIBRARIES += -L../../../Software/generator-framework -lhepmcplus

include ../../../Software/generator-framework/Makefile.inc

LIBRARIES += -lhepmcplus
INCLUDE += $(LCLINCLUDE)
INCLUDE += -DCOMMONDIR='"/home/frederic/Work/gitlab_dihiggs/logbook/2018-10-31-dihiggs-NP-study/../../../Software/generator-framework/"'

COMMONSRC = 
F77SRC = 
COMMONOBJ = $(patsubst %.cc,%.o,$(COMMONSRC)) $(patsubst %.f,%.o,$(F77SRC))

PROGSRC = vbfhh-study.cc
PROGOBJ = $(patsubst %.cc,%.o,$(PROGSRC))


vbfhh-study: vbfhh-study.o ../../../Software/generator-framework/libhepmcplus.a $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES) 


make:
	../../../Software/generator-framework/mkmk.pl ''

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -c $< -o $@
.f.o:         $<
	$(F77) $(FFLAGS) -c $< -o $@


depend:
	makedepend   $(LCLINCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE

vbfhh-study.o: ../../../Software/generator-framework/AnalysisFramework.hh
vbfhh-study.o: ../../../Software/generator-framework/SimpleNTuple.hh
vbfhh-study.o: ../../../Software/generator-framework/SimpleHist.hh
vbfhh-study.o: ../../../Software/generator-framework/SimpleHist2D.hh
vbfhh-study.o: ../../../Software/generator-framework/AveragingHist.hh
vbfhh-study.o: ../../../Software/generator-framework/CorrelationHist.hh
vbfhh-study.o: ../../../Software/generator-framework/AverageAndError.hh
vbfhh-study.o: ../../../Software/generator-framework/NLOHistGeneric.hh
vbfhh-study.o: ../../../Software/generator-framework/CleverStream.hh
vbfhh-study.o: ../../../Software/generator-framework/CommonAnalysis.hh
vbfhh-study.o: ../../../Software/generator-framework/ExtractedEvent.hh
vbfhh-study.o: ../../../Software/generator-framework/addons.hh
vbfhh-study.o: ../../../Software/generator-framework/FlavourHolder.hh
vbfhh-study.o: ../../../Software/generator-framework/FJHepMC.hh
vbfhh-study.o: ../../../Software/generator-framework/FlexiBTagger.hh
vbfhh-study.o: ../../../Software/generator-framework/Calorimeter.hh
vbfhh-study.o: ../../../Software/generator-framework/DriverFromCmdLine.hh
vbfhh-study.o: ../../../Software/generator-framework/AnyDriver.hh
vbfhh-study.o: ../../../Software/generator-framework/LockableMap.hh
