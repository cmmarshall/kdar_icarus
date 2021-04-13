
CXX = g++

CPPFLAGS = -g -fno-inline -O0

ROOTFLAGS = `root-config --cflags --glibs`
ROOTINC = `root-config --cflags`
ROOTLIB = `root-config --glibs`

INCLUDE += -I$(GENIE)/src -I$(LOG4CPP_INC)

LDLIBS += -L$(LOG4CPP_LIB) -llog4cpp
LDLIBS += -L/usr/lib64 -lxml2
LDLIBS += -L$(PYTHIA6) -lPythia6
LDLIBS += -L$(GSL_LIB) -lgsl -lgslcblas
# ROOT Flags are incomplete.
LDLIBS += -L$(ROOTSYS)/lib -lGeom -lEGPythia6
# GENIE libraries
LDLIBS += -L$(GENIE_LIB) \
                              -lGFwAlg \
                              -lGFwEG \
                              -lGFwGHEP \
                              -lGFwInt \
                              -lGFwMsg \
                              -lGFwNtp \
                              -lGFwNum \
                              -lGFwParDat \
                              -lGFwReg \
                              -lGFwUtl \
                              -lGPhAMNGEG \
                              -lGPhAMNGXS \
                              -lGPhBDMEG \
                              -lGPhBDMXS \
                              -lGPhChmXS \
                              -lGPhCmn \
                              -lGPhCohEG \
                              -lGPhCohXS \
                              -lGPhDcy \
                              -lGPhDeEx \
                              -lGPhDfrcEG \
                              -lGPhDfrcXS \
                              -lGPhDISEG \
                              -lGPhDISXS \
                              -lGPhGlwResEG \
                              -lGPhGlwResXS \
                              -lGPhHadnz \
                              -lGPhHadTransp \
                              -lGPhIBDEG \
                              -lGPhIBDXS \
                              -lGPhMEL \
                              -lGPhMNucEG \
                              -lGPhMNucXS \
                              -lGPhNDcy \
                              -lGPhNNBarOsc \
                              -lGPhNuclSt \
                              -lGPhNuElEG \
                              -lGPhNuElXS \
                              -lGPhPDF \
                              -lGPhQELEG \
                              -lGPhQELXS \
                              -lGPhResEG \
                              -lGPhResXS \
                              -lGPhStrEG \
                              -lGPhStrXS \
                              -lGPhXSIg \
                              -lGRwClc \
                              -lGRwFwk \
                              -lGRwIO \
                              -lGTlFlx \
                              -lGTlGeo 

all: flatten

flatten: flatten.o 
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o flatten flatten.cxx $(LDLIBS) 
flatten.o: flatten.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c flatten.cxx -o flatten.o 

clean:
	-rm -f flatten
	-rm -f flatten.o

remake:
	make clean
	make all
