GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O2 -fPIC  -I./ -std=c++11
## to use RooUnfold

lxplus=$(findstring lxplus, $(shell hostname -f) )
mit=$(findstring mit.edu, $(shell hostname -f) )
dtmit=$(findstring dtmit, $(shell hostname -f) )

###################### determine if you are on lxplus or mit or else
ifneq ($(strip $(lxplus)),)
$(info You are on lxplus)
#ROOUNFOLD=$(HOME)/work/RooUnfold/
ROOUNFOLD=/afs/cern.ch/user/a/amarini/public/RooUnfold/
else 
ifneq ($(strip $(dtmit)),)
$(info You are on a dtmit, using lxplus)
ROOUNFOLD=/afs/cern.ch/user/a/amarini/public/RooUnfold/
else 
ifneq ($(strip $(mit)),)
		## ON MIT
$(info You are on a mit server)
ROOUNFOLD=/home/amarini/RooUnfold/
else
$(info You are Not on lxplus nor mit)
ROOUNFOLD=$(HOME)/Downloads/RooUnfold-1.1.1/
endif
endif
endif

#CXXFLAGS += -L$(PWD)/../NeroProducer/Core/bin -lBare -Wl,-rpath=$(PWD)/../NeroProducer/Core/bin -ggdb
CXXFLAGS += -L$(ROOUNFOLD)/ -lRooUnfold -ggdb
CXXFLAGS += -I$(ROOUNFOLD)/src/
DICTFLAGS+= -I$(ROOUNFOLD)/src/
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin
HPPDIR=interface

SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJ=$(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o , $(SRC)  )
HPPLINKDEF=$(patsubst $(SRCDIR)/%.cpp, ../interface/%.hpp , $(SRC)  )

.PHONY: all
all: libBootStrap.so
	$(info, "--- Full compilation --- ")	
	$(MAKE) libBootStrap.so


# check if CMSSW is defined
ifndef CMSSW_BASE
$(info No CMSSSW !!!!)
$(info I ll sleep 3s to let you acknowledge it.)
$(info To avoid sleeping 'touch NO_CMSSW')
$(shell [ -f NO_CMSSW ] || sleep 3s)
else
$(info CMSSW found: $(CMSSW_BASE) )
endif

## if on mac add the -p to the  DICTFLAGS
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
$(info You are compiling on mac)
DICTFLAGS += -p
else 
$(info Your are on a linux machine)
CXXFLAGS +=  -Wl,-rpath=$(ROOUNFOLD) 
endif

# check if Combine is present and compiled 

libBootStrap.so: $(OBJ) Dict | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(SOFLAGS) -o $(BINDIR)/$@ $(OBJ) $(BINDIR)/dict.o

$(OBJ) : $(BINDIR)/%.o : $(SRCDIR)/%.cpp interface/%.hpp | $(BINDIR)
	$(GCC) $(CXXFLAGS) -c -o $(BINDIR)/$*.o $<

.PHONY: Dict
Dict: $(BINDIR)/dict.o

$(BINDIR)/dict.o: $(SRC) | $(BINDIR)
	cd $(BINDIR) && rootcint -v4 -f dict.cc -c -I./ -I../ $(DICTFLAGS)  $(HPPLINKDEF)  ../interface/LinkDef.hpp 
	cd $(BINDIR) && $(GCC) -c -o dict.o $(CXXFLAGS) -I../../ -I../ dict.cc

$(BINDIR):
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm $(BINDIR)/dict*
	-rm $(BINDIR)/libBootStrap.so
	-rmdir $(BINDIR)
