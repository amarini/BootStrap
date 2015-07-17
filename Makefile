GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O2 -fPIC  -I./
## to use RooUnfold
#CXXFLAGS += -L$(PWD)/../NeroProducer/Core/bin -lBare -Wl,-rpath=$(PWD)/../NeroProducer/Core/bin -ggdb
CXXFLAGS += -L/afs/cern.ch/user/a/amarini/work/RooUnfold/ -lRooUnfold -Wl,-rpath=/afs/cern.ch/user/a/amarini/work/RooUnfold/ -ggdb
CXXFLAGS += -I/afs/cern.ch/user/a/amarini/work/RooUnfold/src/
DICTFLAGS+= -I/afs/cern.ch/user/a/amarini/work/RooUnfold/src/
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin
HPPDIR=interface
AUXDIR=aux

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
$(info I ll sleep 3s to let you acknowledge it)
$(shell sleep 3s)
else
$(info CMSSW found: $(CMSSW_BASE) )
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
	mkdir -p $(AUXDIR)

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm $(BINDIR)/dict*
	-rm $(BINDIR)/libBootStrap.so
	-rmdir $(BINDIR)
