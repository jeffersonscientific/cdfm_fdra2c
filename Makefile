# Usage:
#   make mode=p opt=-O
#     parallelization by MPI
#   make mode=s opt=-O
#     serial with no MPI dependencies
# Change these lines to fit your system.
#
# TODO: either set up a locla directory structure or use modules... I think HMVP and the UTILs
# are multi-purpose, so a generalized installation would probably be good.
HMMVPINC = -I $(HMMVP_DIR)/include -I external
HMMVPLIB = -L $(HMMVP_DIR)/lib/ -lhmmvp_$(mode)
UTILINC = -I $(UTIL_DIR)/include
UTILLIB = -L $(UTIL_DIR)/lib -lutil_$(mode)
BLASLIB = -lblas

CPP = $(CXX)
MPICPP = $(MPICXX)
# The rest should not have to be changed.

CPPSRCS = src/main.cpp src/OdeFn.cpp src/Ode.cpp src/StreamDataFile.cpp	\
src/Fdra.cpp src/HmatrixStressFn.cpp #\
/home/camcat/code/ambrad/util/src/ValueSetter.cpp \
/home/camcat/code/ambrad/util/src/KeyValueFile.cpp \
/home/camcat/code/ambrad/util/src/Mpi.cpp \
/home/camcat/code/ambrad/util/src/CodeAnalysis.cpp \

INCLUDE = -I src -I include $(HMMVPINC) $(UTILINC)
LIBS = $(HMMVPLIB) $(UTILLIB) $(BLASLIB) -lgfortran -lgomp
LIBDIRS =
OPTFLAGS = $(opt)
CPPFLAGS = $(OPTFLAGS) #-Wall
LDFLAGS = $(OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJECTS = $(patsubst %.cpp,%.o,$(CPPSRCS))

%.o: %.cpp
ifeq ($(mode),s)
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
else
	$(MPICPP) -DUTIL_MPI $(CPPFLAGS) $(INCLUDE) -c $< -o $@
endif

all: $(OBJECTS)
ifeq ($(mode),s)
	$(CPP) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra
else
	$(MPICPP) $(LDFLAGS) $(LIBFLAGS) -DUTIL_MPI $(LIBDIRS) $(OBJECTS) $(LIBS) -o bin/fdra
endif

add_bc:
	$(CPP) test/add_bc_to_fdra_kvf.cpp $(UTILLIB) $(INCLUDE) -o bin/add_bc_to_fdra_kvf

fdraclang:
	clang -fsyntax-only $(INCLUDE) $(CPPSRCS)

clean:
	rm -f src/*.o hmmvp/*.o bin/fdra
