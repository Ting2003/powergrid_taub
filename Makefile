CC=mpicxx
SRC= util.cpp point.cpp node.cpp circuit.cpp net.cpp parser.cpp vec.cpp \
    main.cpp triplet.cpp algebra.cpp block.cpp mpi_class.cpp transient.cpp sp_graph_table.cpp sp_node.cpp  
#hash_mat.cpp map_mat.cpp 
HDR=$(SRC:.cpp=.h)
OBJ=$(SRC:.cpp=.o) 
BIN=pg
RELEASE=IPGS
CPPFLAGS=
CFLAGS=-lg2c -lrt -Wall -Wextra -pipe -O2 -mfpmath=sse #-msse4.2 -mssse3 -mfpmath=sse -march=native
#CFLAGS=-Wall -Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=native
#CFLAGS=-Wall -g #-Wextra -pipe -O2 -msse4.2 -mssse3 -mfpmath=sse -march=core2
#LDFLAGS=-s -Wl,-O1,-hash-style=gnu
LDFLAGS=

# PACKAGE= /home/aaslyan/Suite/SuiteSparse
PACKAGE= ../../package/SuiteSparse
#./package_ck

UMFPACK=#./umfpack
UMFPACK_LIB_DIR=#$(UMFPACK)/lib
UMFPACK_INC_DIR=#$(UMFPACK)/include
UMFPACK_LIB=#$(UMFPACK_LIB_DIR)/libumfpack.a \
	    $(UMFPACK_LIB_DIR)/libamd.a \
	    $(CHOLMOD_LIB_DIR)/libcholmod.a \
	    $(UMFPACK_LIB_DIR)/libcolamd.a \
            $(UMFPACK_LIB_DIR)/libccolamd.a \
            $(UMFPACK_LIB_DIR)/libcamd.a \
            $(UMFPACK_LIB_DIR)/libmetis.a \
            $(UMFPACK_LIB_DIR)/libgoto2.a

CHOLMOD= $(PACKAGE)/CHOLMOD
CHOLMOD_LIB_DIR=$(CHOLMOD)/Lib
CHOLMOD_INC_DIR=$(CHOLMOD)/Include
CHOLMOD_LIB=$(CHOLMOD_LIB_DIR)/libcholmod.a\
	    $(PACKAGE)/AMD/Lib/libamd.a\
	    $(PACKAGE)/COLAMD/Lib/libcolamd.a\
	    $(PACKAGE)/CCOLAMD/Lib/libccolamd.a\
	    $(PACKAGE)/CAMD/Lib/libcamd.a \
            ../../package/metis-5.1.0/libmetis/libmetis.a \
	    ../../package/lapack-3.1.0/lapack_LINUX.a \
	    ../../package/lapack-3.1.0/blas_LINUX.a\
	    ../../package/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a\
	    /usr/lib/gcc/x86_64-redhat-linux/3.4.6/libg2c.a \
	    /usr/lib64/librt.so

main: $(OBJ)
	@echo "Making project..."
	$(CC) $(LDFLAGS) -o $(BIN) $(OBJ) $(CHOLMOD_LIB)

release: $(OBJ)
	$(CC) $(LDFLAGS) -static -o $(BIN) $(OBJ) $(CHOLMOD_LIB)

test: 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS)\
	-I$(CHOLMOD_INC_DIR) -o test test.cpp $(CHOLMOD_LIB)

all: main
	@echo "Making all..."

%.o: %.cpp  %.h global.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(CHOLMOD_INC_DIR) -c $<  -o $@

.PHONY : clean
clean:
	@echo "Cleaning all..."
	rm -rf *.o $(OBJ) $(DBG) $(BIN) #tags $(CSCOPEFILES)
	rm ./INPUT_FILE/*
