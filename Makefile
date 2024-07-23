
FC       = mpif90
FFLAGS   = -O2
#FFLAGS   = -O2 -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow
LDFLAGS  =
CPP       = -cpp

# metis library
METIS_DIR  = ./submodule/monolis
METIS_INC  =
METIS_LIB  = -L$(METIS_DIR)/lib -lmetis

# mumps library
MUMPS_DIR  = ./submodule/monolis
MUMPS_INC  = -I $(MUMPS_DIR)/include
MUMPS_LIB  = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord

# for ubuntu env
filename = /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so
fileexists = $(shell ls ${filename} | grep ${filename})
ifeq (${fileexists}, ${filename})
  SCALAPACK = -lscalapack-openmpi
else
  SCALAPACK = -lscalapack
endif

BLAS_LIB   = ${SCALAPACK} -llapack -lblas

# monolis library
MONOLIS_DIR= ./submodule/monolis
MONOLIS_INC= -I $(MONOLIS_DIR)/include
MONOLIS_LIB= -L$(MONOLIS_DIR)/lib -lmonolis

LIBS     = $(MONOLIS_LIB) $(MUMPS_LIB) $(METIS_LIB) $(BLAS_LIB)
INCLUDE  = -I ./include $(MONOLIS_INC)
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./obj
BIN_LIST = monolis_eigen
TARGET   = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
SRC_LIST = util.f90 debug.f90 io.f90 element_C3D8.f90 matrix.f90 update.f90 solver.f90 analysis.f90 main.f90
SOURCES  = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
OBJS     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
RM       = rm

$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(CPP) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET) ./include/*.mod

.PHONY: clean
