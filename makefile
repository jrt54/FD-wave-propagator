CC		 = gcc 
CFLAGS=-fopenmp -I.
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = Jeremy
#CLINKER 	 = gcc
OBJ=
DEPS= fd.h utils.c



#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules

naivefd: naivefd.o $(DEPS) 
	 $(CC) -o $@ $^ $(CFLAGS) -lm 
	 rm -f naivefd.o

#utils.o

#naivefd: naivefd.o
#	-${CLINKER}  -o $@ $<
#	-${DSYMUTIL} $@

