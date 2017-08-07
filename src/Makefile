# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/fvtest.app
EC_scheme.o: ./EC_scheme.f90 decimal.o FV_types.o scheme_utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./EC_scheme.f90
FV_Solve.o: ./FV_Solve.f90 decimal.o FV_types.o scheme_utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_Solve.f90
FV_types.o: ./FV_types.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_types.f90
KT_scheme.o: ./KT_scheme.f90 decimal.o FV_types.o scheme_utils.o limiters.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./KT_scheme.f90
SKT_scheme.o: ./SKT_scheme.f90 decimal.o FV_types.o scheme_utils.o limiters.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./SKT_scheme.f90
CU_scheme.o: ./CU_scheme.f90 decimal.o FV_types.o scheme_utils.o limiters.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./CU_scheme.f90
DRCU_scheme.o: ./DRCU_scheme.f90 decimal.o FV_types.o scheme_utils.o limiters.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./DRCU_scheme.f90
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
limiters.o: ./limiters.f90 decimal.o FV_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./limiters.f90
main.o: ./main.f90 test1.o test2.o test3.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
plot.o: ./plot.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
scheme_utils.o: ./scheme_utils.f90 decimal.o FV_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./scheme_utils.f90
test1.o: ./test1.f90 decimal.o plot.o FV_types.o SKT_scheme.o EC_scheme.o FV_Solve.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test1.f90
test2.o: ./test2.f90 decimal.o plot.o FV_types.o KT_scheme.o EC_scheme.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test2.f90
test3.o: ./test3.f90 decimal.o plot.o FV_types.o KT_scheme.o EC_scheme.o FV_Solve.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test3.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./decimal.f90 ./test1.f90 ./main.f90 ./util.f90 ./KT_scheme.f90 ./SKT_scheme.f90 ./CU_scheme.f90 ./DRCU_scheme.f90 ./test3.f90 ./EC_scheme.f90 ./scheme_utils.f90 ./FV_types.f90 ./limiters.f90 ./test2.f90 ./FV_Solve.f90 ./plot.f90
OBJ = decimal.o test1.o main.o util.o KT_scheme.o SKT_scheme.o CU_scheme.o DRCU_scheme.o test3.o EC_scheme.o scheme_utils.o FV_types.o limiters.o test2.o FV_Solve.o plot.o
clean: neat
	-rm -f .../bin/fvtest.app.cppdefs $(OBJ) ../bin/fvtest.app *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
../bin/fvtest.app: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/fvtest.app  $(LDFLAGS)
