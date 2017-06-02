# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/test.app
EC_scheme.o: ./EC_scheme.f90 decimal.o FV_types.o scheme_utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./EC_scheme.f90
FV_Solve.o: ./FV_Solve.f90 decimal.o FV_types.o scheme_utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_Solve.f90
FV_types.o: ./FV_types.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./FV_types.f90
KT_scheme.o: ./KT_scheme.f90 decimal.o FV_types.o scheme_utils.o limiters.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./KT_scheme.f90
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
test1.o: ./test1.f90 decimal.o plot.o FV_types.o KT_scheme.o EC_scheme.o FV_Solve.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test1.f90
test2.o: ./test2.f90 decimal.o plot.o FV_types.o KT_scheme.o EC_scheme.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test2.f90
test3.o: ./test3.f90 decimal.o plot.o FV_types.o KT_scheme.o EC_scheme.o FV_Solve.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./test3.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./FV_Solve.f90 ./KT_scheme.f90 ./plot.f90 ./test1.f90 ./scheme_utils.f90 ./test2.f90 ./EC_scheme.f90 ./decimal.f90 ./test3.f90 ./FV_types.f90 ./limiters.f90 ./util.f90 ./main.f90
OBJ = FV_Solve.o KT_scheme.o plot.o test1.o scheme_utils.o test2.o EC_scheme.o decimal.o test3.o FV_types.o limiters.o util.o main.o
clean: neat
	-rm -f .../bin/test.app.cppdefs $(OBJ) ../bin/test.app *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
../bin/test.app: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/test.app  $(LDFLAGS)
