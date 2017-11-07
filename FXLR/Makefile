FC = gfortran
OBJ = dvode_f90_m.o fxlr.o
OBJ2 = dvode_f90_m.o fxlr1d.o
OPT = -O3 -fopenmp

fxlr: $(OBJ)
	$(FC) $(OPT) $^ -o $@
%.o: %.f90
	$(FC) $(OPT) -c $^
fxlr1d: $(OBJ2)
	$(FC) $(OPT) $^ -o $@
%.o: %.f90
	$(FC) $(OPT) -c $^
.PHONY: clean
clean:
	-rm $(OBJ) *.mod *.a
