FC = gfortran
#FC = ifort
OBJ = dvode_f90_m.o gardner.o
OPT = -O3 -fopenmp

gardner: $(OBJ) $(OBJ1)
	$(FC) $(OPT) $^ -o $@
%.o: %.f90
	$(FC) $(OPT) -c $^
.PHONY: clean
clean:
	-rm $(OBJ) *.mod *.a
