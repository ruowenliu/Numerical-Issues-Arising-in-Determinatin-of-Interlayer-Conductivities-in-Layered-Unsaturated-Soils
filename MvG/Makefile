FC = gfortran
OBJ = dvode_f90_m.o mvg.o
OPT = -O3 -fopenmp

mvg: $(OBJ) $(OBJ1)
	$(FC) $(OPT) $^ -o $@
%.o: %.f90
	$(FC) $(OPT) -c $^
.PHONY: clean
clean:
	-rm $(OBJ) *.mod *.a
