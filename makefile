CFLAGS = -O3 -Wno-write-strings
OPENMPFLAGS = -fopenmp
CC = g++
SRC = main.cpp grid.cpp multigrid.cpp multigridsolver.cpp
OBJ = $(SRC:.cpp = .o)

all: $(OBJ)
	$(CC) $(CFLAGS) -o pcf_solver $(OBJ)

openmp: $(OBJ)
	$(CC) $(CFLAGS) $(OPENMPFLAGS) -o pcf_solver_omp $(OBJ)

clean:
	rm -f core *.o a.out nohup.out