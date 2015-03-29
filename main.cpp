#include "multigrid.h"

#include <iostream>
#include <cstdlib>
#include <csignal>

Grid *grid_ptr;

int *get_input(int argc, char *argv[]);
void signalHandler( int signum );

int main(int argc, char const *argv[])
{
	Multigrid grid(.2, 1.4, 1, 6, 4, 5, 3, 30, 1e-5, 1e-4);
	grid_ptr = &grid;

	signal(SIGINT, signalHandler);

	grid.grid2->sweep(5000);
	// grid.grid2->save("xi.dat");
	grid.interpolate();
	grid.sweep(5000);

	grid.save("xi.dat");

	return 0;
}

int main2 (int argc, char *argv[]) 
{
	Grid grid(.2, 1.4, 1, 41, 17, 30, 1e-5, 1e-4);
	grid_ptr = &grid;

	signal(SIGINT, signalHandler);

	// grid.load("xi.dat");

	grid.sweep(5000);
	// grid.sweep(1);

	grid.save("xi.dat");

	return 0;
}

void signalHandler( int signum )
{
	grid_ptr->save("xi.dat");
	exit(0);  
}
