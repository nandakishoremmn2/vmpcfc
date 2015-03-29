#include "grid.h"

#include <iostream>
#include <cstdlib>
#include <csignal>

Grid *grid_ptr;

int *get_input(int argc, char *argv[]);
void signalHandler( int signum );

int main (int argc, char *argv[]) 
{
	Grid grid(41, 17, 30, 1e-5);
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
