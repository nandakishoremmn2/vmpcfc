#include "grid.h"

#include <iostream>
#include <cstdlib>

int main (int argc, char *argv[]) 
{
	Grid grid(90, 410, 20, 1e-100);

	grid.load("xi.dat");

	grid.sweep(5000);

	grid.save("xi.dat");

	return 0;
}