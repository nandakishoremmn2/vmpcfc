#include "grid.h"

#include <iostream>
#include <cstdlib>

int main (int argc, char *argv[]) 
{
	Grid grid(20, 10, 1e-4);

	grid.save("xi.dat");

	return 0;
}