#include <fstream>
#include <cmath>
#include <cstdlib>

#include "grid.h"

Grid::Grid(int m, int n, int nrtol)
{
	nr = m;
	nt = n;

	NRtol = nrtol;

	xi = allocate(nr, nt);
	temp = allocate(nr, nt);

	r = allocate(nr);
	r1 = allocate(nr-1);
	k = allocate(nr-1);

	t = allocate(nt);
	t1 = allocate(nt-1);
	h = allocate(nt-1);

}

Grid::~Grid()
{
	deallocate(xi, nr, nt);
	deallocate(temp, nr, nt);

	deallocate(r, nr);
	deallocate(r1, nr-1);
	deallocate(k, nr-1);

	deallocate(t, nt);
	deallocate(t1, nt-1);
	deallocate(h, nt-1);

}

real *Grid::allocate(int m)
{
	real *var = new real[m];
	return var;
}
real **Grid::allocate(int m, int n)
{
	real **var = new real*[m];
	for (int i = 0; i < m; ++i)
	{
		var[i] = new real[n];
	}
	return var;
}

void Grid::deallocate(real *var, int m)
{
	delete [] var;
}
void Grid::deallocate(real **var, int m, int n)
{
	for (int i = 0; i < m; ++i)
	{
		delete [] var[i];
	}
	delete [] var;
}