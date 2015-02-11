#include <fstream>
#include <cmath>
#include <cstdlib>

#include "grid.h"

Grid::Grid(int m, int n, int nrtol)
{
	nr = m;
	nt = n;

	NRtol = nrtol;

	gamma = 1.4;
	alpha = gamma / ( gamma - 1 );

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

void Grid::save(char *filename)
{
	std::ofstream outfile(filename);
	for (int i = 0; i < nr; ++i)
	{
		for (int j = 0; j < nt; ++j)
		{
			outfile<<xi[i][j]<<" ";
		}
		outfile<<"\n";
	}
}

void Grid::sweep(int n)
{
	for (int num = 0; num < n; ++num)
	{
		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nt; ++j)
			{
				temp[i][j] = xi[i][j];
				minimize(i, j);
			}
		}
	}
}

void Grid::minimize(int i, int j)
{
	real delta;
	int MAX_NR_ITER = 100;
	int iter = 0;
	do
	{
		delta = get_delta(i, j);
		xi[i][j] = xi[i][j] - delta;

		iter++;

	} while (delta > NRtol && iter < MAX_NR_ITER);
}

real Grid::get_delta(int i, int j)
{
	real g = 0., g_ = 0.; 	// g = dJ/dXij
	real t1, t2;	// Temp variables
	for (int k = 0; k < 4; ++k)
	{
		t1 = ( A[k] * xi[i][j] + B[k] ) * xi[i][j] + C[k];
		t2 = 2. * A[k] * xi[i][j] + B[k];

		g += ( alpha * pow(t1, alpha-1) * t2 + D[k] ) * H[k];
		g_ += ( 2 * A[k] * alpha * pow(t1, alpha - 1) + \
			alpha * (alpha-1) * pow(t1, alpha-2) * t2 ) * H[k];
	}
	return g/g_;
}