#include "multigrid.h"
#include <cstdio>
#include <cmath>

real min(real a, real b)
{
	return a < b ? a : b;
}

Multigrid::Multigrid(int r_density, int t_density, int r_density_min, int t_density_min, real R_max, real nrtol)
	:Grid(pow(2, r_density) + 1, pow(2, t_density) + 1, R_max, nrtol)
{

	level = min(r_density - r_density_min, t_density - t_density_min);

	if (level == 0)
	{
		grid2 = NULL;
	}
	else
	{
		grid2 = new Multigrid(r_density-1, t_density-1, r_density_min, t_density_min, R_max, nrtol);
	}
}

Multigrid::~Multigrid()
{
	delete grid2;
}

void Multigrid::interpolate()
{
	int nr2 = grid2->get_nr();
	int nt2 = grid2->get_nt();
	real **xi2 = grid2->xi;


	for (int i = 0; i < nt2-1; ++i)
	{
		for (int j = 0; j < nr2-1; ++j)
		{
			xi[2*i+1][2*j+1] = ( xi2[i][j] + xi2[i+1][j] + xi2[i][j+1] + xi2[i+1][j+1] ) / 4.0;
		}
	}

	for (int i = 0; i < nt2-1; ++i)
	{
		for (int j = 0; j < nr2; ++j)
		{
			xi[2*i+1][2*j] = ( xi2[i][j] + xi2[i+1][j] ) / 2.0;
		}
	}

	for (int i = 0; i < nt2; ++i)
	{
		for (int j = 0; j < nr2-1; ++j)
		{
			xi[2*i][2*j+1] = ( xi2[i][j] + xi2[i][j+1] ) / 2.0;
		}
	}

	for (int i = 0; i < nt2; ++i)
	{
		for (int j = 0; j < nr2; ++j)
		{
			xi[2*i][2*j] = xi2[i][j];
		}
	}
}

void Multigrid::restrict()
{
	
}