#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "grid.h"

class Multigrid : public Grid
{
public:
	Multigrid(int r_density, int t_density, int r_density_min, int t_density_min, real R_max, real nrtol);
	~Multigrid();

	void interpolate();
	void restrict();

	Multigrid *grid2;

	int level;
	
};

#endif	