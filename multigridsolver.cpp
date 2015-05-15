#include "multigridsolver.h"
#include <cstdio>

MultigridSolver::MultigridSolver(real M, real Gamma, real Tau, int r_density, int t_density, int r_density_min, int t_density_min, real R_max, real nrtol, real Tol)
{
	Grid = new Multigrid(M, Gamma, Tau, r_density, t_density, r_density_min, t_density_min, R_max, nrtol, Tol);
}

MultigridSolver::~MultigridSolver()
{
	delete Grid;
}

void MultigridSolver::solve()
{
	get_initial_solution(Grid);
}

void MultigridSolver::get_initial_solution(Multigrid *grid)
{
	if (grid->grid2 == NULL)
	{
		grid->sweep(5000);
	}
	else
	{
		get_initial_solution(grid->grid2);
		grid->sweep(5000);
	}
}

void MultigridSolver::save(char *filename)
{
	Grid->save(filename);
}
