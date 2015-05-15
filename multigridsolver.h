#ifndef MULTIGRIDSOLVER_H
#define MULTIGRIDSOLVER_H

#include "multigrid.h"

class MultigridSolver
{
public:
	MultigridSolver(real M, real gamma, real tau, int r_density, int t_density, int r_density_min, int t_density_min, real R_max, real nrtol, real tol);
	~MultigridSolver();

	void solve();
	void save(char *filemane);

private:
	Multigrid *Grid;

	void get_initial_solution(Multigrid *grid);

};

#endif	