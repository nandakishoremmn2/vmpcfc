#ifndef GRID_H
#define GRID_H

typedef double real;

class Grid
{
public:

	Grid(real M, real gamma, real tau, int m, int n, real R_max, real nrtol, real tol); 
	~Grid();

	real get_residue();
	real get_delta(int i, int j, real coeff[20]);			// Calculates correction for Newton Raphson iterations
	real get_T(int i, int j);
	void calc_coefficients(int i, int j, real coeff[20]);	// Calculates A, B, C, D, H

	void sweep(int n); 				// Iterates over the grid once
	void minimize(int i, int j); 	// Executes newton raphson iteration on the (i, j)th coordinate

	void apply_boundary_conditions();

	void init_r_and_t();
	void init_xi();

	real *allocate(int m);
	real **allocate(int m, int n);

	void deallocate(real *var, int m);
	void deallocate(real **var, int m, int n);

	void save(char *filename);
	void load(char *filename);

	int get_nr();
	int get_nt();

protected:

	real **xi;
	real **temp;
	real **err; 	// No. of resets happened here

	real *r, *r1, *k;					// R
	int nr;  							// Size of r
	real *t, *t1, *h; 					// Theta
	int nt; 							// Size of t

	real r_max;

	// real A[4], B[4], C[4], D[4], H[4];

	real M_inf, M_inf2;
	real P_inf;

	real gamma, alpha;
	real lambda2, tau;

	real tol; // Tolerance for solution
	real NRtol; // Tolerance for newton raphson iteration

};

#endif // GRID_H