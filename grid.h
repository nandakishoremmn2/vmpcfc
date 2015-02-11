#ifndef GRID_H
#define GRID_H

typedef double real;

class Grid
{
public:

	Grid(int m, int n, int nrtol); 
	~Grid();

	real get_delta(int i, int j);			// Calculates correction for Newton Raphson iterations
	void calc_coefficients(int i, int j);	// Calculates A, B, C, D, H

	void sweep(int n); 				// Iterates over the grid once
	void minimize(int i, int j); 	// Executes newton raphson iteration on the (i, j)th coordinate

	real *allocate(int m);
	real **allocate(int m, int n);

	void deallocate(real *var, int m);
	void deallocate(real **var, int m, int n);

private:

	real **xi;
	real **temp;

	real *r, *r1, *k;					// R
	int nr;  							// Size of r
	real *t, *t1, *h; 					// Theta
	int nt; 							// Size of t

	real A[4], B[4], C[4], D[4], H[4];

	real M_inf;
	real P_inf;

	real NRtol; // Tolerance for newton raphson iteration

};

#endif // GRID_H