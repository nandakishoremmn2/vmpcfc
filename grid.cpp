#include <fstream>
#include <cmath>
#include <cstdlib>

#include "grid.h"

Grid::Grid(real M, real Gamma, real Tau, int m, int n, real R_max, real nrtol, real Tol)
{
	nr = m;
	nt = n;

	r_max = R_max;

	NRtol = nrtol;
	tol = Tol;

	gamma = Gamma;
	alpha = gamma / ( gamma - 1 );

	tau = Tau;
	lambda2 = ( 1 - tau ) / ( 1 + tau );

	M_inf = M;
	M_inf2 = pow(M, 2);

	xi = allocate(nt, nr);
	temp = allocate(nt, nr);
	err = allocate(nt, nr);

	r = allocate(nr);
	r1 = allocate(nr-1);
	k = allocate(nr-1);

	t = allocate(nt);
	t1 = allocate(nt-1);
	h = allocate(nt-1);

	init_r_and_t();
	init_xi();

}

Grid::~Grid()
{
	deallocate(xi, nt, nr);
	deallocate(temp, nt, nr);

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

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			var[i][j] = 0;
		}
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
	for (int i = 0; i < nt; ++i)
	{
		for (int j = 0; j < nr; ++j)
		{
			outfile<<xi[i][j]<<" ";
		}
		outfile<<"\n";
	}
	outfile.close();

	std::ofstream coordsfile("coords.dat");
	for (int i = 0; i < nt; ++i)
	{
		coordsfile<<t[i]<<" ";
	}
	coordsfile<<"\n";
	for (int j = 0; j < nr; ++j)
	{
		coordsfile<<r[j]<<" ";
	}
	coordsfile<<"\n";
	coordsfile<<lambda2;
	coordsfile<<"\n";
	coordsfile<<M_inf2;
	coordsfile.close();
}

void Grid::load(char *filename)
{
	std::ifstream infile(filename);
	for (int i = 0; i < nt; ++i)
	{
		for (int j = 0; j < nr; ++j)
		{
			infile>>xi[i][j];
		}
	}
}

void Grid::sweep(int n)
{
	real res;
	for (int num = 0; num < n; ++num)
	{
		apply_boundary_conditions();
		for (int i = 1; i < nt-1; ++i)
		{
			for (int j = 1; j < nr-1; ++j)
			{
				temp[i][j] = xi[i][j];
				minimize(i, j);
			}
		}
		res = get_residue();
		if( res < tol && num > 5)break;
	}
}

void Grid::minimize(int i, int j)
{
	real delta;
	int MAX_NR_ITER = 10;
	int iter = 0;

	real coeff[20];
	// This need to be done only once prior to NR iterations
	calc_coefficients(i, j, coeff);
	do
	{
		delta = get_delta(i, j, coeff);
		if(delta!=delta)
		{
			// err[i][j]++;
			break;
		}
		
		xi[i][j] = xi[i][j] - delta;

		iter++;

	} while (abs(delta) > NRtol || iter < MAX_NR_ITER);
}

real Grid::get_delta(int i, int j, real coeff[20])
{
	real A[4], B[4], C[4], D[4], H[4];

	for (int k = 0; k < 4; ++k)
	{
		A[k] = coeff[5*k];
		B[k] = coeff[5*k+1];
		C[k] = coeff[5*k+2];
		D[k] = coeff[5*k+3];
		H[k] = coeff[5*k+4];
	}

	real g = 0., g_ = 0.; 	// g = dJ/dXij
	real t1, t2;	// Temp variables
	real s=0.;

	for (int k = 0; k < 4; ++k)
	{
		t1 = ( A[k] * xi[i][j] + B[k] ) * xi[i][j] + C[k];
		t2 = 2. * A[k] * xi[i][j] + B[k];

		g += ( alpha * pow(t1, alpha-1) * t2 + D[k] ) * H[k];
		g_ += ( 2 * A[k] * alpha * pow(t1, alpha - 1) + alpha * (alpha-1) * pow(t1, alpha-2) * t2 ) * H[k];
	}
	return g/g_;
}

real Grid::get_T(int i, int j)
{
	return ( pow(r1[j]*r1[j] + lambda2, 2) - 4*lambda2*r1[j]*r1[j]*pow(cos(t1[i]), 2) ) / pow(r1[j], 4);
}

void Grid::calc_coefficients(int i, int j, real coeff[20])
{
	real A[4], B[4], C[4], D[4], H[4];
	real kk, hh, rr, tt, T, x1, x2;

	T = get_T(i, j); // T^2 not T
	rr = r1[j];
	kk = k[j];
	tt = t1[i];
	hh = h[i];
	x1 = xi[i][j+1] + xi[i+1][j+1] - xi[i+1][j];
	x2 = xi[i+1][j] + xi[i+1][j+1] - xi[i][j+1];

	A[0] = -.5 * (gamma-1) * M_inf2 / T * ( 1/(kk*kk) + 1/(rr*rr*hh*hh) ) / 4;
	B[0] = .5 * (gamma-1) * M_inf2 / T * ( cos(tt)/kk - sin(tt)/(rr*hh) 
		+ x1 / (2*kk*kk) + x2 / (2*hh*hh*rr*rr) 
		);
	C[0] = 1 + .5 * (gamma-1) * M_inf2 / T * ( T - 1 
		- cos(tt) * x1 / kk + sin(tt) * x2 / ( rr * hh ) 
		- x1*x1 / ( 4*kk*kk ) - x2*x2 / ( 4*rr*rr*hh*hh ) 
		);
	D[0] = gamma * M_inf2 / T * ( \
		- .5 * ( rr*rr - 1 ) / ( rr*rr ) * cos(tt) / kk
		+ .5 * ( rr*rr - 1 ) / ( rr*rr*rr ) * sin(tt) / hh
		);
	H[0] = rr*T*hh*kk;


	T = get_T(i-1, j);
	rr = r1[j];
	kk = k[j];
	tt = t1[i-1];
	hh = h[i-1];
	x1 = xi[i-1][j+1] + xi[i][j+1] - xi[i-1][j];
	x2 = xi[i][j+1] - xi[i-1][j] - xi[i-1][j+1];

	A[1] = -.5 * (gamma-1) * M_inf2 / T * ( 1/(kk*kk) + 1/(rr*rr*hh*hh) ) / 4;
	B[1] = .5 * (gamma-1) * M_inf2 / T * ( cos(tt)/kk + sin(tt)/(rr*hh) 
		+ x1 / (2*kk*kk) - x2 / (2*hh*hh*rr*rr) 
		);
	C[1] = 1 + .5 * (gamma-1) * M_inf2 / T * ( T - 1 
		- cos(tt) * x1 / kk + sin(tt) * x2 / ( rr * hh ) 
		- x1*x1 / ( 4*kk*kk ) - x2*x2 / ( 4*rr*rr*hh*hh ) 
		);
	D[1] = gamma * M_inf2 / T * ( 
		- .5 * ( rr*rr - 1 ) / ( rr*rr ) * cos(tt) / kk
		- .5 * ( rr*rr - 1 ) / ( rr*rr*rr ) * sin(tt) / hh
		);
	H[1] = rr*T*hh*kk;

	T = get_T(i, j-1);
	rr = r1[j-1];
	kk = k[j-1];
	tt = t1[i];
	hh = h[i];
	x1 = xi[i+1][j] - xi[i][j-1] - xi[i+1][j-1];
	x2 = xi[i+1][j-1] + xi[i+1][j] - xi[i][j-1];

	A[2] = -.5 * (gamma-1) * M_inf2 / T * ( 1/(kk*kk) + 1/(rr*rr*hh*hh) ) / 4;
	B[2] = .5 * (gamma-1) * M_inf2 / T * ( -cos(tt)/kk - sin(tt)/(rr*hh) 
		- x1 / (2*kk*kk) + x2 / (2*hh*hh*rr*rr) 
		);
	C[2] = 1 + .5 * (gamma-1) * M_inf2 / T * ( T - 1 
		- cos(tt) * x1 / kk + sin(tt) * x2 / ( rr * hh ) 
		- x1*x1 / ( 4*kk*kk ) - x2*x2 / ( 4*rr*rr*hh*hh ) 
		);
	D[2] = gamma * M_inf2 / T * ( 
		+ .5 * ( rr*rr - 1 ) / ( rr*rr ) * cos(tt) / kk
		+ .5 * ( rr*rr - 1 ) / ( rr*rr*rr ) * sin(tt) / hh
		);
	H[2] = rr*T*hh*kk;

	T = get_T(i-1, j-1);
	rr = r1[j-1];
	kk = k[j-1];
	tt = t1[i-1];
	hh = h[i-1];
	x1 = xi[i-1][j] - xi[i-1][j-1] - xi[i][j-1];
	x2 = xi[i][j-1] - xi[i-1][j-1] - xi[i-1][j];

	A[3] = -.5 * (gamma-1) * M_inf2 / T * ( 1/(kk*kk) + 1/(rr*rr*hh*hh) ) / 4;
	B[3] = .5 * (gamma-1) * M_inf2 / T * ( -cos(tt)/kk + sin(tt)/(rr*hh) 
		- x1 / (2*kk*kk) - x2 / (2*hh*hh*rr*rr) 
		);
	C[3] = 1 + .5 * (gamma-1) * M_inf2 / T * ( T - 1 
		- cos(tt) * x1 / kk + sin(tt) * x2 / ( rr * hh ) 
		- x1*x1 / ( 4*kk*kk ) - x2*x2 / ( 4*rr*rr*hh*hh ) 
		);
	D[3] = gamma * M_inf2 / T * ( 
		+ .5 * ( rr*rr - 1 ) / ( rr*rr ) * cos(tt) / kk
		- .5 * ( rr*rr - 1 ) / ( rr*rr*rr ) * sin(tt) / hh
		);
	H[3] = rr*T*hh*kk;

	for (int k = 0; k < 4; ++k)
	{
		coeff[5*k] = A[k];
		coeff[5*k+1] = B[k];
		coeff[5*k+2] = C[k];
		coeff[5*k+3] = D[k];
		coeff[5*k+4] = H[k];
	}

}

void Grid::init_r_and_t()
{
	// Init R
	real dr1 = ( 1. - 1./r_max ) / ( nr - 1 );
	real dr = 1.;

	for (int i = 0; i < nr; ++i)
	{
		r[i] = 1./dr;
		dr -= dr1;
	}
	for (int i = 0; i < nr-1; ++i)
	{
		k[i] = r[i+1] - r[i];
		r1[i] = r[i] + .5*k[i];
	}

	// Init theta
	real t_min = M_PI*0, t_max = M_PI;
	real dt = ( t_max - t_min ) / ( nt - 3 );

	for (int i = 1; i < nt-1; ++i)
	{
		t[i] = t_min + ( i-1 ) * dt;
	}
	for (int i = 1; i < nt-2; ++i)
	{
		h[i] = t[i+1] - t[i];
	}

	// Adding extra Theta boundary values on both ends
	h[0] = h[1];
	t[0] = t[1] - h[0];
	h[nt-2] = h[nt-3];
	t[nt-1] = t[nt-2] + h[nt-2];

	for (int i = 0; i < nt-1; ++i)
	{
		t1[i] = t[i] + .5*h[i];
	}

}

void Grid::apply_boundary_conditions()
{
	for (int j = 0; j < nr; ++j)
	{
		xi[0][j] = xi[2][j];
		xi[nt-1][j] = xi[nt-3][j];
	}
	for (int i = 0; i < nt; ++i)
	{
		xi[i][nr-1] = 1/r[nr-1]*cos(t[i]);
		xi[i][0] = ( ( pow(k[0] + k[1], 2) * xi[i][1] - k[0]*k[0]*xi[i][2] ) / k[1] \
			+ k[0] * ( k[0] + k[1] ) * -cos(t[i]) 
			) / ( k[1] + 2*k[0] );
	}
}

real abs(real num)
{
	return num > 0 ? num : 0;
}

real Grid::get_residue()
{
	real tmp = 0;
	for (int i = 1; i < nt-1; ++i)
	{
		for (int j = 1; j < nr-1; ++j)
		{
			if(abs(xi[i][j]-temp[i][j]) > tmp)
			{
				tmp = abs(xi[i][j]-temp[i][j]);
			}
		}
	}
	return tmp;
}

void Grid::init_xi()
{
	for (int i = 0; i < nt; ++i)
	{
		for (int j = 0; j < nr; ++j)
		{
			xi[i][j] = cos(t[i])/r[j];
		}
	}
}

int Grid::get_nr()
{
	return nr;
}

int Grid::get_nt()
{
	return nt;
}