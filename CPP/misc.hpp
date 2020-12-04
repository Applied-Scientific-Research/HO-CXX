#pragma once

cell_sides i2f[] = { west,east,south,north }; //(i2f[ijp] is the local side index (south, east, ... CCW) of the current cell: south =0:ibnd=0,idir=1; east=1:ibnd=1, idir=0, ...
unsigned int nbr[] = { east, south, west, north }; //ijpm = nbr[ijp]: neighbor of ijp face, returned in the same ijp ordering
double sgn[] = { -1., 1., -1., 1. }; //to calculate the upwind flux fronm the two neighboring fluxes at the common face

inline double minus_one_to_power(int n) {
	return (n % 2 == 0 ? 1. : -1.);
}

double Legendre(int n, double x) {
	//evaluates the order n Legendre polynomial at point xsi = x
	double LegendreP;
	double xn, xnm1, xnm2;
	if (!n) LegendreP = 1.;
	else if (n == 1) LegendreP = x;
	else {
		xnm2 = 1.;
		xnm1 = x;
		for (int i = 2; i <= n; ++i) {
			xn = ((2. * i - 1.) * x * xnm1 - (i - 1.) * xnm2) / i;
			xnm2 = xnm1;
			xnm1 = xn;
		}	
		LegendreP = xn;
	}
	return LegendreP;
}

void Gauss_solver(int n, double** LHS, double* RHS, double* x) {
	//Gauss Elimination solver without pivoting solves the system [LHS]*{x} = {RHS} to find x. The size of x and Rhs are n;  LHS is n by n
	double c;
	for (int k = 0; k < n - 1; ++k) {
		for (int i = k + 1; i < n; ++i) {
			c = LHS[i][k] / LHS[k][k];
			LHS[i][k] = 0.;
			RHS[i] = RHS[i] - c * RHS[k];
			for (int j = k + 1; j < n; ++j)
				LHS[i][j] -= c * LHS[k][j];
		}
	}

	//Back - substitution
	x[n - 1] = RHS[n - 1] / LHS[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		c = 0.;
		for (int j = i + 1; j < n; ++j)
			c += LHS[i][j] * x[j];
		x[i] = (RHS[i] - c) / LHS[i][i];
	}
}