#pragma once
inline double minus_one_to_power(int n) {
	return (n % 2 == 0 ? 1. : -1);
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

