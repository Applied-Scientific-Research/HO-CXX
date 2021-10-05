/*
 * MathHelpers.hpp - Useful mathematical helper functions
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mohammad Hajit <mhajit@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

cell_sides i2f[] = { west,east,south,north }; //(i2f[ijp] is the local side index (south, east, ... CCW) of the current cell: south =0:ibnd=0,idir=1; east=1:ibnd=1, idir=0, ...
cell_sides f2i[] = { north,east,west,south }; //(f2i[ijp] is the inverse of i2f, maps the senw to ijp
double sgn[] = { -1., 1., -1., 1. }; //to calculate the upwind flux from the two neighboring fluxes at the common face
double RK4_coeff[] = {0., 0.5, 0.5, 1.};  //the coefficients for the explicit RK4 time integration
double RK2_coeff[] = {0., 1.}; //the coefficients for the explicit RK2 time integration

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
		xn = -1.;
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

double DOT(Cmpnts2 a, Cmpnts2 b) {
	return (a.x * b.x + a.y * b.y);
}

template <typename T>
inline T sign(T x) { return ((x > 0) ? T(1) : T(-1)); }

int are_intersecting_infinite_line(Cmpnts2 point, double slope, Cmpnts2 node1, Cmpnts2 node2) { //checks if a ray casted from point with slope slope intersects the line node1-node2. it uses the original line equation
	// the drawback of the method is that it does not see the casted ray as half line, rather it sees an infinite line. SO the below approach just tells us if the infinite line intersects the segment or not
    // Every point (x,y) that is substituted in the ray equation (m(x-xpoint)+ypoint): 1) if zero: is on the ray, 2) if positive: is on one side and 3) negative: on the othet side of the ray
    //So lets insert node1 and node2 into the ray equation.
    double d1 = slope * (node1.x - point.x) + point.y - node1.y;
    double d2 = slope * (node2.x - point.x) + point.y - node2.y;

    // If d1 and d2 both have the same sign, they are both on the same side of the ray and in that case no intersection is possible.
    // if any of them is zero then the point is located on the node1-node2 line, which in our case, we can consider it to be inside the element.
	//note: since, I limited the slope to be between -1 and 1, it means the x component of the


	if (std::fabs(d1 * d2) < 1.e-6 || d1 * d2 < 0) return 1;
	else return 0;
}

int are_intersecting(Cmpnts2 point, double slope, Cmpnts2 node1, Cmpnts2 node2) { //checks if a ray casted from point with slope slope intersects the line node1-node2.
	//The casted ray moves from the point and goes in the slope direction. This is the difference between this method and the method are_intersecting_infinite_line above.
	// here I used the parametric equations for both the lines. Then solved the system to get the t and s. t is the coeff for ray and s is the coeff for segment:
	// vector(r) = vector(rp) + t.vector(V), where V is the directional vector for the ray. V = (deltaX, DeltaY). I pick DeltaX=1, so DeltaY=slope. For a valid intersection t>=0
	// Also, vector(u) = vector(node1) + s.vector(U), where U is node2-node1= (Xnode2-Xnode1, Ynode2-Ynode1). Then intersection is valid when 0=<s<=1.
	// The lines intersect when vector(r) = vector(u), which gives:
	// s = (slope*(Xnode1-Xp)+Yp-Ynode1)/(Ynode2-Ynode1 - slope*(Xnode2-Xnode1)); t = Xnode1-Xp + s * (Xnode2-Xnode1)
	Cmpnts2 Delta;
	Delta.subtract(node2, node1);
	double s = (slope * (node1.x - point.x) + point.y - node1.y) / (Delta.y - slope * Delta.x);
	double t = node1.x - point.x + s * Delta.x;

	if (t>-1.e-8 && s >-1.e-8 && s<1.+1.e-8) return 1;
	else return 0;
}
