Point(1) = {0, 0.5, 0, 1.0};
//+
Point(2) = {4, 0.5, 0, 1.0};
//+
Point(3) = {4, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {4, 0, 0, 1.0};
//+
Point(6) = {30, 0, 0, 1.0};
//+
Point(7) = {30, 1, 0, 1.0};
//+
Point(8) = {30, 0.5, 0, 1.0};
//+
Point(9) = {30, 0.5, 0, 1.0};
//+
Point(10) = {6, 1, 0, 1.0};
//+
Point(11) = {6, 0.5, 0, 1.0};
//+
Point(12) = {6, 0, 0, 1.0};
//+
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {2, 3};
//+
Line(5) = {2, 5};
//+
Line(6) = {5, 12};
//+
Line(7) = {12, 11};
//+
Line(8) = {11, 2};
//+
Line(9) = {11, 10};
//+
Line(10) = {10, 3};
//+
Line(11) = {12, 6};
//+
Line(12) = {6, 8};
//+
Line(13) = {8, 11};
//+
Line(14) = {8, 7};
//+
Line(15) = {7, 10};
//+
//+
Curve Loop(1) = {11, 12, 13, -7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, 15, -9, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -4, -8, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, 7, 8, 5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, 3, 1, 4};
//+
Plane Surface(5) = {5};
//+
//+
Transfinite Surface {5} = {1, 2, 3, 4};
//+
Transfinite Surface {4} = {5, 12, 11, 2};
//+
Transfinite Surface {3} = {2, 11, 10, 3};
//+
Transfinite Surface {1} = {12, 6, 8, 11};
//+
Transfinite Surface {2} = {11, 8, 7, 10};
//+
Recombine Surface {5, 3, 4, 2, 1};
//+ Used for K = 5+
Transfinite Curve {3, 4, 9, 14, 5, 7, 12} = 6 Using Progression 1;
//+
Transfinite Curve {-2, 1} = 6 Using Progression .75;
//+
Transfinite Curve {10, 8, -6} = 6 Using Progression 1;
//+
Transfinite Curve {15, 13, -11} = 13 Using Progression .8;
//+ Used for K = 3+
//Transfinite Curve {3, 4, 9, 14, 5, 7, 12} = 9 Using Progression 1;
//+
//Transfinite Curve {-2, 1} = 9 Using Progression .85;
//+
//Transfinite Curve {10, 8, -6} = 9 Using Progression 1;
//+
//Transfinite Curve {15, 13, -11} = 21 Using Progression .875;
//+
Extrude {0, 0, 1} {
  Surface{1, 2, 3, 4, 5};
  Layers{1};
  Recombine;
}
//+
Physical Surface("Inlet") = {116};
//+
Physical Surface("Outlet") = {46, 28};
//+
Physical Surface("Top") = {50, 68, 112};
//+
Physical Surface("Bottom_Out") = {24, 90};
//+
Physical Surface("Bottom_In") = {120};
//+
Physical Surface("Step") = {102};
//+
Physical Surface("Periodic_Back") = {1, 2, 4, 3, 5};
//+
Physical Surface("Periodic_Front") = {59, 37, 81, 103, 125};
