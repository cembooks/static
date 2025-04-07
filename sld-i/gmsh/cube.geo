/******************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

r = 5; //@1

d1 = 0.05;
a = 0.2;
b = 0.4;
d2 = 0.8;
d3 = 2.0;

a3 = a / Sqrt(3);
b3 = b / Sqrt(3);

Point(0) = { 0, 0, 0};

Point(1) = { d1, d1, d1};
Point(2) = { d1,-d1, d1};

Point(3) = { a3, a3, a3};
Point(4) = { a3,-a3, a3};

Point(5) = { b3, b3, b3};
Point(6) = { b3,-b3, b3};

Point(7) = { d2, d2, d2};
Point(8) = { d2,-d2, d2};

Point(9) = { d3, d3, d3};
Point(10) = { d3,-d3, d3};

Line(1) = {1, 3};
Line(2) = {3, 5};
Line(3) = {5, 7};
Line(4) = {7, 9};

Line(5) = {2, 4};
Line(6) = {4, 6};
Line(7) = {6, 8};
Line(8) = {8, 10};

Line(9) = {1, 2};

Circle(10) = {3, 0, 4};

Circle(11) = {5, 0, 6};

Line(12) = {7, 8};
Line(13) = {9, 10};

Line Loop(1) = {5, -10, -1, 9};
Plane Surface(1) = {1};

Line Loop(2) = {6, -11, -2, 10};
Plane Surface(2) = {2};

Line Loop(3) = {7, -12, -3, 11};
Plane Surface(3) = {3};

Line Loop(4) = {8, -13, -4, 12};
Plane Surface(4) = {4};

Q1[] = Symmetry {0, 1, 1, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {0, -1, 1, 0} {Duplicata {Surface{:};}};

ll = newll; Line Loop(ll) = {9, -18, 57, -37};
Plane Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {10, 16, -55, 35};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {11, 21, -60, 40};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {12, 26, -65, 45};
Plane Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {13, 31, -70, 50};
Plane Surface(ll) = {ll};

Surface Loop(1) = {1, 33, 53, 14, 71, 72};
Volume(1) = {1};

Surface Loop(2) = {2, 19, 58, 38, 72, 73};
Volume(2) = {2};

Surface Loop(3) = {3, 24, 63, 43, 73, 74};
Volume(3) = {3};

Surface Loop(4) = {4, 29, 68, 48, 74, 75};
Volume(4) = {4};

V1[] = Symmetry {1, 1, 0, 0} {Duplicata {Volume{:};}};
V2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Volume{:};}};
V3[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3, 4};}};
V4[] = Symmetry {1, 0, -1, 0} {Duplicata {Volume{1, 2, 3, 4};}};

Surface Loop(5) = {71, 97, 341, 217, 461, 581};
Volume(5) = {5};

Physical Volume(0) = {5};

Physical Volume(11) = {1};
Physical Volume(12) = {2};
Physical Volume(13) = {3};
Physical Volume(14) = {4};

Physical Volume(21) = {560};
Physical Volume(22) = {591};
Physical Volume(23) = {622};
Physical Volume(24) = {653};

Physical Volume(31) = {320};
Physical Volume(32) = {351};
Physical Volume(33) = {382};
Physical Volume(34) = {413};

Physical Volume(41) = {440};
Physical Volume(42) = {471};
Physical Volume(43) = {502};
Physical Volume(44) = {533};

Physical Volume(51) = {76};
Physical Volume(52) = {107};
Physical Volume(53) = {138};
Physical Volume(54) = {169};

Physical Volume(61) = {196};
Physical Volume(62) = {227};
Physical Volume(63) = {258};
Physical Volume(64) = {289};

// To apply the homogeneous Neumann boundary conditions on the vertical
// boundaries uncomment the following two lines ...

// Physical Surface(1) = {559, 679};
// Physical Surface(0) = {75, 439, 195, 315};

// ... and comment out the next line.

Physical Surface(0) = {75, 439, 195, 315};
Physical Surface(1) = {559, 679};

Recombine Surface "*";
Recombine Volume "*";

Transfinite Volume "*";
Transfinite Surface "*";
Transfinite Line "*" = 2*r;

Transfinite Line {1, 5, 80, 198, 17, 83, 210, 34} = r;
Transfinite Line {2, 6, 22, 39, 111, 119, 236, 229} = r;
Transfinite Line {3, 7, 27, 44, 142, 150, 267, 260} = r;
Transfinite Line {-4, -8, 32, -49, 173, -181, -291, 298} = r  Using Progression 0.85;

