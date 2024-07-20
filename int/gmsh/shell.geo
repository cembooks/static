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

r = 3; //@1

a = 0.3;
b = 1.0;
d = 0.65;

a2 = a / Sqrt(2);
b2 = b / Sqrt(2);
d2 = d / Sqrt(2);

a3 = a / Sqrt(3);
b3 = b / Sqrt(3);
d3 = d / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { a3, a3, a3};
Point(3) = { d3, d3, d3};
Point(4) = { b3, b3, b3};

Point(5) = { 0, a, 0};
Point(6) = { 0, d, 0};
Point(7) = { 0, b, 0};

Point(8) = { 0, a2, a2};
Point(9) = { 0, d2, d2};
Point(10) = { 0, b2, b2};

Point(11) = { a2, a2, 0};
Point(12) = { d2, d2, 0};
Point(13) = { b2, b2, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};

Line(3) = {5, 6};
Line(4) = {6, 7};

Line(5) = {8, 9};
Line(6) = {9, 10};

Line(7) = {11, 12};
Line(8) = {12, 13};

Circle(9) = {2, 1, 11};
Circle(10) = {11, 1, 5};
Circle(11) = {5, 1, 8};
Circle(12) = {8, 1, 2};

Circle(13) = {3, 1, 12};
Circle(14) = {12, 1, 6};
Circle(15) = {6, 1, 9};
Circle(16) = {9, 1, 3};

Circle(17) = {4, 1, 13};
Circle(18) = {13, 1, 7};
Circle(19) = {7, 1, 10};
Circle(20) = {10, 1, 4};

Curve Loop(1) = {9, 10, 11, 12};
Surface(1) = {1};

Curve Loop(2) = {13, 14, 15, 16};
Surface(2) = {2};

Curve Loop(3) = {17, 18, 19, 20};
Surface(3) = {3};

Curve Loop(4) = {12, 1, -16, -5};
Plane Surface(4) = {4};

Curve Loop(5) = {16, 2, -20, -6};
Plane Surface(5) = {5};

Curve Loop(6) = {9, 7, -13, -1};
Plane Surface(6) = {6};

Curve Loop(7) = {13, 8, -17, -2};
Plane Surface(7) = {7};

Curve Loop(8) = {10, 3, -14, -7};
Plane Surface(8) = {8};

Curve Loop(9) = {14, 4, -18, -8};
Plane Surface(9) = {9};

Curve Loop(10) = {11, 5, -15, -3};
Plane Surface(10) = {10};

Curve Loop(11) = {15, 6, -19, -4};
Plane Surface(11) = {11};

Surface Loop(1) = {4, 6, 8, 10, 1, 2};
Volume(1) = {1};

Surface Loop(2) = {5, 7, 9, 11, 2, 3};
Volume(2) = {2};

V1[] = Symmetry {1, 0, 0, 0} {Duplicata {Volume{1, 2};}};
V2[] = Symmetry {0, 0, 1, 0} {Duplicata {Volume{1, 2, V1[0], V1[1]};}};
V3[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, 2, V1[0], V1[1], V2[0],
V2[1], V2[2], V2[3]};}};
V4[] = Symmetry {0, 1,-1, 0} {Duplicata {Volume{1, 2, V1[0], V1[1], V2[0],
V2[1], V2[2], V2[3], V3[0], V3[1], V3[2], V3[3], V3[4], V3[5], V3[6], V3[7]};}};
V5[] = Symmetry {1, 1, 0, 0} {Duplicata {Volume{1, 2, V1[0], V1[1], V2[0],
V2[1], V2[2], V2[3]};}};
V6[] = Symmetry {1, -1, 0, 0} {Duplicata {Volume{1, 2, V1[0], V1[1], V2[0],
V2[1], V2[2], V2[3]};}};

Physical Volume(1) = {Volume{:}} ;

// The physical surfaces must be numbered in the software.

// Rotate the shell to get exact approximation on the x axis. This is
// needed for projecting the potential on the segment between [a,0,0]
// and [b,0,0].
V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

Recombine Surface "*";
Recombine Volume "*";

Transfinite Surface "*";
Transfinite Volume "*";

Transfinite Line "*" = r;

