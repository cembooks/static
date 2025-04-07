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

d1 = 0.1;
a = 1.2;
d2 = 2.0;
d3 = 3.0;

a3 = a / Sqrt(3);
d23 = d2 / Sqrt(3);
d33 = d3 / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { d1, d1, d1};
Point(3) = {-d1, d1, d1};

Point(4) = { a3, a3, a3};
Point(5) = {-a3, a3, a3};

Point(6) = { d23, d23, d23};
Point(7) = {-d23, d23, d23};

Point(8) = { d33, d33, d33};
Point(9) = {-d33, d33, d33};

Line(1) = {2, 4};
Line(2) = {4, 6};

Line(3) = {3, 5};
Line(4) = {5, 7};

Line(5) = {2, 3};

Circle(6) = {4, 1, 5};
Circle(7) = {6, 1, 7};

Line(8) = {6, 8};
Line(9) = {7, 9};

Circle(10) = {8, 1, 9};

Curve Loop(1) = {1, 6, -3, -5};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 7, -4, -6};
Plane Surface(2) = {2};

Curve Loop(3) = {8, 10, -9, -7};
Plane Surface(3) = {3};

Q1[] = Symmetry {1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};

Curve Loop(4) = {5, 15, -43, 28};
Surface(4) = {4};

Curve Loop(5) = {6, -13, 41, -26};
Surface(5) = {5};

Curve Loop(6) = {7, -18, 46, -31};
Surface(6) = {6};

Curve Loop(7) = {10, -23, 51, -36};
Surface(7) = {7};

Surface Loop(1) = {1, 11, 39, 24, 4, 5};
Volume(1) = {1};

Surface Loop(2) = {2, 16, 44, 29, 5, 6};
Volume(2) = {2};

Surface Loop(3) = {3, 21, 49, 34, 6, 7};
Volume(3) = {3};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1, 2, 3};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, 2, 3, V1[0], V1[1], V1[2]};}};
V3[] = Symmetry {-1, 0, 1, 0} {Duplicata {Volume{1, 2, 3};}};
V4[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3};}};

Surface Loop(4) = {4, 162, 255, 73, 344, 433};
Volume(4) = {4};

Physical Volume(1) = {Volume{:}};

Physical Surface(2) = {7, 411, 322, 500, 140, 229};

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;

