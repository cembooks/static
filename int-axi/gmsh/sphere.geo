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

r = 7; //@1

a = 0.3;
b = 1.0;
d = 0.65;
h = (d-a)/2;

Point(0) = { 0, 0, 0 };

Point(1) = { 0, a, 0 };
Point(2) = { 0, d, 0 };
Point(3) = { 0, b, 0 };

Point(4) = { a / Sqrt(2), a / Sqrt(2), 0 };
Point(5) = { d / Sqrt(2), d / Sqrt(2), 0 };
Point(6) = { b / Sqrt(2), b / Sqrt(2), 0 };

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};

Circle(5) = { 1, 0, 4};
Circle(6) = { 2, 0, 5};
Circle(7) = { 3, 0, 6};

Line Loop(1) = {5, 3, -6, -1};
Plane Surface(1) = {1};

Line Loop(2) = {6, 4, -7, -2};
Plane Surface(2) = {2};

Q1[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{:};}};

Physical Surface(100) = Surface{:};

Physical Line(1) = {5, 9, 19, 29};
Physical Line(3) = {7, 16, 26, 36};
Physical Line(0) = {2, 1, 22, 27};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

