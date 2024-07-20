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

a = 0.4;
d1 = 0.6;
d2 = 0.8;
b = 1.0;

Point(0) = { 0, 0, 0 };
Point(1) = { 0, a, 0 };
Point(2) = { 0, d1, 0 };
Point(3) = { d1 / Sqrt(2), d1 / Sqrt(2), 0 };
Point(4) = { a  / Sqrt(2), a  / Sqrt(2), 0 };
Point(5) = { d2 / Sqrt(2), d2 / Sqrt(2), 0 };
Point(6) = { 0, d2, 0 };
Point(7) = { b / Sqrt(2), b / Sqrt(2), 0 };
Point(8) = { 0, b, 0 };

Line(1) = {1, 2};
Line(2) = {3, 4};
Circle(3) = { 2, 0, 3};
Circle(4) = { 4, 0, 1};

Line(5) = {2, 6};
Line(6) = {5, 3};
Circle(7) = { 6, 0, 5};

Line(8) = {6, 8};
Line(9) = {7, 5};
Circle(10) = { 8, 0, 7};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {6, 7, 5, -3};
Plane Surface(2) = {2};

Line Loop(3) = {9, 10, 8, -7};
Plane Surface(3) = {3};

Q1[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{:};}};

Physical Surface(100) = Surface{:};

Physical Line(1) = {4, 15, 30, 45};
Physical Line(3) = {10, 25, 40, 55};
Physical Line(0) = {1, 5, 8, 27, 34, 39};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

