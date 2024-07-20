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

a = 0.3;
b = 1.0;
d = 0.65;

Point(1) = { 0, 0, 0};
Point(2) = { a / Sqrt(2), a / Sqrt(2), 0};
Point(3) = { d / Sqrt(2), d / Sqrt(2), 0};
Point(4) = { b / Sqrt(2), b / Sqrt(2), 0};
Point(5) = { 0, a, 0};
Point(6) = { 0, d, 0};
Point(7) = { 0, b, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};

Line(3) = {5, 6};
Line(4) = {6, 7};

Circle(5) = {2, 1, 5};
Circle(6) = {3, 1, 6};
Circle(7) = {4, 1, 7};

Line Loop(1) = {-5, 1, 6, -3};
Plane Surface(1) = {1};

Line Loop(2) = {-6, 2, 7, -4};
Plane Surface(2) = {2};

Q1[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{1, 2};}};
Q2[] = Symmetry { 1, 0, 0, 0} {Duplicata {Surface{1, 2, Q1[0], Q1[1]};}};
Q3[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{1, 2, Q1[0], Q1[1], Q2[0],
Q2[1], Q2[2], Q2[3]};}};

Physical Surface(1) = Surface{:};

// The physical surfaces must be numbered in the software.

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

