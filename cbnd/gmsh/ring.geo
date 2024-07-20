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

a = 0.5;
b = 1.0;

Point(1) = { 0, 0, 0};
Point(2) = { a, 0, 0};
Point(3) = { b, 0, 0};
Point(4) = { 0, a, 0};
Point(5) = { 0, b, 0};

Line(1) = {2, 3};
Line(2) = {5, 4};

Circle(9) = {4, 1, 2};
Circle(10) = {3, 1, 5};

Line Loop(1) = {1, 9, 10, 2};
Plane Surface(1) = {1};

Q1[] = Symmetry {1, 0, 0, 0} {Duplicata {Surface{1};}};
Q2[] = Symmetry {0, 1, 0, 0} {Duplicata {Surface{1, Q1[0]};}};

Physical Surface(1) = {1, Q1[0], Q2[0], Q2[1]};
Physical Line(1) = {9, 15, 20, 25};
Physical Line(3) = {10, 13, 18, 23};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

