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

d1 = 0.25;
a = 1.0;

Point(0) = { 0, 0, 0};
Point(1) = { d1, 0 , 0};
Point(2) = { d1, d1, 0};
Point(3) = { 0, d1, 0};

Point(4) = {a, 0, 0};
Point(5) = {a/Sqrt(2), a/Sqrt(2), 0};
Point(6) = { 0, a, 0};

Line(1) = {1, 4};
Line(2) = {2, 5};
Line(3) = {3, 6};

Line(4) = {1, 2};
Line(5) = {2, 3};

Line(6) = {0, 1};
Line(7) = {3, 0};

Circle(8) = {4, 0, 5};
Circle(9) = {5, 0, 6};

Line Loop(1) = {1, 8, -2, -4};
Plane Surface(1) = {1};

Line Loop(2) = {2, 9, -3, -5};
Plane Surface(2) = {2};

Line Loop(3) = {4, 5, 7, 6};
Plane Surface(3) = {3};

Q1[] = Symmetry {1, 0, 0, 0} {Duplicata {Surface{1, 2, 3};}};
Q2[] = Symmetry {0, 1, 0, 0} {Duplicata {Surface{1, 2, 3, Q1[0], Q1[1], Q1[2]};}};

Physical Surface(1) = { 1, 2, 3, Q1[0], Q1[1], Q1[2],
Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Q2[5]};

Physical Line(1) = {8, 9, 17, 12};
Physical Line(2) = {27, 32, 47, 42};

// Physical Line(1) = {8, 9, 17, 12, 27, 32, 47, 42};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

