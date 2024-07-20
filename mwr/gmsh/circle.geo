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

d1 = 0.2;
a = 0.5;
b = 1.0;

Point(0) = { 0, 0, 0};
Point(1) = { d1, 0 , 0};
Point(2) = { d1, d1, 0};
Point(3) = { 0, d1, 0};

Point(4) = {a, 0, 0};
Point(5) = {a/Sqrt(2), a/Sqrt(2), 0};
Point(6) = { 0, a , 0};

Point(7) = { b, 0, 0};
Point(8) = {b/Sqrt(2), b/Sqrt(2), 0};
Point(9) = { 0, b, 0};

Line(1) = {1, 4};
Line(2) = {4, 7};

Line(3) = {2, 5};
Line(4) = {5, 8};

Line(5) = {3, 6};
Line(6) = {6, 9};

Line(7) = {1, 2};
Line(8) = {2, 3};

Circle(9) = {4, 0, 5};
Circle(10) = {5, 0, 6};

Circle(11) = {7, 0, 8};
Circle(12) = {8, 0, 9};

Line(13) = {0, 1};
Line(14) = {3, 0};

Line Loop(1) = {1, 9, -3, -7};
Plane Surface(1) = {1};

Line Loop(2) = {2, 11, -4, -9};
Plane Surface(2) = {2};

Line Loop(3) = {3, 10, -5, -8};
Plane Surface(3) = {3};

Line Loop(4) = {4, 12, -6, -10};
Plane Surface(4) = {4};

Line Loop(5) = {13, 7, 8, 14};
Plane Surface(5) = {5};

Q1[] = Symmetry {1, 0, 0, 0} {Duplicata {Surface{1, 2, 3, 4, 5};}};
Q2[] = Symmetry {0, 1, 0, 0} {Duplicata {Surface{1, 2, 3, 4, 5, Q1[0], Q1[1],
Q1[2], Q1[3], Q1[4]};}};

Physical Surface(1) = Surface{:};

//Physical Surface(1) = {1, 3, 25, 15, 62, 72, 47, 37, 5, 35, 82, 57, 5};
//Physical Surface(2) = {2, 4, 30, 20, 67, 77, 52, 42};

Physical Line(1) = {11, 12, 32, 22, 44, 54, 79, 69};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

