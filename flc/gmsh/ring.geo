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

a = 0.4;
b = 1.0;
d_1 = 0.6;
d_2 = 0.8;

Point(1) = { 0, 0, 0};
Point(2) = { a / Sqrt(2), a / Sqrt(2), 0};
Point(3) = { d_1 / Sqrt(2), d_1 / Sqrt(2), 0};
Point(4) = { d_2 / Sqrt(2), d_2 / Sqrt(2), 0};
Point(5) = { b / Sqrt(2), b / Sqrt(2), 0};
Point(6) = { 0, a, 0};
Point(7) = { 0, d_1, 0};
Point(8) = { 0, d_2, 0};
Point(9) = { 0, b, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};

Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 9};

Circle(7) = {2, 1, 6};
Circle(8) = {3, 1, 7};
Circle(9) = {4, 1, 8};
Circle(10) = {5, 1, 9};

Line Loop(1) = {1, 8, -4, -7};
Plane Surface(1) = {1};

Line Loop(2) = {2, 9, -5, -8};
Plane Surface(2) = {2};

Line Loop(3) = {3, 10, -6, -9};
Plane Surface(3) = {3};

Q1[] = Rotate { {0,0,1}, {0,0,0},  -Pi/4}{ Duplicata{ Surface{:};} };
Q2[] = Rotate { {0,0,1}, {0,0,0}, -Pi/2}{ Duplicata{ Surface{:};} };
Q3[] = Rotate { {0,0,1}, {0,0,0}, -Pi}{ Duplicata{ Surface{:};} };

Physical Surface(100) = Surface{:};

// The physical surfaces are numbered in the software.

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

