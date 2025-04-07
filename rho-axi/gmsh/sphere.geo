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

d = 0.2;
b = 1.0;
a = 0.5;

Point(0) = { 0, 0, 0 };
Point(1) = { 0, d, 0 };
Point(2) = { 0, a, 0 };
Point(3) = { a/Sqrt(2), a/Sqrt(2), 0 };
Point(4) = { d, d, 0 };
Point(5) = { b/Sqrt(2), b/Sqrt(2), 0 };
Point(6) = { 0, b, 0 };
Point(7) = { d, 0, 0 };
Point(8) = { a, 0, 0 };
Point(9) = { b, 0, 0 };


Line(1) = {1, 2};
Line(2) = {3, 4};
Circle(3) = { 2, 0, 3};

Line(5) = {2, 6};
Line(6) = {5, 3};
Circle(7) = { 6, 0, 5};

Line(8) = {1, 4};

Line(9) = {7, 8};
Circle(10) = { 3, 0, 8};

Line(11) = {8, 9};
Circle(12) = { 5, 0, 9};

Line(13) = {4, 7};

Line(14) = {0, 1};

Line(15) = {7, 0};

Line Loop(1) = {1, 3, 2, -8};
Plane Surface(1) = {1};

Line Loop(2) = {5, 7, 6, -3};
Plane Surface(2) = {2};

Line Loop(3) = {-13, -2, 10, -9};
Plane Surface(3) = {3};

Line Loop(4) = {-10, -6, 12, -11};
Plane Surface(4) = {4};

Line Loop(5) = {14, 8, 13, 15};
Plane Surface(5) = {5};

Q1[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{:};}};

Physical Surface(100) = {1, 2, 3, 4, 16, 21, 26, 31, 5, 36};

Physical Line(1) = {7, 12, 23, 34};
Physical Line(0) = {1, 5, 17, 22, 14, 37};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

