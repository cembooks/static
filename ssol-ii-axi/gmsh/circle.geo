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
b = 0.6;
d1 = 0.1;
d2 = 1.0;
d3 = 2.0;

Point(0) = { 0, 0, 0 };
Point(1) = { 0, d1, 0 };
Point(2) = { 0, a, 0 };
Point(3) = { a/Sqrt(2), a/Sqrt(2), 0 };
Point(4) = { d1, d1, 0 };
Point(5) = { b/Sqrt(2), b/Sqrt(2), 0 };
Point(6) = { 0, b, 0 };
Point(7) = { d1, 0, 0 };
Point(8) = { a, 0, 0 };
Point(9) = { b, 0, 0 };
Point(10) = { 0, d2, 0 };
Point(11) = { d2/Sqrt(2), d2/Sqrt(2), 0 };
Point(12) = { d2, 0, 0 };
Point(13) = { 0, d3, 0 };
Point(14) = { d3/Sqrt(2), d3/Sqrt(2), 0 };
Point(15) = { d3, 0, 0 };

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

Line(16) = {6, 10};
Line(17) = {5, 11};
Line(18) = {9, 12};

Circle(19) = { 10, 0, 11};
Circle(20) = { 11, 0, 12};

Line(21) = {10, 13};
Line(22) = {11, 14};
Line(23) = {12, 15};

Circle(24) = { 13, 0, 14};
Circle(25) = { 14, 0, 15};

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

Line Loop(6) = {16, 19, -17, -7};
Plane Surface(6) = {6};

Line Loop(7) = {17, 20, -18, -12};
Plane Surface(7) = {7};

Line Loop(8) = {21, 24, -22, -19};
Plane Surface(8) = {8};

Line Loop(9) = {22, 25, -23, -20};
Plane Surface(9) = {9};

Q1[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{:};}};

Physical Surface(100) = Surface{:};

Physical Line(4) = {24, 25, 68, 63};
Physical Line(1) = {21, 16, 5, 1, 14, 47, 27, 32, 52, 62};

//Physical Line(1) = {24, 25, 68, 63, 21, 16, 5, 1, 14, 47, 27, 32, 52, 62};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

Transfinite Line {-21, -22, -23, 64, -62} = 2*r Using Progression 0.95;

