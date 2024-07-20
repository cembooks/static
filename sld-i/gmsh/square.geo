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

r = 10; //@1

d1 = 0.07;
a = 0.2;
b = 0.4;
d2 = 0.8;
d3 = 2.0;

Point(0) = { 0, 0, 0};

Point(1) = { d1, 0 , 0};
Point(2) = { d1, d1, 0};
Point(3) = { 0, d1, 0};

Point(4) = {a, 0, 0};
Point(5) = {a/Sqrt(2), a/Sqrt(2), 0};
Point(6) = { 0, a, 0};

Point(7) = { b, 0, 0};
Point(8) = { b/Sqrt(2), b/Sqrt(2), 0};
Point(9) = { 0, b, 0};

Point(10) = { d2, 0, 0};
Point(11) = { d2, d2, 0};
Point(12) = { 0, d2, 0};

Point(13) = { d3, 0, 0};
Point(14) = { d3, d3, 0};
Point(15) = { 0, d3, 0};

Line(1) = {1, 4};
Line(2) = {4, 7};
Line(3) = {7, 10};

Line(4) = {2, 5};
Line(5) = {5, 8};
Line(6) = {8, 11};

Line(7) = {3, 6};
Line(8) = {6, 9};
Line(9) = {9, 12};

Line(10) = {10, 13};
Line(11) = {11, 14};
Line(12) = {12, 15};

Line(13) = {1, 2};
Line(14) = {2, 3};

Circle(15) = {4, 0, 5};
Circle(16) = {5, 0, 6};

Circle(17) = {7, 0, 8};
Circle(18) = {8, 0, 9};

Line(19) = {10, 11};
Line(20) = {11, 12};

Line(21) = {13, 14};
Line(22) = {14, 15};

Line Loop(1) = {1, 15, -4, -13};
Plane Surface(1) = {1};

Line Loop(2) = {2, 17, -5, -15};
Plane Surface(2) = {2};

Line Loop(3) = {3, 19, -6, -17};
Plane Surface(3) = {3};

Line Loop(4) = {10, 21, -11, -19};
Plane Surface(4) = {4};

Line Loop(5) = {4, 16, -7, -14};
Plane Surface(5) = {5};

Line Loop(6) = {5, 18, -8, -16};
Plane Surface(6) = {6};

Line Loop(7) = {6, 20, -9, -18};
Plane Surface(7) = {7};

Line Loop(8) = {11, 22, -12, -20};
Plane Surface(8) = {8};

Q1[] = Symmetry {1, 0, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {0, 1, 0, 0} {Duplicata {Surface{:};}};

l = newl;  Line(l) = {0, 1};
l = newl;  Line(l) = {0, 3};
l = newl;  Line(l) = {0, 16};
l = newl;  Line(l) = {0, 152};

ll = newll; Line Loop(ll) = {13, 14, -140, 139};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {140, 47, 27, -141};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-142, 141, -105, -125};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {65, -139, 142, 85};
s = news; Plane Surface(s) = {ll};

Physical Surface(100) = Surface{:};

Physical Line(1) = {21, 22, 60, 40, 118, 138, 98, 78};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

Transfinite Line {-10, -11, -12, 41, -39, 119, 99, 79} = r Using Progression 0.85;

