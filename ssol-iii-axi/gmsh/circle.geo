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

a1 = 0.3;
b1 = 0.6;
a2 = 0.9;
b2 = 1.2;

d1 = 0.1;
d2 = 1.5; 
d3 = 3.0;

Point(0) = { 0, 0, 0 };
Point(1) = { 0, d1, 0 };
Point(2) = { 0, a1, 0 };
Point(3) = { a1/Sqrt(2), a1/Sqrt(2), 0 };
Point(4) = { d1, d1, 0 };
Point(5) = { b1/Sqrt(2), b1/Sqrt(2), 0 };
Point(6) = { 0, b1, 0 };
Point(7) = { d1, 0, 0 };
Point(8) = { a1, 0, 0 };
Point(9) = { b1, 0, 0 };

Point(10) = { 0, a2, 0 };
Point(11) = { 0, b2, 0 };
Point(12) = { 0, d2, 0 };

Point(13) = { a2/Sqrt(2), a2/Sqrt(2), 0 };
Point(14) = { b2/Sqrt(2), b2/Sqrt(2), 0 };
Point(15) = { d2/Sqrt(2), d2/Sqrt(2), 0 };

Point(16) = { a2, 0, 0 };
Point(17) = { b2, 0, 0 };
Point(18) = { d2, 0, 0 };

Point(19) = { 0, d3, 0 };
Point(20) = { d3/Sqrt(2), d3/Sqrt(2), 0  };
Point(21) = { d3, 0, 0 };

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
Line(17) = {10, 11};
Line(18) = {11, 12};

Line(19) = {5, 13};
Line(20) = {13, 14};
Line(21) = {14, 15};

Line(22) = {9, 16};
Line(23) = {16, 17};
Line(24) = {17, 18};

Circle(25) = { 10, 0, 13};
Circle(26) = { 11, 0, 14};
Circle(27) = { 12, 0, 15};

Circle(28) = { 13, 0, 16};
Circle(29) = { 14, 0, 17};
Circle(30) = { 15, 0, 18};

Line(31) = {12, 19};
Line(32) = {15, 20};
Line(33) = {18, 21};

Circle(34) = { 19, 0, 20};
Circle(35) = { 20, 0, 21};

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

Line Loop(6) = {16, 25, -19, -7};
Plane Surface(6) = {6};

Line Loop(7) = {17, 26, -20, -25};
Plane Surface(7) = {7};

Line Loop(8) = {18, 27, -21, -26};
Plane Surface(8) = {8};

Line Loop(9) = {19, 28, -22, -12};
Plane Surface(9) = {9};

Line Loop(10) = {20, 29, -23, -28};
Plane Surface(10) = {10};

Line Loop(11) = {21, 30, -24, -29};
Plane Surface(11) = {11};

Line Loop(12) = {31, 34, -32, -27};
Plane Surface(12) = {12};

Line Loop(13) = {32, 35, -33, -30};
Plane Surface(13) = {13};

Q1[] = Symmetry { 0, 1, 0, 0} {Duplicata {Surface{:};}};

//Physical Surface(100) = Surface{:};

Physical Surface(1) = {5, 56, 1, 3, 46, 36, 6, 9, 76, 61, 8, 11, 86, 71, 12, 13,
96, 91};

Physical Surface(2) = {2, 4, 51, 41};
Physical Surface(15) = {7, 10, 81, 66};


Physical Line(4) = {34, 35, 98, 93};
Physical Line(1) = {31, 18, 17, 16, 5, 1, 14, 57, 37, 42, 62, 67, 72, 92};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

Transfinite Line {-31, -32, -33, 94, -92} = 2*r Using Progression 0.95;

