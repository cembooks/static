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

r = 3; //@1

a = 0.4;
b = 1.0;
d_1 = 0.6;
d_2 = 0.8;

a2 = a / Sqrt(2);
b2 = b / Sqrt(2);
d2_1 = d_1 / Sqrt(2);
d2_2 = d_2 / Sqrt(2);

a3 = a / Sqrt(3);
b3 = b / Sqrt(3);
d3_1 = d_1 / Sqrt(3);
d3_2 = d_2 / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { a3, a3, a3};
Point(3) = { d3_1, d3_1, d3_1};
Point(4) = { d3_2, d3_2, d3_2};
Point(5) = { b3, b3, b3};

Point(6) = { 0, a, 0};
Point(7) = { 0, d_1, 0};
Point(8) = { 0, d_2, 0};
Point(9) = { 0, b, 0};

Point(10) = { 0, a2, a2};
Point(11) = { 0, d2_1, d2_1};
Point(12) = { 0, d2_2, d2_2};
Point(13) = { 0, b2, b2};

Point(14) = { a2, a2, 0};
Point(15) = { d2_1, d2_1, 0};
Point(16) = { d2_2, d2_2, 0};
Point(17) = { b2, b2, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};

Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 9};

Line(7) = {10, 11};
Line(8) = {11, 12};
Line(9) = {12, 13};

Line(10) = {14, 15};
Line(11) = {15, 16};
Line(12) = {16, 17};

Circle(13) = {2, 1, 14};
Circle(14) = {14, 1, 6};
Circle(15) = {6, 1, 10};
Circle(16) = {10, 1, 2};

Circle(17) = {3, 1, 15};
Circle(18) = {15, 1, 7};
Circle(19) = {7, 1, 11};
Circle(20) = {11, 1, 3};

Circle(21) = {4, 1, 16};
Circle(22) = {16, 1, 8};
Circle(23) = {8, 1, 12};
Circle(24) = {12, 1, 4};

Circle(25) = {5, 1, 17};
Circle(26) = {17, 1, 9};
Circle(27) = {9, 1, 13};
Circle(28) = {13, 1, 5};

Curve Loop(1) = {13, 14, 15, 16};
Surface(1) = {1};

Curve Loop(2) = {17, 18, 19, 20};
Surface(2) = {2};

Curve Loop(3) = {21, 22, 23, 24};
Surface(3) = {3};

Curve Loop(4) = {25, 26, 27, 28};
Surface(4) = {4};


Curve Loop(5) = {16, 1, -20, -7};
Plane Surface(5) = {5};

Curve Loop(6) = {20, 2, -24, -8};
Plane Surface(6) = {6};

Curve Loop(7) = {24, 3, -28, -9};
Plane Surface(7) = {7};


Curve Loop(8) = {13, 10, -17, -1};
Plane Surface(8) = {8};

Curve Loop(9) = {17, 11, -21, -2};
Plane Surface(9) = {9};

Curve Loop(10) = {21, 12, -25, -3};
Plane Surface(10) = {10};


Curve Loop(11) = {14, 4, -18, -10};
Plane Surface(11) = {11};

Curve Loop(12) = {18, 5, -22, -11};
Plane Surface(12) = {12};

Curve Loop(13) = {22, 6, -26, -12};
Plane Surface(13) = {13};


Curve Loop(14) = {15, 7, -19, -4};
Plane Surface(14) = {14};

Curve Loop(15) = {19, 8, -23, -5};
Plane Surface(15) = {15};

Curve Loop(16) = {23, 9, -27, -6};
Plane Surface(16) = {16};

Surface Loop(1) = {5, 8, 11, 14, 1, 2};
Volume(1) = {1};

Surface Loop(2) = {6, 9, 12, 15, 2, 3};
Volume(2) = {2};

Surface Loop(3) = {7, 10, 13, 16, 3, 4};
Volume(3) = {3};

V1 = Rotate { {0,1,0}, {0,0,0}, -Pi/2} { Duplicata{Volume{:};} };
V2 = Rotate { {0,1,0}, {0,0,0}, -Pi} { Duplicata{ Volume{:};} };
V3 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Volume{:};} };
V4 = Rotate { {0,0,1}, {0,0,0},  Pi} { Duplicata{ Volume{:};} };
V5 = Rotate { {1,0,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{1, 2, 3, V1(0), V1(1), V1(2), V2(0), V2(1), V2(2), V2(3), V2(4), V2(5)};} };
V6 = Rotate { {1,0,0}, {0,0,0}, -Pi/2} { Duplicata{ Volume{1, 2, 3, V1(0), V1(1), V1(2), V2(0), V2(1), V2(2), V2(3), V2(4), V2(5)};} };

Physical Volume(1) = {Volume{:}} ;

// The physical surfaces are numbered in the software.

// Rotate the shell to get exact approximation on the x axis. This is
// needed for projecting the potential on the segment between [a,0,0]
// and [b,0,0].
V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

Recombine Surface "*";
Recombine Volume "*";

Transfinite Surface "*";
Transfinite Volume "*";

Transfinite Line "*" = r;

