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

d1 = 0.1;
a = 0.9;
b = 1.2;
d2 = 2.0;
d3 = 3.0;

a3 = a / Sqrt(3);
b3 = b / Sqrt(3);
d23 = d2 / Sqrt(3);
d33 = d3 / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { d1, d1, d1};
Point(3) = {-d1, d1, d1};

Point(4) = { a3, a3, a3};
Point(5) = {-a3, a3, a3};

Point(6) = { b3, b3, b3};
Point(7) = {-b3, b3, b3};

Point(8) = { d23, d23, d23};
Point(9) = {-d23, d23, d23};

Point(10) = { d33, d33, d33};
Point(11) = {-d33, d33, d33};

Line(1) = {2, 4};
Line(2) = {4, 6};
Line(3) = {6, 8};

Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 9};

Line(7) = {2, 3};

Circle(8) = {4, 1, 5};
Circle(9) = {6, 1, 7};
Circle(10) = {8, 1, 9};

Line(11) = {9, 11};
Line(12) = {8, 10};

Circle(13) = {10, 1, 11};

Curve Loop(1) = {1, 8, -4, -7};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 9, -5, -8};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 10, -6, -9};
Plane Surface(3) = {3};

Curve Loop(4) = {12, 13, -11, -10};
Plane Surface(4) = {4};

Q1[] = Symmetry {1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};

Curve Loop(5) = {7, 18, -56, 36};
Plane Surface(5) = {5};

Curve Loop(6) = {8, -16, 54, -34};
Surface(6) = {6};

Curve Loop(7) = {9, -21, 59, -39};
Surface(7) = {7};

Curve Loop(8) = {10, -26, 64, -44};
Surface(8) = {8};

Curve Loop(9) = {13, -31, 69, -49};
Surface(9) = {9};

Surface Loop(1) = {1, 14, 52, 32, 5, 6};
Volume(1) = {1};

Surface Loop(2) = {2, 19, 57, 37, 6, 7};
Volume(2) = {2};

Surface Loop(3) = {3, 24, 62, 42, 7, 8};
Volume(3) = {3};

Surface Loop(4) = {4, 29, 67, 47, 8, 9};
Volume(4) = {4};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1, 2, 3, 4};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, 2, 3, 4, V1[0], V1[1], V1[2], V1[3]};}};
V3[] = Symmetry {-1, 0, 1, 0} {Duplicata {Volume{1, 2, 3, 4};}};
V4[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3, 4};}};

Surface Loop(5) = {5, 91, 335, 211, 455, 575};
Volume(5) = {5};

Physical Volume(1) = {Volume{:}};

Physical Surface(1) = {9, 189, 433, 309, 553, 673};

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;
