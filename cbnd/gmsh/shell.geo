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

a3 = a / Sqrt(3);
b3 = b / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { a3, a3, a3};
Point(3) = {-a3, a3, a3};

Point(4) = { b3, b3, b3};
Point(5) = {-b3, b3, b3};

Line(1) = {2, 4};
Line(2) = {5, 3};

Circle(6) = {3, 1, 2};
Circle(7) = {4, 1, 5};

Curve Loop(1) = {1, 7, 2, 6};
Plane Surface(1) = {1};

S2 = Rotate { {1,0,0}, {0,0,0}, -Pi/2} { Duplicata{ Surface{1};} };

Circle(13) = {2, 1, 6};
Circle(14) = {7, 1, 4};

Circle(15) = {16, 1, 3};
Circle(16) = {5, 1, 12};

Curve Loop(2) = {9, 14, -1, 13};
Plane Surface(2) = {2};

Curve Loop(3) = {-2, 16, 11, 15};
Plane Surface(3) = {3};

Curve Loop(4) = {-6, -15, 12, -13};
Surface(4) = {4};

Curve Loop(5) = {-7, -14, 10, -16};
Surface(5) = {5};

Surface Loop(1) = {1, 2, 8, 3, 4, 5};
Volume(1) = {1};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, V1[0]};}};
V3[] = Symmetry {1, 1, 0, 0} {Duplicata {Volume{1};}};
V4[] = Symmetry {-1, 1, 0, 0} {Duplicata {Volume{1};}};

Physical Volume(1) = {1};
Physical Volume(2) = {75};
Physical Volume(3) = {44};
Physical Volume(4) = {17};
Physical Volume(5) = {129};
Physical Volume(6) = {102};

Physical Surface(1) = {4, 96, 65, 38, 150, 123};
Physical Surface(3) = {5, 101, 70, 43, 155, 128};

// Rotate the shell to get exact approximation on the x axis. This is
// needed for projecting the potential on the segment between [a,0,0]
// and [b,0,0].
V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{1, 75, 44, 17, 129, 102}; };
V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{1, 75, 44, 17, 129, 102}; };

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;

