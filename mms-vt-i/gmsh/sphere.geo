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

d1 = 0.5;
a = 1.0;

a3 = a / Sqrt(3);
d13 = d1 / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { d13, d13, d13};
Point(3) = {-d13, d13, d13};

Point(4) = { a3, a3, a3};
Point(5) = {-a3, a3, a3};

Line(1) = {2, 4};
Line(2) = {3, 5};
Line(3) = {2, 3};

Circle(4) = {4, 1, 5};

Curve Loop(1) = {1, -2, -3, 4};
Plane Surface(1) = {1};

Q1[] = Symmetry {1, 1, 0, 0} {Duplicata {Surface{1};}};
Q2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{1, Q1[0]};}};

Curve Loop(4) = {3, 9, -19, 14};
Surface(4) = {4};

Curve Loop(6) = {4, -7, 17, -12};
Surface(6) = {6};

Surface Loop(1) = {1, 5, 10, 15, 4, 6};
Volume(1) = {1};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, V1[0]};}};
V3[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1};}};
V4[] = Symmetry {-1, 0, 1, 0} {Duplicata {Volume{1};}};

Surface Loop(4) = {4, 41, 99, 68, 126, 153};
Volume(4) = {4};

Physical Volume(0) = {4};
Physical Volume(1) = {1};
Physical Volume(2) = {20};
Physical Volume(3) = {47};
Physical Volume(4) = {78};
Physical Volume(5) = {105};
Physical Volume(6) = {132};

Physical Surface(1) = {6, 46, 158};
Physical Surface(2) = {104, 131, 73};

// Physical Surface(2) = {6, 46, 158, 104, 131, 73};

// Rotate the shell to get exact approximation on the x axis. This is
// needed for projecting the potential on the segment between [a,0,0]
// and [b,0,0].
// V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
// V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = 2*r - 1;

