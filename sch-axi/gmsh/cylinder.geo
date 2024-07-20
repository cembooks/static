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

q = 0.0;
b = 1.0;
a = 0.5;

Point(1) = { q , -0.25, 0 };
Point(2) = { a , -0.25, 0 };
Point(3) = { a ,  0.25, 0 };
Point(4) = { q ,  0.25, 0 };
Point(5) = { b ,  0.25, 0 };
Point(6) = { b ,  -0.25, 0 };

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 6};
Line(6) = {6, 5};
Line(7) = {5, 3};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, -2};
Plane Surface(2) = {2};

Physical Surface(100) = {1,2};

Physical Line(1) = {6};
Physical Line(0) = {1, 3, 4, 5, 7};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

