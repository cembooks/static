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

a = 0.4;
d1 = 0.6;
d2 = 0.8;
b = 1.0;

Point(1) = { a , -0.1, 0 };
Point(2) = { d1 , -0.1, 0 };
Point(3) = { d1 ,  0.1, 0 };
Point(4) = { a ,  0.1, 0 };
Point(5) = { d2 ,  0.1, 0 };
Point(6) = { d2 ,  -0.1, 0 };
Point(7) = { b ,  0.1, 0 };
Point(8) = { b ,  -0.1, 0 };

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 6};
Line(6) = {6, 5};
Line(7) = {5, 3};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, -2};
Plane Surface(2) = {2};

Line Loop(3) = {8, 9, 10, -6};
Plane Surface(3) = {3};

Physical Surface(100) = {1,2,3};

Physical Line(1) = {4};
Physical Line(3) = {9};
Physical Line(0) = {1, 3, 5 , 7, 8, 10};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

