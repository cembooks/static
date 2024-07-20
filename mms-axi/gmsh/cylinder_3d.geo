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

r = 2;//@1

a = 0.4;
b = 1.0;
z = -1.0;

Point(0) = { 0, 0, z };
Point(1) = { 0, a, z };
Point(2) = { a, a, z };
Point(3) = { b/Sqrt(2), b/Sqrt(2), z };
Point(4) = { 0, b, z };
Point(5) = { a, 0, z };
Point(6) = { b, 0, z };

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 5};
Line(4) = {5, 0};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line(5) = {1, 4};
Line(6) = {3, 2};
Circle(7) = { 4, 0, 3};

Line Loop(2) = {-2, 5, 7, 6};
Plane Surface(2) = {2};

Line(8) = {6, 5};
Circle(9) = { 3, 0, 6};

Line Loop(3) = {-3, -6, 9, 8};
Plane Surface(3) = {3};

Q1[] = Rotate { {0,0,1}, {0,0,0},  -Pi/2} { Duplicata{ Surface{1,2,3};} };
Q2[] = Rotate { {0,0,1}, {0,0,0},  -Pi} { Duplicata{ Surface{1,2,3};} };
Q3[] = Rotate { {0,0,1}, {0,0,0},  -3*Pi/2} { Duplicata{ Surface{1,2,3};} };

Extrude {0.0, 0.0, 2.0} {
  Surface{1, 2, 3,
  Q1[0], Q1[1], Q1[2],
  Q2[0], Q2[1], Q2[2],
  Q3[0], Q3[1], Q3[2]};
Layers{2*r}; Recombine;}

Physical Surface(1) = {
1, 2, 3,
Q1[0], Q1[1], Q1[2],
Q2[0], Q2[1], Q2[2],
Q3[0], Q3[1], Q3[2],
75,  97,  119,
273, 295, 317,
207, 229, 251,
141, 163, 185};

Physical Surface(2) = {92, 114, 158, 180, 224, 246, 290, 312};

Physical Volume(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;

