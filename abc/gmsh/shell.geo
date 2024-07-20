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

r = 1; //@1
m = 1; //@2

R = 5.000000e-02;
x0 = 1.000000e-01;

Rmid = 1.5 * (R + x0);

Rinfty = 2 * m * Rmid;

r_in3 = x0 / Sqrt(3);
r_out3 = Rmid / Sqrt(3);
r_inf3 = Rinfty / Sqrt(3);

// First, we construct the top segment of the shell --------------------------
Point(1) = {0,0,0};

Point(2) = { r_in3, r_in3, r_in3};
Point(3) = {-r_in3, r_in3, r_in3};
Point(4) = { r_out3, r_out3, r_out3};
Point(5) = {-r_out3, r_out3, r_out3};

Point(6) = { r_inf3, r_inf3, r_inf3};
Point(7) = {-r_inf3, r_inf3, r_inf3};

Line(1) = {2, 4};
Line(2) = {3, 5};

Line(3) = {4, 6};
Line(4) = {5, 7};

Circle(5) = {2, 1, 3};
Circle(6) = {4, 1, 5};
Circle(7) = {6, 1, 7};

Curve Loop(1) = {-5, 1, 6, -2};
Plane Surface(1) = {1};

Curve Loop(2) = {-6, 3, 7, -4};
Plane Surface(2) = {2};

S2 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{1, 2};} };
S3 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{S2[0], S2[1]};} };
S4 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{S3[0], S3[1]};} };

Curve Loop(3) = {-5, 9, 19, 29};
Surface(3) = {3};

Curve Loop(4) = {6, 11, 21, 31};
Surface(4) = {4};

Curve Loop(5) = {7, 16, 26, 36};
Surface(5) = {5};

Surface Loop(1) = {1, 8, 18, 28, 3, 4};
Volume(1) = {1}; // The top segment of the shell

Surface Loop(2) = {2, 13, 23, 33, 4, 5};
Volume(2) = {2}; // The top segment of the shell

// Next, we replicate the top segment five times by rotation -----------------
V2 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{1, 2};} };
V3 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{V2[0], V2[1]};} };
V4 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{V3[0], V3[1]};} };

V5 = Rotate { {1,0,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{1, 2};} };
V6 = Rotate { {1,0,0}, {0,0,0}, -Pi/2} { Duplicata{ Volume{1, 2};} };

// Assign an index to a volume ----------------------------------------------
Physical Volume(101) = {37};
Physical Volume(102) = {269};
Physical Volume(103) = {153};
Physical Volume(104) = {211};
Physical Volume(105) = {1};
Physical Volume(106) = {95};

Physical Volume(107) = {68};
Physical Volume(108) = {300};
Physical Volume(109) = {184};
Physical Volume(110) = {242};
Physical Volume(111) = {2};
Physical Volume(112) = {126};

// Assign boundary IDs ------------------------------------------------------
Physical Surface(1) = {3, 232, 116, 290, 174, 58};
Physical Surface(2) = {5, 268, 152, 326, 94, 210};

//Delete Embedded "*";

// Finally, we mesh the results ----------------------------------------------
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = 3*r;

Transfinite Line {-3, -4, 17, 27, -71, 83, 141, -129} = (2*m+1)*r Using Progression 0.92;

//Mesh 3;

