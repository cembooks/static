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

Point(1) = {x0, 0, 0};

Point(2) = {x0, R, 0};
Point(3) = {x0+R/Sqrt(2), R/Sqrt(2), 0};
Point(4) = {x0+R, 0, 0};
Point(5) = {2*x0, 0 , 0};
Point(6) = {2*x0, x0, 0};
Point(7) = {x0, x0, 0};

Point(8) = {2*x0, 2*x0, 0};
Point(9) = {x0, 2*x0, 0};

Point(10) = {Rinfty, 0, 0};
Point(11) = {Rinfty*2/Sqrt(5), Rinfty/Sqrt(5), 0};
Point(12) = {Rinfty/Sqrt(2), Rinfty/Sqrt(2), 0};
Point(13) = {Rinfty/Sqrt(5), Rinfty*2/Sqrt(5), 0};
Point(14) = {0, Rinfty, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 2};
Line(7) = {6, 3};

Line Loop(1) = {2, 3, 4, 7};
Plane Surface(1) = {1};

Line Loop(2) = {-7, 5, 6, 1};
Plane Surface(2) = {2};

Line(8) = {6, 8};
Line(9) = {8, 9};
Line(10) = {9, 7};

Line Loop(3) = {-5, 8, 9, 10};
Plane Surface(3) = {3};

Q1[] = Rotate { {0,0,1}, {x0,0,0},  Pi/2} { Duplicata{ Surface{1, 2};} };
Q2[] = Rotate { {0,0,1}, {x0,0,0},  Pi} { Duplicata{ Surface{1, 2, Q1[0], Q1[1]};} };
Q3[] = Translate {-2*x0, 0, 0} { Duplicata{ Surface{1, 2, Q1[0], Q1[1], Q2[0], Q2[1],
Q2[2], Q2[3] };} };

Q4[] = Translate {-x0, 0, 0} { Duplicata{ Surface{3};} };
Q5[] = Translate {-2*x0, 0, 0} { Duplicata{ Surface{3};} };
Q6[] = Translate {-3*x0, 0, 0} { Duplicata{ Surface{3};} };

Q7[] = Translate {0, -3*x0, 0} { Duplicata{ Surface{3, Q4, Q5, Q6};} };

Line(1001) = {5, 10};
Line(1002) = {6, 11};
Line(1003) = {8, 12};
Line(1004) = {9, 13};
Line(1005) = {199, 14};

Circle(1006) = {10, 37, 11};
Circle(1007) = {11, 37, 12};
Circle(1008) = {12, 37, 13};
Circle(1009) = {13, 37, 14};

Line Loop(1000) = {-4, 1001, 1006, -1002};
Plane Surface(1000) = {1000};

Line Loop(1001) = {-8, 1002, 1007, -1003};
Plane Surface(1001) = {1001};

Line Loop(1002) = {-9, 1003, 1008, -1004};
Plane Surface(1002) = {1002};

Line Loop(1003) = {-84, 1004, 1009, -1005};
Plane Surface(1003) = {1003};

Q8[] = Rotate { {0,0,1}, {0,0,0}, Pi/2} { Duplicata{ Surface{1000, 1001, 1002, 1003};} };
Q9[] = Rotate { {0,0,1}, {0,0,0}, Pi} { Duplicata{ Surface{1000, 1001, 1002,1003, Q8[0] , Q8[1], Q8[2], Q8[3]};} };

Physical Surface(100) = { 1, 2, Q1[0], Q1[1], Q2[0], Q2[1], Q2[2], Q2[3], Q3[0],
Q3[1], Q3[2], Q3[3], Q3[4], Q3[5], Q3[6], Q3[7], 3, Q4, Q5, Q6, Q7[0], Q7[1],
Q7[2], Q7[3], 1000, 1001, 1002, 1003, Q8[0], Q8[1], Q8[2], Q8[3], Q9[0], Q9[1],
Q9[2], Q9[3], Q9[4], Q9[5], Q9[6], Q9[7] };

Physical Line(2) = {1006, 1007, 1008, 1009, 1013, 1018, 1023, 1028, 1033,
1038, 1043, 1048, 1053, 1058, 1063, 1068};
Physical Line(3) = {42, 50, 52, 60, 62, 70, 72, 80};
Physical Line(1) = {2, 1, 12, 20, 22, 30, 32, 40};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = 3*r;

Transfinite Line {-1001, -1002, -1003, -1004, -1005, 1014, 1019, 1024, 1029, 1034,
	1039, 1044, 1049, 1054, 1059, 1064} = (2*m+1)*r Using Progression 0.92;

