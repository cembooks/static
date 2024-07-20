r = 5; //@1
m = 1;

d1 = 0.2;
a = 0.5;
b = m * 1.0;

gridsize = 250e-03;

Point(0) = { 0, 0, 0, gridsize};
Point(1) = { d1, 0 , 0, gridsize};
Point(2) = { d1, d1, 0, gridsize};
Point(3) = { 0, d1, 0, gridsize};

Point(4) = {a, 0, 0, gridsize};
Point(5) = {a/Sqrt(2), a/Sqrt(2), 0, gridsize};
Point(6) = { 0, a , 0, gridsize};

Point(7) = { b, 0, 0, gridsize};
Point(8) = {b/Sqrt(2), b/Sqrt(2), 0, gridsize};
Point(9) = { 0, b, 0, gridsize};

Line(1) = {1, 4};
Line(2) = {4, 7};

Line(3) = {2, 5};
Line(4) = {5, 8};

Line(5) = {3, 6};
Line(6) = {6, 9};

Line(7) = {1, 2};
Line(8) = {2, 3};

Circle(9) = {4, 0, 5};
Circle(10) = {5, 0, 6};

Circle(11) = {7, 0, 8};
Circle(12) = {8, 0, 9};

Line Loop(1) = {1, 9, -3, -7};
Plane Surface(1) = {1};

Line Loop(2) = {2, 11, -4, -9};
Plane Surface(2) = {2};

Line Loop(3) = {3, 10, -5, -8};
Plane Surface(3) = {3};

Line Loop(4) = {4, 12, -6, -10};
Plane Surface(4) = {4};

Q1[] = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{1, 2, 3, 4};} };

Q2[] = Rotate { {0,0,1}, {0,0,0}, Pi }
{ Duplicata{ Surface{1,2,3,4,Q1[0],Q1[1],Q1[2],Q1[3]};} };

l = newl;  Line(l) = {0, 1};
l = newl;  Line(l) = {0, 3};
l = newl;  Line(l) = {0, 55};
l = newl;  Line(l) = {0, 114};

ll = newll; Line Loop(ll) = {7, 8, -71, 70};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {71, -17, -27, -72};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-73, 72, -36, -46};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-66, -70, 73, -56};
s = news; Plane Surface(s) = {ll};

Physical Surface(100) = { 1, 2, 3, 4, 75, 77, 79, 81, Q1[0], Q1[1], Q1[2], Q1[3],
Q2[0], Q2[1], Q2[2], Q2[3], Q2[4], Q2[5], Q2[6],Q2[7]};

Physical Line(1) = {11, 12, 20, 30, 39, 49, 59, 69};

//Physical Line(1) = {11, 12, 20, 30};
//Physical Line(2) = {39, 49, 59, 69};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

