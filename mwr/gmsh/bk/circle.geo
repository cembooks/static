r = 5; //@1

d1 = 0.1;
a1 = 0.3;
a2 = 0.7;
b = 1.0;

Point(0) = { 0, 0, 0};
Point(1) = { d1, 0 , 0};
Point(2) = { d1, d1, 0};
Point(3) = { 0, d1, 0};

Point(4) = {a1, 0, 0};
Point(5) = {a1/Sqrt(2), a1/Sqrt(2), 0};
Point(6) = { 0, a1 , 0};

Point(7) = {a2, 0, 0};
Point(8) = {a2/Sqrt(2), a2/Sqrt(2), 0};
Point(9) = { 0, a2, 0};

Point(10) = {b, 0, 0};
Point(11) = {b/Sqrt(2), b/Sqrt(2), 0};
Point(12) = { 0, b, 0};

Line(1) = {1, 4};
Line(2) = {4, 7};
Line(3) = {7, 10};

Line(4) = {2, 5};
Line(5) = {5, 8};
Line(6) = {8, 11};

Line(7) = {3, 6};
Line(8) = {6, 9};
Line(9) = {9, 12};

Line(10) = {1, 2};
Line(11) = {2, 3};

Circle(12) = {4, 0, 5};
Circle(13) = {5, 0, 6};

Circle(14) = {7, 0, 8};
Circle(15) = {8, 0, 9};

Circle(16) = {10, 0, 11};
Circle(17) = {11, 0, 12};

Line Loop(1) = {1, 12, -4, -10};
Plane Surface(1) = {1};

Line Loop(2) = {2, 14, -5, -12};
Plane Surface(2) = {2};

Line Loop(3) = {3, 16, -6, -14};
Plane Surface(3) = {3};

Line Loop(4) = {4, 13, -7, -11};
Plane Surface(4) = {4};

Line Loop(5) = {5, 15, -8, -13};
Plane Surface(5) = {5};

Line Loop(6) = {6, 17, -9, -15};
Plane Surface(6) = {6};

Q1[] = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{1, 2, 3, 4, 5, 6};} };
Q2[] = Rotate { {0,0,1}, {0,0,0}, Pi}
{ Duplicata{ Surface{1,2,3,4,5,6,Q1[0],Q1[1],Q1[2],Q1[3],Q1[4],Q1[5]};} };

l = newl;  Line(l) = {0, 1};
l = newl;  Line(l) = {0, 3};
l = newl;  Line(l) = {0, 76};
l = newl;  Line(l) = {0, 171};

ll = newll; Line Loop(ll) = {10, 11, -106, 105};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {106, -22, -37, -107};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-108, 107, -51, -66};
s = news; Plane Surface(s) = {ll};

ll = newll; Line Loop(ll) = {-96, -105, 108, -81};
s = news; Plane Surface(s) = {ll};

Physical Surface(1) = {Surface{:}} ;

Physical Line(1) = {16, 17, 30, 45, 59, 74, 89, 104};

Recombine Surface "*";

Transfinite Surface "*";
Transfinite Line "*" = r;

