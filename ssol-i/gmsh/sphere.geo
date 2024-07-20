r = 5; //@1

d1 = 0.2;
a = 0.5;
b = 1.0;

a3 = a / Sqrt(3);
d13 = d1 / Sqrt(3);
b3 = b / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { d13, d13, d13};
Point(3) = {-d13, d13, d13};

Point(4) = { a3, a3, a3};
Point(5) = {-a3, a3, a3};

Point(6) = { b3, b3, b3};
Point(7) = {-b3, b3, b3};

Line(1) = {2, 4};
Line(2) = {4, 6};

Line(3) = {3, 5};
Line(4) = {5, 7};

Line(5) = {2, 3};

Circle(6) = {4, 1, 5};
Circle(7) = {6, 1, 7};

Curve Loop(1) = {1, 6, -3, -5};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 7, -4, -6};
Plane Surface(2) = {2};

Q1[] = Symmetry {1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};

Curve Loop(4) = {5, 12, -30, 20};
Surface(4) = {4};

Curve Loop(5) = {6, -10, 28, -18};
Surface(5) = {5};

Curve Loop(6) = {7, -15, 33, -23};
Surface(6) = {6};

Surface Loop(1) = {1, 4, 16, 5, 8, 26};
Volume(1) = {1}; 

Surface Loop(2) = {5, 31, 6, 13, 2, 21};
Volume(2) = {2};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1, 2};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, 2, V1[0], V1[1]};}};
V3[] = Symmetry {-1, 0, 1, 0} {Duplicata {Volume{1, 2};}};
V4[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2};}};

Surface Loop(4) = {4, 40, 98, 160, 218, 261};
Volume(4) = {4}; 

Physical Volume(100) = {4};

Physical Volume(10) = {1};
Physical Volume(11) = {2};

Physical Volume(20) = {34};
Physical Volume(21) = {65};

Physical Volume(30) = {92};
Physical Volume(31) = {123};

Physical Volume(40) = {212};
Physical Volume(41) = {243};

Physical Volume(50) = {255};
Physical Volume(51) = {286};

Physical Volume(60) = {154};
Physical Volume(61) = {185};

Physical Surface(2) = {6, 196, 134, 76, 254, 297};

// Rotate the shell to get exact approximation on the x axis. This is 
// needed for projecting the potential on the segment between [a,0,0] 
// and [b,0,0]. 
// V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
// V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;

//+
Show "*";
//+
Show "*";
