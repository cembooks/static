r = 5; //@1

d1 = 0.2;
a1 = 0.3;
a2 = 0.5;
b = 1.0;

d13 = d1 / Sqrt(3);
a13 = a1 / Sqrt(3);
a23 = a2 / Sqrt(3);
b3 = b / Sqrt(3);

Point(1) = {0,0,0};

Point(2) = { d13, d13, d13};
Point(3) = {-d13, d13, d13};

Point(4) = { a13, a13, a13};
Point(5) = {-a13, a13, a13};

Point(6) = { a23, a23, a23};
Point(7) = {-a23, a23, a23};

Point(8) = { b3, b3, b3};
Point(9) = {-b3, b3, b3};

Line(1) = {2, 4};
Line(2) = {4, 6};
Line(3) = {6, 8};

Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 9};

Line(7) = {2, 3};

Circle(8) = {4, 1, 5};
Circle(9) = {6, 1, 7};
Circle(10) = {8, 1, 9};

Curve Loop(1) = {1, 8, -4, -7};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 9, -5, -8};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 10, -6, -9};
Plane Surface(3) = {3};

Q1[] = Symmetry {1, 1, 0, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Surface{:};}};

Curve Loop(4) = {7, 15, -43, 28};
Plane Surface(4) = {4};

Curve Loop(5) = {8, -13, 41, -26};
Surface(5) = {5};

Curve Loop(6) = {9, -18, 46, -31};
Surface(6) = {6};

Curve Loop(7) = {10, -23, 51, -36};
Surface(7) = {7};

Surface Loop(1) = {1, 24, 39, 11, 4, 5};
Volume(1) = {1};

Surface Loop(2) = {2, 16, 44, 29, 5, 6};
Volume(2) = {2};

Surface Loop(3) = {3, 21, 49, 34, 6, 7};
Volume(3) = {3};

V1[] = Symmetry {0, 1, -1, 0} {Duplicata {Volume{1, 2, 3};}};
V2[] = Symmetry {0, 1, 1, 0} {Duplicata {Volume{1, 2, 3, V1[0], V1[1], V1[2]};}};
V3[] = Symmetry {-1, 0, 1, 0} {Duplicata {Volume{1, 2, 3};}};
V4[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3};}};

Surface Loop(4) = {73, 433, 255, 344, 4, 162};
Volume(4) = {4}; 

Physical Volume(1) = {Volume{:}};

Physical Surface(1) = {7, 140, 322, 229, 411, 500};

// Rotate the shell to get exact approximation on the x axis. This is 
// needed for projecting the potential on the segment between [a,0,0] 
// and [b,0,0]. 
// V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
// V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

// Finally, we mesh the results ----------------------------------------------
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r; 

