
r = 5; //@1

d1 = 0.1;
a1 = 0.3;
b1 = 0.6;
a2 = 0.9;
b2 = 1.2;
d2 = 2.4;

a13 = a1/Sqrt(3);
b13 = b1/Sqrt(3);

a23 = a2/Sqrt(3);
b23 = b2/Sqrt(3);

Point(0) = { 0, 0, 0};

Point(1) = { d1, d1, d1};
Point(2) = { d1,-d1, d1};

Point(3) = { a13, a13, a13};
Point(4) = { a13,-a13, a13};

Point(5) = { b13, b13, b13};
Point(6) = { b13,-b13, b13};

Point(7) = { a23, a23, a23};
Point(8) = { a23,-a23, a23};

Point(9) = { b23, b23, b23};
Point(10) = { b23,-b23, b23};

Point(11) = { d2, d2, d2};
Point(12) = { d2,-d2, d2};

Line(1) = {1, 3};
Line(2) = {3, 5};
Line(3) = {5, 7};
Line(4) = {7, 9};
Line(5) = {9, 11};

Line(7) = {2, 4};
Line(8) = {4, 6};
Line(9) = {6, 8};
Line(10) = {8, 10};
Line(11) = {10, 12};

Line(13) = {1, 2};

Circle(14) = {3, 0, 4};
Circle(15) = {5, 0, 6};
Circle(16) = {7, 0, 8};
Circle(17) = {9, 0, 10};

Circle(18) = {11, 0, 12};

Line Loop(1) = {7, -14, -1, 13};
Plane Surface(1) = {1};

Line Loop(2) = {8, -15, -2, 14};
Plane Surface(2) = {2};

Line Loop(3) = {9, -16, -3, 15};
Plane Surface(3) = {3};

Line Loop(4) = {10, -17, -4, 16};
Plane Surface(4) = {4};

Line Loop(5) = {11, -18, -5, 17};
Plane Surface(5) = {5};

Q1[] = Symmetry {0, 1, 1, 0} {Duplicata {Surface{:};}};
Q2[] = Symmetry {0, -1, 1, 0} {Duplicata {Surface{:};}};

ll = newll; Line Loop(ll) = {13, -23, 72, -47};
Plane Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {14, 21, -70, 45};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {15, 26, -75, 50};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {16, 31, -80, 55};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {17, 36, -85, 60};
Surface(ll) = {ll};

ll = newll; Line Loop(ll) = {18, 41, -90, 65};
Surface(ll) = {ll};

Surface Loop(1) = {1, 43, 68, 19, 91, 92};
Volume(1) = {1};

Surface Loop(2) = {2, 48, 73, 24, 92, 93};
Volume(2) = {2};

Surface Loop(3) = {3, 53, 78, 29, 93, 94};
Volume(3) = {3};

Surface Loop(4) = {4, 58, 83, 34, 94, 95};
Volume(4) = {4};

Surface Loop(5) = {5, 63, 88, 39, 95, 96};
Volume(5) = {5};

V1[] = Symmetry {1, 1, 0, 0} {Duplicata {Volume{:};}};
V2[] = Symmetry {-1, 1, 0, 0} {Duplicata {Volume{:};}};
V3[] = Symmetry {1, 0, 1, 0} {Duplicata {Volume{1, 2, 3, 4, 5};}};
V4[] = Symmetry {1, 0, -1, 0} {Duplicata {Volume{1, 2, 3, 4, 5};}};

Surface Loop(6) = {91, 118, 424, 269, 575, 726};
Volume(6) = {6};

Physical Volume(1) = {Volume{:}};

Physical Surface(1) = {96, 247, 556, 401, 710, 861};

Recombine Surface "*";
Recombine Volume "*";

Transfinite Volume "*";
Transfinite Surface "*";
Transfinite Line "*" = r;

