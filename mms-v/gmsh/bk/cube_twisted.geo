r = 7; //@1

a = 1.0;

Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
Point(3) = {0, a, 0};
Point(4) = {a, a, 0};
Point(5) = {0, 0, a};
Point(6) = {a, 0, a};
Point(7) = {0, a, a};
Point(8) = {a, a, a};

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};

Line(5) = {3, 4};
Line(6) = {4, 8};
Line(7) = {8, 7};
Line(8) = {7, 3};

Line(9) = {6, 8};
Line(10) = {5, 7};
Line(11) = {1, 3};
Line(12) = {2, 4};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {-6, -5, -8, -7};
Plane Surface(2) = {2};

Line Loop(3) = {7, -10, -3, 9};
Plane Surface(3) = {3};

Line Loop(4) = {-12, -1, 11, 5};
Plane Surface(4) = {4};

Line Loop(5) = {-9, -2, 12, 6};
Plane Surface(5) = {5};

Line Loop(6) = {10, 8, -11, -4};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Q1[] = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Volume{1};} };
Q2[] = Rotate { {0,0,1}, {0,0,0},  Pi} { Duplicata{ Volume{1, Q1[0]};} };
Q3[] = Rotate { {1,0,0}, {0,0,0},  Pi} { Duplicata{ Volume{1, Q1[0], Q2[0],	Q2[1]};} };

Physical Volume(100) = Volume{:};

// Physical Surface(1) = {3, 24, 51, 82, 5, 77, 114, 192, 2, 34, 207, 161};
// Physical Surface(2) = {104, 135, 166, 197, 46, 92, 99, 145, 19, 61, 130, 176};

Physical Surface(2) = {3, 24, 51, 82, 5, 77, 114, 192, 2, 34, 207, 161, 104, 135, 166, 197, 46, 92, 99, 145, 19, 61, 130, 176};

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r;

Mesh 3;

