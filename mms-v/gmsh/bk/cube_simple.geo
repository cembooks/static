r = 7; //@1

a = 1.0;

Point(1) = {-a, -a, -a};
Point(2) = { a, -a, -a};
Point(3) = { a,  a, -a};
Point(4) = {-a,  a, -a};
Point(5) = {-a, -a,  a};
Point(6) = { a, -a,  a};
Point(7) = { a,  a,  a};
Point(8) = {-a,  a,  a};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Line Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};

Line Loop(4) = {2, 11, -6, -10};
Plane Surface(4) = {4};

Line Loop(5) = {3, 12, -7, -11};
Plane Surface(5) = {5};

Line Loop(6) = {4, 9, -8, -12};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume(100) = {1};

//Physical Surface(1) = {2, 4, 5};
//Physical Surface(2) = {1, 3, 6};

Physical Surface(2) = {1, 2, 3, 4, 5, 6};

Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = 2*r - 1;

Mesh 3;

