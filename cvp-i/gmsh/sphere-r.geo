r = 5; //@1

d1 = 0.1;
a1 = 0.3;
a2 = 0.5;
b = 1.0;

d13 = d1 / Sqrt(3);
a13 = a1 / Sqrt(3);
a23 = a2 / Sqrt(3);
b3 = b / Sqrt(3);

// First, we construct the top segment of the shell -------------------------- 
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

S2 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{1,2,3};} };
S3 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{S2[0],S2[1],S2[2]};} };
S4 = Rotate { {0,0,1}, {0,0,0},  Pi/2} { Duplicata{ Surface{S3[0],S3[1],S3[2]};} };

Curve Loop(4) = {-7, 15, 29 , 43};
Plane Surface(4) = {4};

Curve Loop(5) = {8, 13, 27, 41};
Surface(5) = {5};

Curve Loop(6) = {9, 18, 32, 46};
Surface(6) = {6};

Curve Loop(7) = {10, 23, 37, 51};
Surface(7) = {7}; 

Surface Loop(1) = {1, 11, 25, 39, 4, 5};
Volume(1) = {1};

Surface Loop(2) = {2, 16, 30, 44, 5, 6};
Volume(2) = {2};

Surface Loop(3) = {3, 21, 35, 49, 6, 7};
Volume(3) = {3};


// Next, we replicate the top segment five times by rotation -----------------
V2 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{1, 2, 3};} };
V3 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{V2[0],V2[1],V2[2]};} };
V4 = Rotate { {0,1,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{V3[0],V3[1],V3[2]};} };
V5 = Rotate { {1,0,0}, {0,0,0},  Pi/2} { Duplicata{ Volume{1, 2, 3};} };
V6 = Rotate { {1,0,0}, {0,0,0}, -Pi/2} { Duplicata{ Volume{1, 2, 3};} };

Surface Loop(4) = {73, 429, 251, 340, 4, 162};
Volume(4) = {4}; // The top segment of the shell

Physical Volume(1) = {Volume{:}};

// Assign boundary IDs -----------------------------------------------------
Physical Surface(1) = {7, 140, 229, 318, 407, 496};

// Rotate the shell to get exact approximation on the x axis. This is 
// needed for projecting the potential on the segment between [a,0,0] 
// and [b,0,0]. 
//V0 = Rotate { {0,0,1}, {0,0,0}, -Pi/4} { Volume{:}; };
//V0 = Rotate { {0,1,0}, {0,0,0}, 0.195913276015303635085*Pi} { Volume{:}; };

// Finally, we mesh the results ----------------------------------------------
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";

Transfinite Line "*" = r; 

//Mesh 3;

