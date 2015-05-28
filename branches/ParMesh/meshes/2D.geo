cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {0, 1, 0, 1};
Point(4) = {1, 1, 0, 1};
Line(1) = {3, 4};
Line(2) = {4, 2};
Line(3) = {2, 1};
Line(4) = {1, 3};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};


Transfinite Surface {6};
Recombine Surface {6};

Physical Line(2000)={1,2,3,4};
Physical Surface(3300)={6};
