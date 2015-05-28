lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 1, 0, lc};
Point(4) = {0.70711, 0.70711, 0, lc};
Point(5) = {2, 0, 0, lc};
Point(6) = {0, 2, 0, lc};
Point(7) = {2, 2, 0, lc};

Line(1) = {3, 6};
Line(2) = {6, 7};
Line(3) = {7, 5};
Line(4) = {5, 2};
Line(5) = {7, 4};

Circle(6) = {3, 1, 4};
Circle(7) = {4, 1, 2};

Line Loop(8) = {2, 5, -6, 1};
Plane Surface(9) = {8};
Line Loop(10) = {3, 4, -7, -5};
Plane Surface(11) = {10};

Transfinite Line {2, 1, 6, 5, 7, 4, 3} = 10 Using Progression 1;

Transfinite Surface {9};
Transfinite Surface {11};

Recombine Surface {9};
Recombine Surface {11};
