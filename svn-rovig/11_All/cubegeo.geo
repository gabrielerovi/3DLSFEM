lc=1;

Point(1) = {0, 0, 0,lc};
Point(2) = {1, 0, 0,lc};
Point(3) = {1, 1, 0,lc};
Point(4) = {0, 1, 0,lc};
Point(5) = {0, 0, 1,lc};
Point(6) = {1, 0, 1,lc};
Point(7) = {1, 1, 1,lc};
Point(8) = {0, 1, 1,lc};
Line(1) = {8, 7};
Line(2) = {7, 6};
Line(3) = {6, 5};
Line(4) = {5, 8};
Line(5) = {3, 2};
Line(6) = {2, 1};
Line(7) = {1, 4};
Line(8) = {4, 3};
Line(9) = {3, 7};
Line(10) = {2, 6};
Line(11) = {8, 4};
Line(12) = {5, 1};

Line Loop(13) = {8, 5, 6, 7};
Plane Surface(14) = {13};


Line Loop(15) = {1, -9, -8, -11};
Plane Surface(16) = {15};
Line Loop(17) = {9, 2, -10, -5};
Plane Surface(18) = {17};

Line Loop(19) = {3, 12, -6, 10};
Plane Surface(20) = {19};
Line Loop(21) = {12, 7, -11, -4};
Plane Surface(22) = {21};
Line Loop(23) = {2, 3, 4, 1};
Plane Surface(24) = {-23};
Surface Loop(25) = {24, 14, 16, 18, 20, 22};
Volume(999) = {25};


/*
Physical Volume(9999) = {999};

Physical Surface(814) = {14};
Physical Surface(816) = {16};
Physical Surface(818) = {18};
Physical Surface(820) = {20};
Physical Surface(822) = {22};
Physical Surface(824) = {24};
*/