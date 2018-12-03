lc=100;
ref=1;
ref1=ref;
ref2=ref;
ref3=ref;


Point(1) = {0, 0, 0,lc};
Point(2) = {ref1, 0, 0,lc};
Point(3) = {0, ref2, 0,lc};
Point(4) = {0, 0, ref3,lc};

/*
Point(1) = {0, 1, 3,lc};
Point(2) = {2, 2, 4,lc};
Point(3) = {3, 4, 9,lc};
Point(4) = {4, 11, 12,lc};
*/


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Line(4) = {2,4};
Line(5) = {3, 4};
Line(6) = {4, 1};


Line Loop(7) = {1,4,6};
Plane Surface(8) = {7};

Line Loop(9) = {2,5,-4};
Plane Surface(10) = {9};

Line Loop(11) = {3,-6,-5};
Plane Surface(12) = {11};

Line Loop(13) = {1,2,3};
Plane Surface(14) = {13};

Surface Loop(15) = {8,10,12,14};
Volume(999) = {15};


/*
Physical Volume(9999) = {999};

Physical Surface(814) = {14};
Physical Surface(816) = {16};
Physical Surface(818) = {18};
Physical Surface(820) = {20};
Physical Surface(822) = {22};
Physical Surface(824) = {24};
*/
