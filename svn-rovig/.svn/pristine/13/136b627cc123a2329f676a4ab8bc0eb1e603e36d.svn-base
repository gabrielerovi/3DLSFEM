hfine = 0.025;
hcoarse=0.12;
radius=0.5;
Point(1) = {0, 0, radius, hcoarse};
Point(2) = {radius, 0, radius, hcoarse};
Point(3) = {-radius, 0, radius, hcoarse};
Point(4) = {0, radius, radius, hcoarse};
Point(5) = {0, -radius, radius, hcoarse};
Point(7) = {0, 0, 0, hfine};
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};
Circle(7) = {3, 1, 7};
Circle(8) = {7, 1, 2};
Circle(11) = {5, 1, 7};
Circle(12) = {7, 1, 4};

Line(13) = {1,2};
Line(14) = {4,1};
Line(15) = {3,1};
Line(16) = {5,1};


Line Loop(14) = {2, 7, 12};
Ruled Surface(14) = {14};
Line Loop(20) = {3, 11, -7};
Ruled Surface(20) = {20};
Line Loop(22) = {4, -8, -11};
Ruled Surface(22) = {22};
Line Loop(28) = {1, -12, 8};
Ruled Surface(28) = {28};

Line Loop(29) = {13,1,14};
Plane Surface(29) = {29};

Line Loop(30) = {2,15,-14};
Plane Surface(30) = {30};

Line Loop(31) = {15,-16,-3};
Plane Surface(31) = {31};

Line Loop(32) = {13,-4,16};
Plane Surface(32) = {32};

Surface Loop(33) = {14, 20, 22, 28,29,30,31,32};


Volume(33) = {33};
Physical Surface(700)={14,20,22,28};
Physical Surface(800)={29,30,31,32};
Physical Volume(0)={33};
