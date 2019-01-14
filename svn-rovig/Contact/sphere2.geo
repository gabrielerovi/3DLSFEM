hfine = 0.0025;
hfinetop = 0.05;
hcoarsetop=0.1;
hcoarsemiddle=0.15;
radius=0.5;
zz=radius*0.9;
phi=Acos(zz/radius);

Point(1) = {radius, 0, radius, hcoarsetop};
Point(2) = {-radius, 0, radius, hcoarsetop};
Point(3) = {0, radius, radius, hcoarsetop};
Point(4) = {0, -radius, radius, hcoarsetop};

Point(5) = {radius*Sin(phi), 0, radius-zz, hcoarsemiddle};
Point(6) = {-radius*Sin(phi), 0, radius-zz, hcoarsemiddle};
Point(7) = {0, radius*Sin(phi), radius-zz, hcoarsemiddle};
Point(8) = {0, -radius*Sin(phi), radius-zz, hcoarsemiddle};

Point(9) = {0, 0, radius, hfinetop};
Point(10) = {0, 0, 0, hfine};
Point(11) = {0, 0, radius-zz, hfine};


Circle(1) = {1, 9, 3};
Circle(2) = {3, 9, 2};
Circle(3) = {2, 9, 4};
Circle(4) = {4, 9, 1};


Circle(5) = {1, 9, 5};
Circle(6) = {2, 9, 6};
Circle(7) = {3, 9, 7};
Circle(8) = {4, 9, 8};

Circle(9) = {5, 11, 7};
Circle(10) = {7, 11, 6};
Circle(11) = {6, 11, 8};
Circle(12) = {8, 11, 5};

Circle(13) = {5, 9, 10};
Circle(14) = {6, 9, 10};
Circle(15) = {7, 9, 10};
Circle(16) = {8, 9, 10};




Line(17) = {1,9};
Line(18) = {3,9};
Line(19) = {2,9};
Line(20) = {4,9};


Line Loop(21) = {1,18,-17};
Ruled Surface(21) = {21};

Line Loop(22) = {2,19,-18};
Ruled Surface(22) = {22};

Line Loop(23) = {3,20,-19};
Ruled Surface(23) = {23};

Line Loop(24) = {4,17,-20};
Ruled Surface(24) = {24};

Physical Surface(700)={21,22,23,24};





Line Loop(25) = {1,7,-9,-5};
Ruled Surface(25) = {25};

Line Loop(26) = {2,6,-10,-7};
Ruled Surface(26) = {26};

Line Loop(27) = {3,8,-11,-6};
Ruled Surface(27) = {27};

Line Loop(28) = {4,5,-12,-8};
Ruled Surface(28) = {28};

Physical Surface(800)={25,26,27,28};






Line Loop(29) = {9,15,-13};
Ruled Surface(29) = {29};

Line Loop(30) = {10,14,-15};
Ruled Surface(30) = {30};

Line Loop(31) = {11,16,-14};
Ruled Surface(31) = {31};

Line Loop(32) = {12,13,-16};
Ruled Surface(32) = {32};

Physical Surface(900)={29,30,31,32};

Surface Loop(333) = {21,22,23,24,25,26,27,28,29,30,31,32};

Volume(333) = {333};

Physical Volume(0)={333};


/*
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

Physical Surface(800)={29,30,31,32};
Physical Volume(0)={33};
*/
