Point(1) = {0,0,0,0.1,10};
Point(2) = {1,0,0,0.1,10};
Point(3) = {1,1,0,0.1,10};
Point(4) = {0,1,0,0.1,10};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Extrude  {0,0,1}
{Surface{6};
}

Surface Loop(29) = {15,6,19,23,27,28};
Volume(30) = {29};
Transfinite Line {8,9,10,11,22,13,2,1,3,14,4,18} = 21;
Progression=1;
Transfinite Surface {28} = {1,10,6,4};
Transfinite Surface {27} = {14,2,1,10};
Transfinite Surface {19} = {14,2,3,5};
Transfinite Surface {15} = {3,5,6,4};
Transfinite Surface {23} = {6,5,14,10};
Transfinite Surface {6} = {2,3,4,1};
Recombine Surface{28,19,15,6,27,23};
Transfinite Volume{1} = {10,14,5,6,1,2,3,4};