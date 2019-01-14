 /*Point(1) = {-100, 100, 0, 1e+22};
  Point(2) = {100, 100, 0, 1e+22};
  Point(3) = {100, -100, 0, 1e+22};
  Point(4) = {-100, -100, 0, 1e+22};
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};
    Line Loop(6) = {4, 1, 2, 3};
  Plane Surface(6) = {6};
  Physical Surface("bottom") = {1};
   Physical Surface("right") = {2}; 
   Physical Surface("top") = {3};  
  Physical Surface("left") = {4};*/

  // Inputs
  squareSide = 1; //m
  meshThickness = squareSide / 1; 
  gridsize = squareSide / 4;
 
        // Geometry
  Point(1) = {0,0, 0, gridsize};
  Point(2) = {squareSide,0, 0, gridsize};
  Point(3) = {squareSide, squareSide, 0, gridsize};
  Point(4) = {0, squareSide, 0, gridsize};
  Line(1) = {1, 2};       // bottom line
  Line(2) = {2, 3};       // right line
  Line(3) = {3, 4};       // top line
  Line(4) = {4, 1};       // left line
  Line Loop(5) = {1, 2, 3, 4};  
  Plane Surface(6) = {5};
 
        //Transfinite surface:
  Transfinite Surface {6};
  //Recombine Surface {6};

