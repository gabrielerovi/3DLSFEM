function b = isPointInTriangle(x,y, node)
%ISPOINTINTRIANGLE Test if a point is located inside a triangle
%
%   B = isPointInTriangle(POINT, V1, V2, V3)
%   POINT is a 1-by-2 row vector containing coordinates of the test point,
%   V1, V2 and V3 are 1-by-2 row vectors containing coordinates of triangle
%   vertices. The function returns 1 is the point is inside or on the
%   boundary of the triangle, and 0 otherwise.
%
%   B = isPointInTriangle(POINT, VERTICES)
%   Specifiy the coordinates of vertices as a 3-by-2 array.
%
%   If POINT contains more than one row, the result B has as many rows as
%   the input POINT.
%
%
%   Example
%     % vertices of the triangle
%     p1 = [0 0];
%     p2 = [10 0];
%     p3 = [5 10];
%     tri = [p1;p2;p3];
%     % check if points are inside
%     isPointInTriangle([0 0], tri)
%     ans =
%         1
%     isPointInTriangle([5 5], tri)
%     ans =
%         1
%     isPointInTriangle([10 5], tri)
%     ans =
%         0
%     % check for an array of points
%     isPointInTriangle([0 0;1 0;0 1], tri)
%     ans =
%         1
%         1
%         0
%
%   See also
%   polygons2d, isPointInPolygon, isCounterClockwise
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-05-16,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.
point(1)=x;
point(2)=y;

p1(1)=node(1,1) ; p1(2)= node(1,2);
p2(1)=node(2,1) ; p2(2)= node(2,2);
p3(1)=node(3,1) ; p3(2)= node(3,2);

% if triangle vertices are given as a single array, extract vertices
if nargin == 2
    p2 = p1(2, :);
    p3 = p1(3, :);
    p1 = p1(1, :);
end

% check triangle orientation
isDirect = isCounterClockwise(p1, p2, p3);

big=100000.0;    

% check location of point with respect to each side
AreaT=area(node);
AreaTp(1)=area([p1;p2;point]);
AreaTp(2)=area([p3;p1;point]);
AreaTp(3)=area([p2;p3;point]);

% if isDirect
%     b12 = isCounterClockwise(p1, p2, point) >= -epsilon12;
%     b23 = isCounterClockwise(p2, p3, point) >= -epsilon23;
%     b31 = isCounterClockwise(p3, p1, point) >= -epsilon13;
% else
%     b12 = isCounterClockwise(p1, p2, point) <= epsilon12;
%     b23 = isCounterClockwise(p2, p3, point) <= epsilon23;
%     b31 = isCounterClockwise(p3, p1, point) <= epsilon13;
% end

       

diffarea= abs(AreaT-sum(AreaTp));

% combines the 3 results
%b = b12 & b23 & b31;

epsilon=AreaT/big;
if(diffarea<epsilon)
b=1;
else
    b=0;
end
end
