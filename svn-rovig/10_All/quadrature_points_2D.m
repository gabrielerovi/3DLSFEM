
function [q_point,weights,area]=quadrature_points_2D(N,node)
% This function evaluates \iint_K f(x,y) dxdy using
% the Gaussian quadrature of order N where K is a
% triangle with vertices (x1,y1), (x2,y2) and (x3,y3).
xw = TriGaussPoints(N); % get quadrature points and weights
% find number of Gauss points 
NP=length(xw(:,1));

res = isCounterClockwise(node(1,:), node(2,:), node(3,:));

if(res==1)
x1=node(1,1) ; y1= node(1,2);
x2=node(2,1) ; y2= node(2,2);
x3=node(3,1) ; y3= node(3,2);
else
x1=node(1,1) ; y1= node(1,2);
x3=node(2,1) ; y3= node(2,2);
x2=node(3,1) ; y2= node(3,2); 
end

% calculate the area of the triangle 

area=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0;

for j = 1:NP
q_point(j,1) = x1*(1-xw(j,1)-xw(j,2))+x2*xw(j,1)+x3*xw(j,2) ;
q_point(j,2) = y1*(1-xw(j,1)-xw(j,2))+y2*xw(j,1)+y3*xw(j,2) ;
end
weights=xw(:,3);
end
