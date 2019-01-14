
function [q_point,weights,Volume]=quadrature_points_3D(N,node)
% This function evaluates \iint_K f(x,y) dxdy using
% the Gaussian quadrature of order N where K is a
% triangle with vertices (x1,y1), (x2,y2) and (x3,y3).
xw = TriGaussPoints3D(N); % get quadrature points and weights
% find number of Gauss points 
NP=length(xw(:,1));

x=node(:,1);
y=node(:,2);
z=node(:,3);

% calculate the area of the triangle 
Volume=VolumeTetrahedronAndNormalsigns(node); 

for j = 1:NP
q_point(j,1) = x(1) * (1-xw(j,1)-xw(j,2)-xw(j,3)-xw(j,4)) + x(2) * xw(j,1) + x(3) * xw(j,2) + x(4) * xw(j,3);
q_point(j,2) = y(1) * (1-xw(j,1)-xw(j,2)-xw(j,3)-xw(j,4)) + y(2) * xw(j,1) + y(3) * xw(j,2) + y(4) * xw(j,3);
q_point(j,3) = z(1) * (1-xw(j,1)-xw(j,2)-xw(j,3)-xw(j,4)) + z(2) * xw(j,1) + z(3) * xw(j,2) + z(4) * xw(j,3);
end
% q_point(:,1) = x(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + x(2) * xw(:,1) + x(3) * xw(:,2) + x(4) * xw(:,3);
% q_point(:,2) = y(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + y(2) * xw(:,1) + y(3) * xw(:,2) + y(4) * xw(:,3);
% q_point(:,3) = z(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + z(2) * xw(:,1) + z(3) * xw(:,2) + z(4) * xw(:,3);


weights=xw(:,end);
end