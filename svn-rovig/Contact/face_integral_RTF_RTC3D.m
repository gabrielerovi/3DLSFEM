
function z =face_integral_RTF_RTC3D(kC,face,nodeC,N,normalF,p0C)
% This function evaluates \iint_F f(x,y,z) dxdydz using
% the Gaussian quadrature of order N where K is a
% tetrahedron with vertices (x1,y1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4).
xw = TriGaussPoints(N); % get quadrature points and weights
% find number of Gauss points 
NP=length(xw(:,1));
dim=3;
face_per_elem=4;

% calculate the length of the segment 

Area=AreaTriangle(face);


% nodes_triangleABC=[0 0 0; 1 0 0 ; 0 1 0];
% 
% for j = 1:NP
% point=xw(j,1:2);
% point(3)=0;
% alpha(j,:)= BarycentricCoordinatesTriangle(point, nodes_triangleABC);
% q_point(j,:)=sparse(1,3);
% for kk=1:3
% q_point(j,:)=q_point(j,:)+alpha(j,kk) * face(kk,:) ;
% end
% end


J=[face(2,1)-face(1,1), face(3,1)-face(1,1);
   face(2,2)-face(1,2), face(3,2)-face(1,2);
   face(2,3)-face(1,3), face(3,3)-face(1,3);];

bJ=[face(1,1);face(1,2);face(1,3)];

for j=1:NP
  qpointref=[xw(j,1);xw(j,2)];
  q_point(j,:)=bJ+J* qpointref; 
end

[RT0_basisC] = phiRT3Dcell(q_point,nodeC);
                

    
    
    for qp=1:NP
          f(qp)= ( RT0_basisC{kC}(qp,:)*normalF );
    end
    
    % 0.5 L because the reference element is [-1,1]
    z=Area * f* xw(:,end);



end
