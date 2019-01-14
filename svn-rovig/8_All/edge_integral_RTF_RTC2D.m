
function z =edge_integral_RTF_RTC2D(k_coarse,segment,nodeC,N,normal,p0C)
% This function evaluates \iint_K f(x,y) dxdy using
% the Gaussian quadrature of order N where K is a
% triangle with vertices (x1,y1), (x2,y2) and (x3,y3).
xw = MonoGaussPoints(N); % get quadrature points and weights
% find number of Gauss points 
NP=length(xw(:,1));

x1=segment(1,1) ; y1= segment(1,2);
x2=segment(2,1) ; y2= segment(2,2);

% calculate the length of the segment 

L=sqrt( (x1-x2)^2 + (y1-y2)^2 );
big_number=100000.0;
eps=L/big_number;

z =0.0;
%if x1==x2, then the segment is vertical
% therefore we have a line integral in 1D (y instead of x)

for j = 1:NP
q_point(j,1)= x1 +  0.5 * (xw(j,1)+1) * (x2-x1) ;
q_point(j,2)= y1 +  0.5 * (xw(j,1)+1) * (y2-y1) ;
bool(j)=fix(isPointInTriangle(q_point(j,1),q_point(j,2), nodeC));
end


[RT0_basisC,RT0_divergenceC] = phiRT2D(q_point,nodeC);
                

    for j=1:NP
     f(j)=bool(j) * ( RT0_basisC(k_coarse,j,1)*normal(1)+RT0_basisC(k_coarse,j,2)*normal(2) );
    end
    
    % 0.5 L because the reference element is [-1,1]
    z=0.5 * L.* f* xw(:,2);



end
