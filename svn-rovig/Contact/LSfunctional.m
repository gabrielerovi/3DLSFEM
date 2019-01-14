function [J]=LSfunctional(x,f1,f2,alpha,beta, coeff_equilibrium,coeff_constitutive,coeff_asymmetry,mesh,qrule)

L=size(mesh);
L=L(1);


N=mesh{L}.N;
NE=mesh{L}.NE;

stress1=x( 1:NE );
stress2=x( 1 + NE : 2 * NE);
disp1=  x( 1 + 2 * NE : 2 * NE + N);
disp2=  x(1 + 2 * NE + N : end);

% initialize the functional J
J=zeros(3,1);


for t=1:mesh{L}.NT
    % dof vertices of the triangle
    elem=mesh{L}.elem(t,:);
    % dofs edge of the triangle
    elemE=mesh{L}.elemE(t,:);
    %node coordinates
    node=mesh{L}.node(elem,:);
    
    % recover local unknowns
    stress1_loc=stress1(elemE);
    stress2_loc=stress2(elemE);
    disp1_loc=disp1(elem);
    disp2_loc=disp2(elem);
    
    % compute quadrature points and weights
    [q_point,weights,area]=quadrature_points_2D(qrule,node);
    number_of_qp=length(q_point(:,1));
    % compute gradient of displacement
    grad = P1grad2D  (node);
    %compute basis and divergence of RT functions
    [RT0_basis,RT0_divergence]=phiRT2Dcell(q_point,node);
    
    
    for qp=1:number_of_qp
        for dd=1:2
        sigma1(qp,dd)=RT0_basis{1}(qp,dd)* stress1_loc(1)+RT0_basis{2}(qp,dd)* stress1_loc(2)+RT0_basis{3}(qp,dd)* stress1_loc(3);
        sigma2(qp,dd)=RT0_basis{1}(qp,dd)* stress2_loc(1)+RT0_basis{2}(qp,dd)* stress2_loc(2)+RT0_basis{3}(qp,dd)* stress2_loc(3);        
        end
        epsxx(qp)= (grad(1,1)*disp1_loc(1)+grad(2,1)*disp1_loc(2)+grad(3,1)*disp1_loc(3));
        epsyy(qp)= (grad(1,2)*disp2_loc(1)+grad(2,2)*disp2_loc(2)+grad(3,2)*disp2_loc(3));
        epsxy(qp)=0.5 * ( epsxx(qp) + epsyy(qp) );
        divsigma1(qp)= sum(RT0_divergence .* stress1_loc);
        divsigma2(qp)= sum(RT0_divergence .* stress2_loc);
    end

    
         
         
for qp=1:number_of_qp
  value(qp)= weights(qp)* ...
  (coeff_equilibrium * ( (divsigma1(qp) + f1(q_point(qp,1),q_point(qp,2)) )^2 + ...
                         (divsigma2(qp) + f2(q_point(qp,1),q_point(qp,2)) )^2));% + ...   
end
J(1)=J(1)+sum(value)*area;

for qp=1:number_of_qp
  value(qp)= weights(qp)* ...                     
   (coeff_constitutive * ( +( beta * sigma1(qp,1)+alpha * (sigma1(qp,1)+sigma2(qp,2)) - epsxx(qp))^2 ...
                           +( beta * sigma2(qp,2)+alpha * (sigma1(qp,1)+sigma2(qp,2)) - epsyy(qp))^2 ...
                           +( beta * sigma1(qp,2)                                     - epsxy(qp))^2 ...
                           +( beta * sigma2(qp,1)                                     - epsxy(qp))^2 ));

J(2)=J(2)+sum(value)*area;


for qp=1:number_of_qp
  value(qp)= weights(qp)* ...                     
   (coeff_asymmetry * ( 2*( sigma1(qp,2) - sigma2(qp,1))^2 ));%
end
 J(3)=J(3)+sum(value)*area;                       
end


end


