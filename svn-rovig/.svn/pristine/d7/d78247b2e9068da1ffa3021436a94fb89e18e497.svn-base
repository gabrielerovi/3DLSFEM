function x= gauss_seidel_monolithic(components,M_lev,b_lev,x,C_lev,L_lev,mesh,E_ord,NDtoRT,smoothing_steps)

% For each bc edge/node, we want to assemble the edge/node solution solving a 2x2 system
% we then check if we must impose bc

N_label=mesh.N_label;
E_label=mesh.E_label;
N=mesh.N;
NE=mesh.NE;
E_int=NE-L_lev(1);
N_int=N-(L_lev(3)-L_lev(2));
N_tot=2*NE+2*N;
N_bc=L_lev(4);
contact=0;

for KK=1:smoothing_steps
x_old=x;

% we are dealing with edge elements

if(L_lev(1)>0)
     cont=1;  
    type_of_dof=2;
    % we consider the first point
    ii=1;
    A=[M_lev(ii,ii)   M_lev(ii,ii+1);
       M_lev(ii+1,ii) M_lev(ii+1,ii+1)];
    b(1,1)=b_lev(ii)   -  x_old(ii+2:end)'* M_lev(ii,ii+2:end)' ;
    b(2,1)=b_lev(ii+1) -  x_old(ii+2:end)'* M_lev(ii+1,ii+2:end)' ;

    U_xy=solve_2x2system(A,b);
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(cont), contact);  
    coeff = RT_dirichlet_coeff(E_ord(ii), mesh);
    x(2*ii-1)=coeff*tmp(1);
    x(2*ii)=coeff*tmp(2);    
    
    
for ii=2:L_lev(1)
    cont=cont+1;
    A=[M_lev(ii,ii)   M_lev(ii,ii+1);
       M_lev(ii+1,ii) M_lev(ii+1,ii+1)];
    b(1)=b_lev(ii)   -  x(1:ii-1)'* M_lev(ii,1:ii-1)'   - x_old(ii+2:end)'* M_lev(ii,ii+2:end)' ;
    b(2)=b_lev(ii+1) -  x(1:ii-1)'* M_lev(ii+1,1:ii-1)' - x_old(ii+2:end)'* M_lev(ii+1,ii+2:end)' ;

    U_xy=solve_2x2system(A,b);
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(cont), contact);   
    coeff = RT_dirichlet_coeff(E_ord(ii), mesh);
    x(2*ii-1)=tmp(1)/coeff;
    x(2*ii)=tmp(2)/coeff; 
end


end


if(L_lev(3)>L_lev(2))
 type_of_dof=1;
 cont=0;   
    
for ii=L_lev(2)+1:L_lev(3)
    cont=cont+1;
    
    A=[M_lev(ii,ii)   M_lev(ii,ii+1);
       M_lev(ii+1,ii) M_lev(ii+1,ii+1)];
    b(1,1)=b_lev(ii)   -  x(1:ii-1)'* M_lev(ii,1:ii-1)'   - x_old(ii+2:end)'* M_lev(ii,ii+2:end)' ;
    b(2,1)=b_lev(ii+1) -  x(1:ii-1)'* M_lev(ii+1,1:ii-1)' - x_old(ii+2:end)'* M_lev(ii+1,ii+2:end)' ;

    U_xy=solve_2x2system(A,b);
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(cont), contact);
    x(2*ii-1 - L_lev(2))=tmp(1);
    x(2*ii - L_lev(2) )=tmp(2);
end


end

%  for i=L_lev(4)+1:L_lev(end)
% if(i==480)
% a=1
% end
%  x(i)=1/M_lev(i,i) * (b_lev(i) - x(1:i-1)'*M_lev(i,1:i-1)' -x_old(i+1:end)'*M_lev(i,i+1:end)' );
%  end
 
 
for ii=1:L_lev(4)
   M_lev(ii,:) =zeros(1,N_tot);
   M_lev(ii,ii)=1;
   b_lev(ii)=x(ii);
end

x=M_lev\b_lev';


 res=b_lev'-M_lev*x;

normres(1)=norm(res);
tmp1=1:2:(2*L_lev(1)-1);
tmp2=1:2:(2*E_int-1);
tmp2=tmp2+N_bc;
index_sigma1=cat(2,tmp1,tmp2);
tmp1=2:2:(2*L_lev(1));
tmp2=2:2:(2*E_int);
tmp2=tmp2+N_bc;
index_sigma2=cat(2,tmp1,tmp2);
res1=res(index_sigma1);
res2=res(index_sigma2);

tmp1=L_lev(2)+1:2:(2*(L_lev(3)-L_lev(2) )-1+L_lev(2));
tmp2=L_lev(6)+1:2:(2*(L_lev(7)-L_lev(6))-1+L_lev(6));
index_disp1=cat(2,tmp1,tmp2);
disp1=x(index_disp1);
tmp1=L_lev(2)+2:2:(2*(L_lev(3)-L_lev(2))+L_lev(2));
tmp2=L_lev(6)+2:2:(2*(L_lev(7)-L_lev(6))+L_lev(6));
index_disp2=cat(2,tmp1,tmp2);
disp2=x(index_disp2);
%  
% rhs1=C_lev(1:2:(2*N-1), 2*N+1:3*N ) *disp1 + C_lev( 1:2:(2*N-1), 3*N+1:4*N )*disp2;
% rhs2=C_lev(2:2:(2*N),   2*N+1:3*N ) *disp1 + C_lev( 2:2:(2*N),   3*N+1:4*N )*disp2;
% res_proj1=NDtoRT' *  res1 - rhs1;
% res_proj2=NDtoRT' *  res2 -rhs2;
% res_proj=zeros(2*N,1);
% res_proj(1:2:(2*N-1))=res_proj1;
% res_proj(2:2:(2*N))=res_proj2;
% 
% C=C_lev(1:2*N,1:2*N);
% 
% correction=zeros(2*N,1);
% 
% smoothing_steps=3;
% 
% correction=gauss_seidel_const_bc(components,C,res_proj,correction,smoothing_steps,mesh);
% 
% x(index_sigma1)= x(index_sigma1) + NDtoRT*correction(1:2:(2*N-1));
% x(index_sigma2)  = x(index_sigma2)   + NDtoRT*correction(2:2:2*N);
% 
% res=b_lev'-M_lev*x;
% 
% normres(2)=norm(res);


end