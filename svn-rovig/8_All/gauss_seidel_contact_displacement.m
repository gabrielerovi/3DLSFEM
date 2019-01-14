function[x]=gauss_seidel_contact_displacement(A,b,x,smoothing_steps,mesh)

L=size(mesh);
L=L(1);
n=length(b);
nhalf=0.5*n;
type_of_dof=1;
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);
N_bc=mesh{L}.N_bc;


if(n>1)
for kk=1:smoothing_steps
x_old=x;
ii=1;
jj=ii+nhalf; 
   
A_loc=[A(ii,ii) A(ii,jj);
       A(jj,ii) A(jj,jj);];
b_loc=[A(ii,ii) - x_old(ii+1:nhalf)'*A(ii,ii+1:nhalf)' -  x_old(jj+1:n)'*A(jj,jj+1:n)';
       A(jj,jj) - x_old(ii+1:nhalf)'*A(ii,ii+1:nhalf)'-   x_old(jj+1:n)'*A(jj,jj+1:n)'; ]
x([ii,jj])=A_loc\b_loc; 
if(N_bc(ii)>0 && abs(dirichlet(N_bc(ii),3)-1)<0.00001)
    dof=ii;
    disp=x([ii,jj]);        
   [x([ii;jj]),Nconstraint]=ContactConstraint(mesh{L},disp,dof,type_of_dof,parameters);
end

for ii=2:nhalf
jj=ii+nhalf;   
A_loc=[A(ii,ii) A(ii,jj);
       A(jj,ii) A(jj,jj);];
b_loc=[A(ii,ii)-x(1:ii-1)'*A(ii,1:ii-1)' - x(nhalf+1:jj-1)'*A(jj,nhalf+1:jj-1)' -x_old(ii+1:nhalf)'*A(ii,ii+1:nhalf)' -x_old(jj+1:n)'*A(jj,jj+1:n)';
       A(jj,jj)-x(1:ii-1)'*A(ii,1:ii-1)' - x(nhalf+1:jj-1)'*A(jj,nhalf+1:jj-1)' -x_old(ii+1:nhalf)'*A(ii,ii+1:nhalf)' -x_old(jj+1:n)'*A(jj,jj+1:n)'; ]
disp=(A_loc\b_loc); 
x(ii)=disp(1); x(jj)=disp(2); 
if(N_bc(ii)>0 && abs(dirichlet(N_bc(ii),3)-1)<0.00001)
    dof=ii;
    disp=x([ii,jj]);        
   [x([ii,jj]),Nconstraint]=ContactConstraint(mesh{L},disp,dof,type_of_dof,parameters);
end

end

end
else
    x(1)=b(1)/A(1,1);
end

end
    
    