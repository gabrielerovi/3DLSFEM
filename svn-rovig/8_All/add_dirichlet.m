function y =add_dirichlet(mesh,x,type_of_dof,component)
cont=0;
L=size(mesh);
L=L(1);
% dof=node
if (type_of_dof==1)
    [yy,w]=dirichlet_bc_vector(mesh{L}.N_bc,type_of_dof,mesh);
for n=1:mesh{L}.N
    if(mesh{L}.N_dirichlet(n)>0)
        y(n)=yy(n,component);
        %y(n)=dirichlet_bc_vector(mesh{L}.N_bc(n),type_of_dof);
    else
        cont=cont+1;
        y(n)=x(cont);
    end
end
% dof=edge
elseif(type_of_dof==2)
    
    [yy,w]=dirichlet_bc_vector(mesh{L}.E_bc,type_of_dof,mesh{L});
 for e=1:mesh{L}.NE
    if(mesh{L}.E_dirichlet(e)>0)
        y(e)=yy(e,component);
        %y(e)=dirichlet_bc_vector(mesh{L}.E_bc(e),type_of_dof,mesh{L});
    else
        cont=cont+1;
        y(e)=x(cont);
    end
 end   
 
elseif(type_of_dof==3)
 for f=1:mesh{L}.NF
    if(mesh{L}.F_dirichlet(f)>0)
        [yy,w]=dirichlet_bc_vector(mesh{L}.F_bc(f),type_of_dof);
        y(f)=yy(component);
        %y(f)=dirichlet_bc_vector(mesh{L}.F_bc(f),type_of_dof);
    else
        cont=cont+1;
        y(f)=x(cont);
    end
 end   

end
x=y;

end