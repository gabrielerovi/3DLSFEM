function  ContactDisplacement(mesh,x,parameters)




L=size(mesh);
L=L(1);

N=mesh{L}.N;
NE=mesh{L}.NE;
NT=mesh{L}.NT;
E_bc=mesh{L}.E_bc;

u1=x(1+2*NE:2*NE+N);
u2=x(1+N+2*NE:end);

type_of_dof=2;
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);

kk=[3;1;2];
%initialize contact pressure


cont=0;
for nn=mesh{L}.N_contact
    cont=cont+1;
       u_x(cont)=u1(nn);
       u_y(cont)=u2(nn);
       coor(cont)=mesh{L}.node(nn,1);

end





figure
plot(coor,u_x,'ro','MarkerSize',20,'MarkerFaceColor','red');
hold on
plot(coor,u_y,'bo','MarkerSize',20,'MarkerFaceColor','blue');
l=legend('u_x','u_y');
l.FontSize=46;
xx=xlabel('X');
yy=ylabel('U');
set(gca,'fontsize',48)

end