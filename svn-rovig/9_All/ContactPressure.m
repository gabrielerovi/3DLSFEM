function p = ContactPressure(mesh,x,parameters)




L=size(mesh);
L=L(1);

N=mesh{L}.N;
NE=mesh{L}.NE;
NT=mesh{L}.NT;
E_bc=mesh{L}.E_bc;
node_tot=mesh{L}.node;
edge_tot=mesh{L}.edge;
sigma1=x(1:NE);
sigma2=x(1+NE:2*NE);

type_of_dof=2;
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);

kk=[3;1;2];
%initialize contact pressure


p = zeros(NE,1);
cont=0;
for ee=mesh{L}.E_contact

       % the edge ee belongs to the element t
       t=cell2mat(mesh{L}.E_to_T{ee});
       % nodes dof of the element
       elem=mesh{L}.elem(t,:);
       % edges of the element
       elemE=mesh{L}.elemE(t,:);
       % nodes coordinate of the element
       node=mesh{L}.node(elem,:);  
       % find the local edge and its nodes
       opposite_node=kk(find(elemE==ee));
       edge=setdiff([1,2,3],opposite_node);
       % midpoint of the edge on contact boundary
       midpoint = 0.5 * (node(edge(1),[1,2])+node(edge(2),[1,2]));
       % shape functions of the element t defined in the midpoint
       phi=phiRT2D(midpoint,node(:,[1,2]));       
       % define local stress tensor in the midpoint
       LocalStress=zeros(2,2);
       
       for et=1:3
           % we compute the stress in the midpoint of the edge, considering
           % the sum of the shape functions on the whole element t, multiplied with the respective
           % coefficients. Indeed, since we are not multiplying by the
           % normal, each edge play a role in the definition of the stress.
           LocalStress=LocalStress+full(  [ phi(et,1,1) * sigma1(elemE(et)), phi(et,1,2) * sigma1(elemE(et));
                                           phi(et,1,1) * sigma2(elemE(et)), phi(et,1,2) * sigma2(elemE(et));]);
       end    
           LocalStress=0.5 * (LocalStress+LocalStress');
           principalstresses=eig(LocalStress);
           p(ee) = mean( principalstresses ) ;
           
           
           
           
           
           % here we compute the pressure as the normal component of:
           % sigma n

            normal=normal_contact(mean(node_tot(edge_tot(ee,:),:)),parameters);
            %mesh{L}.normal_edge_contact{ee};                 

            phi=phiRT2D(midpoint,node);
            
            phidotn=phi_dot_n(mesh{L},ee);
            p(ee) = phidotn * ( normal(1) * sigma1(ee) +  normal(2) * sigma2(ee) );
            cont=cont+1;
            coord(cont,[1,2])=midpoint;
            pressure(cont,1)=p(ee);
end

pressure=-pressure;
figure
plot(coord(:,1),pressure,'ro','MarkerSize',20,'MarkerFaceColor','red');
l.FontSize=46;
xx=xlabel('X');
yy=ylabel('P');
set(gca,'fontsize',48)
ylim([0,max(pressure)+0.005]);
xlim([-0.1, 0.1]);

end