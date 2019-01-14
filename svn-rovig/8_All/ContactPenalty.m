function [A,b1,b2] =ContactPenalty(A,L,b1,b2,mesh,parameters)

E_bc=mesh.E_bc;
NE=mesh.NE;
gap=parameters.gap;
C_contact=parameters.C_contact;
% we are only interested in contact nodes, edges. 
% so we look at the 3th column of dirichlet (1 = contact, 0 = no contact)
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(2);

N=mesh.N;
NE=mesh.NE;

for ee=1:NE
    % E_bc(ee)>0 : I am no the boundary
    if(E_bc(ee)>0)
        % dirichlet(E_bc(ee),3)==1: the boundary belongs to GammaC
        if(dirichlet(E_bc(ee),3)==1)
            edge=mesh.edge(ee,:);
            node=mesh.node(edge,:);
            normal_edge=mesh.normal_edge{ee};
            normal_node(1,:)=mesh.normal_node{edge(1)};
            normal_node(2,:)=mesh.normal_node{edge(2)};
            phidotn=phi_dot_n(mesh,ee);
            Lside=norm( node(1,:)-node(2,:));
            
            
            % int ( un )^2, where we have 4 dofs for each edge, 2 for the
            % first node, 2 for the second
            % the shape functions are (1-x/L) and (x/L)
            %therefore if we have u and v on the same node (same shape
            %function), the integral is:
            % int  (1-x/L)^2 = int (x/L)^2 = L/3
            % Otherwise:
            % int  (1-x/L) x /L =  L/6
            
            for kk=1:2
               A{L,3,3}(edge(kk),edge(kk))  =  A{L,3,3}(edge(kk),edge(kk)) + C_contact * (Lside/3.0) * normal_node(kk,1) * normal_node(kk,1);
               A{L,3,4}(edge(kk),edge(kk))  =  A{L,3,4}(edge(kk),edge(kk)) + C_contact * (Lside/3.0) * normal_node(kk,1) * normal_node(kk,2);
               A{L,4,3}(edge(kk),edge(kk))  =  A{L,4,3}(edge(kk),edge(kk)) + C_contact * (Lside/3.0) * normal_node(kk,2) * normal_node(kk,1);
               A{L,4,4}(edge(kk),edge(kk))  =  A{L,4,4}(edge(kk),edge(kk)) + C_contact * (Lside/3.0) * normal_node(kk,2) * normal_node(kk,2);
            end
            
            for k=1:2
                if(k==1)
                    ii=1;
                    jj=2;
                else
                    ii=2;
                    jj=1;
                end
            A{L,3,3}(edge(ii),edge(jj))  =  A{L,3,3}(edge(ii),edge(jj)) +  C_contact * (Lside/6.0) * normal_node(ii,1) * normal_node(jj,1);
            A{L,3,4}(edge(ii),edge(jj))  =  A{L,3,4}(edge(ii),edge(jj)) +  C_contact * (Lside/6.0) * normal_node(ii,1) * normal_node(jj,2);
            A{L,4,3}(edge(ii),edge(jj))  =  A{L,4,3}(edge(ii),edge(jj)) +  C_contact * (Lside/6.0) * normal_node(ii,2) * normal_node(jj,1);
            A{L,4,4}(edge(ii),edge(jj))  =  A{L,4,4}(edge(ii),edge(jj)) +  C_contact * (Lside/6.0) * normal_node(ii,2) * normal_node(jj,2); 
            end

            qrule=4;
            xw = MonoGaussPoints(qrule);
            node1=node(1,[1,2]);
            node2=node(2,[1,2]);
            for kk=1:length(xw(:,1))
            xx(kk)= 0.5 * ( node2(1) + node1(1) ) + xw(kk,1) * ( node2(1) - node1(1) );
            yy(kk)= 0.5 * ( node2(2) + node1(2) ) + xw(kk,1) * ( node2(2) - node1(2) );
            convective_vector=[xx(kk)-node1(1),yy(kk)-node1(2)];
            s=norm(convective_vector);
            % g(1-x/L)  (here x is the distance from node1, of the actual (x,y)
            rhsintegrand1(kk,1)= gap(xx(kk),yy(kk)) * (1-s/Lside );  
            % g * x/L  
            rhsintegrand2(kk,1)= gap(xx(kk),yy(kk)) *  s/Lside ;     
            end
            
            % we consider the two nodes, edge(1), edge(2)
            % for fixed node i=1,2 and edge, the basis function is fixed
            % Then depending on the component (b1,b2), we have 
            % normal_x or normal_y
            
            b1(edge(1)) = b1(edge(1)) + C_contact * Lside * xw(:,2)'* rhsintegrand1 * normal_node(1,1);
            b2(edge(1)) = b2(edge(1)) + C_contact * Lside * xw(:,2)'* rhsintegrand1 * normal_node(2,1);
            b1(edge(2)) = b1(edge(2)) + C_contact * Lside * xw(:,2)'* rhsintegrand2 * normal_node(1,1);
            b2(edge(2)) = b2(edge(2)) + C_contact * Lside * xw(:,2)'* rhsintegrand2 * normal_node(2,1);
            
            
            
                        % int (sigma n n)^2, where sigma n n = constant function
            % we have only two dofs for each edge (1 and 2th raws of stress
            % tensor)
            A{L,1,1}(ee,ee) = A{L,1,1}(ee,ee) +  C_contact * Lside * phidotn * normal_edge(1) * phidotn * normal_edge(1);
            A{L,1,2}(ee,ee) = A{L,1,2}(ee,ee) +  C_contact * Lside * phidotn * normal_edge(1) * phidotn * normal_edge(2);
            A{L,2,1}(ee,ee) = A{L,2,1}(ee,ee) +  C_contact * Lside * phidotn * normal_edge(2) * phidotn * normal_edge(1);
            A{L,2,2}(ee,ee) = A{L,2,2}(ee,ee) +  C_contact * Lside * phidotn * normal_edge(2) * phidotn * normal_edge(2);
            
        end
  
         
           
            
        end
    
        
    end
    
end








