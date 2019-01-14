function [A,b1,b2] =ComplementarityCondition(A,L,b1,b2,mesh,parameters)

E_bc=mesh.E_bc;
NE=mesh.NE;
gap=parameters.gap;
C_contact=parameters.C_contact;
% we are only interested in contact nodes, edges. 
% so we look at the 3th column of dirichlet (1 = contact, 0 = no contact)
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(2);

N=mesh.N;
NE=mesh.NE;


    
% for ee=1:NE
%     if(E_bc(ee)>0)
%         if(dirichlet(E_bc(ee),3)==1)
%            edge=mesh.edge(ee,:);
%            for jj=1:2
%            for kk=1:4
%                
%            A{L,jj,kk}(ee,:)=0;
%            for ii=1:2
%            if(mesh.N_dirichlet(edge(ii))==0)
%            A{L,jj+2,kk}(edge(ii),:)=0;
%            end
%            end
%            
%            end
%            end
%            
%            
%         end
%         
%     end
%     
% end





for ee=1:NE
    % E_bc(ee)>0 : I am no the boundary
    if(E_bc(ee)>0)
        % dirichlet(E_bc(ee),3)==1: the boundary belongs to GammaC
        if(dirichlet(E_bc(ee),3)==1)
            edge=mesh.edge(ee,:);
            node=mesh.node(edge,:);
            
            % the normals we are going to use are the ones of the obstacle
            normal_edge=mesh.normal_edge_contact{ee};
            normal_node(1,:)=mesh.normal_node_contact{edge(1)};
            normal_node(2,:)=mesh.normal_node_contact{edge(2)};
            
            phidotn=phi_dot_n(mesh,ee);
            
            Lside=norm( node(1,:)-node(2,:));
            
            % n_obstacle (sigma n_body) =constant, u n_obstacle = linear
            % therefore we can compute exact integrals
            
            % <phidton(S1 n_obst1 + S2 n_obst2 ), u1 n_obst1+ u2 n_obst2- g>
            % A13=  <phidton S1 n_obst1, u1 n_obst1>
            % A14= <phidton S1 n_obst1 , u2 n_obst2>
            % A23= <phidton S2 n_obst2, u1 n_obst1>
            % A24= <phidton  S2 n_obst2 ,  u2 n_obst2>
            % F1 = <phidton S1 n_obst1,g>
            % F2 =  <phidton S2 n_obst2 , g>
            
            % WE ASSUME THAT INT_[A,B] = INT_[A+EPS,B-EPS], SO WE REMOVE
            % DISCONTINITUI ON THE BORDER
            
            % SO HERE NORMAL NODE  AND NORMAL EDGE SHOULD BE THE SAME,
            % MAYBE
            for kk=1:2
%             v11(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(1) * normal_node(kk,1);
%             v12(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(1) * normal_node(kk,2);
%             v21(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(2) * normal_node(kk,1);
%             v22(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(2) * normal_node(kk,2);
            v11(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(1) * normal_edge(1);
            v12(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(1) * normal_edge(2);
            v21(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(2) * normal_edge(1);
            v22(kk) = C_contact * 0.5 * Lside * phidotn * normal_edge(2) * normal_edge(2);
            
            end
            
  
            % g =whatever, so we need to compute accurate integrals
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
            rhsintegrand(kk,1)=gap_function(xx(kk),yy(kk)) ;  
            end
    
            
            F1  = C_contact * Lside * phidotn * normal_edge(1) * xw(:,2)'* rhsintegrand;
            F2  = C_contact * Lside * phidotn * normal_edge(2) * xw(:,2)'* rhsintegrand;
            
            for kk=1:2
            A{L,1,3}(ee,edge(kk))= A{L,1,3}(ee,edge(kk)) + v11(kk);
            A{L,1,4}(ee,edge(kk))= A{L,1,4}(ee,edge(kk)) + v12(kk);
            A{L,2,3}(ee,edge(kk))= A{L,2,3}(ee,edge(kk)) + v21(kk);
            A{L,2,4}(ee,edge(kk))= A{L,2,4}(ee,edge(kk)) + v22(kk);
            
            if(mesh.N_dirichlet(edge(kk))==0)
            A{L,3,1}(edge(kk),ee)= A{L,3,1}(edge(kk),ee) + v11(kk);
            A{L,3,2}(edge(kk),ee)= A{L,3,2}(edge(kk),ee) + v21(kk);
            A{L,4,1}(edge(kk),ee)= A{L,4,1}(edge(kk),ee) + v12(kk);
            A{L,4,2}(edge(kk),ee)= A{L,4,2}(edge(kk),ee) + v22(kk);
            end        
           
            end
            b1(ee)=b1(ee)+F1;
            b2(ee)=b2(ee)+F2;
           
            
        end
    
        
    end
    
end








end