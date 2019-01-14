function [M_lev,b_lev, C_lev,L_lev,go_back,N_ord,E_ord,RTCtoRTF,P1CtoP1F,NDCtoNDF,NDtoRT]=create_levels2D ...
          (A11,A12,A13,A14,A22,A23,A24,A33,A34,A44,C11,C12,C13,C14,C22,C23,C24,b1,b2,b3,b4,mesh,RTCtoRTF,P1CtoP1F,NDCtoNDF,NDtoRT)

L=size(mesh);

L=L(1);
M_lev=cell(L,1);
C_lev=cell(L,1);
A11_lev=cell(L,1);
A12_lev=cell(L,1);
A13_lev=cell(L,1);
A14_lev=cell(L,1);

A22_lev=cell(L,1);
A23_lev=cell(L,1);
A24_lev=cell(L,1);

A33_lev=cell(L,1);
A34_lev=cell(L,1);
A44_lev=cell(L,1);

C11_lev=cell(L,1);
C12_lev=cell(L,1);
C22_lev=cell(L,1);

A11_lev{L}=A11;
A12_lev{L}=A12;
A13_lev{L}=A13;
A14_lev{L}=A14;

A22_lev{L}=A22;
A23_lev{L}=A23;
A24_lev{L}=A24;

A33_lev{L}=A33;
A34_lev{L}=A34;
A44_lev{L}=A44;

C11_lev{L}=C11;
C12_lev{L}=C12;
C13_lev{L}=C13;
C14_lev{L}=C14;

C22_lev{L}=C22;
C23_lev{L}=C23;
C24_lev{L}=C24;

ND_ors=cell(L,1);
N_bc=cell(L,1);
E_bc=cell(L,1);
N_ord=cell(L,1);
E_ord=cell(L,1);

% for each level
for lev=L:-1:1
    
    % total number of edge
    NE=mesh{lev}.NE;
    % total number of nodes
    N=mesh{lev}.N;
    
    % number of edges that are on the boundary
    E_L_remove=length(mesh{lev}.E_remove);
    % number of nodes that are on the boundary    
    N_L_remove=length(mesh{lev}.N_remove); 

    % number of edges that are internal
    E_L_internal=NE-E_L_remove;
    % number of edges that are internal
    N_L_internal=N-N_L_remove;
    
    %  edges dofs that are on the boundary    
    E_beginning=mesh{lev}.E_remove;
    %  nodes dofs that are on the boundary    
    N_beginning=mesh{lev}.N_remove;
    
    % we are going to order the dofs in this way:
    % 1) edge bc, first component
    % 2) edge bc, second component
    % 3) node bc, first component
    % 4) node bc, second component
    % 5) edge internal, first component
    % 6) edge internal, second component
    % 7) node internal, first component
    % 8) node internal, second component    
    
    % save in L_lev the previous lengths
    L_lev{lev}(1)=E_L_remove;
    L_lev{lev}(2)=L_lev{lev}(1)+E_L_remove;
    L_lev{lev}(3)=L_lev{lev}(2)+N_L_remove;
    L_lev{lev}(4)=L_lev{lev}(3)+N_L_remove;
    L_lev{lev}(5)=L_lev{lev}(4)+E_L_internal;
    L_lev{lev}(6)=L_lev{lev}(5)+E_L_internal;
    L_lev{lev}(7)=L_lev{lev}(6)+N_L_internal;
    L_lev{lev}(8)=L_lev{lev}(7)+N_L_internal;
    
    % L_bc is the total length of the bc-vector
    L_bc=2*E_L_remove+2*N_L_remove;
    
    E1=1:E_L_internal; 
    N1=1:N_L_internal;
        
    E_order=1:mesh{lev}.NE;
    N_order=1:mesh{lev}.N;

    E1=E1+L_bc;
    E2=E1+E_L_internal;

    N1=N1+L_bc+2*E_L_internal;
    N2=N1+N_L_internal;
    
    ND_length=mesh{lev}.ND_lengths;
     
    E_order=cat(2,E_beginning,setdiff(E_order,E_beginning));
    N_order=cat(2,N_beginning,setdiff(N_order,N_beginning));
    ND_order=mesh{lev}.ND_nodes_reorder;
    
    E_ord{lev}=E_order;
    ND_ord{lev}=N_order;
    
    E_bc_loc_to_glob1=[];
    E_bc_loc_to_glob2=[];
    N_bc_loc_to_glob3=[];
    N_bc_loc_to_glob4=[];
    if(E_L_remove>0)
    E_bc_loc_to_glob1=1:2:(2*E_L_remove-1);
    E_bc_loc_to_glob2=2:2:(2*E_L_remove);
    end
    E_tmp=1:2:(2*E_L_internal-1);
    E_tmp=E_tmp+L_bc;
    E_bc_loc_to_glob1=cat(2,E_bc_loc_to_glob1,E_tmp);
    E_tmp=2:2:(2*E_L_internal);
    E_tmp=E_tmp+L_bc;
    E_bc_loc_to_glob2=cat(2,E_bc_loc_to_glob2,E_tmp);
    
    if(N_L_remove>0)
    N_bc_loc_to_glob3=1:2:(2*N_L_remove-1);
    N_bc_loc_to_glob4=2:2:(2*N_L_remove);    

    end
    N_bc_loc_to_glob3=N_bc_loc_to_glob3+2*E_L_remove;
    N_bc_loc_to_glob4=N_bc_loc_to_glob4+2*E_L_remove;
    N_tmp=1:2:(2*N_L_internal-1);
    N_tmp=N_tmp+L_bc+2*E_L_internal;
    N_bc_loc_to_glob3=cat(2,N_bc_loc_to_glob3,N_tmp);
    N_tmp=2:2:(2*N_L_internal);
    N_tmp=N_tmp+L_bc+2*E_L_internal;
    N_bc_loc_to_glob4=cat(2,N_bc_loc_to_glob4,N_tmp);    

    
    ND_bc_loc_to_glob1=1:2:(2*N-1);
    ND_bc_loc_to_glob2=2:2:(2*N);
    
    
    
    
    
    

    NDtoRT{lev}=NDtoRT{lev}(E_order,ND_order);
    
    if(lev==L)
    % change the order of the fine components
    if(L>1)
    RTCtoRTF{lev-1}=RTCtoRTF{lev-1}(E_order,:);
    P1CtoP1F{lev-1}=P1CtoP1F{lev-1}(N_order,:); 
    NDCtoNDF{lev-1}=NDCtoNDF{lev-1}(ND_order,:);
    end
    
    
 
      for ii=1:E_L_remove
      go_back(E_order(ii) )=2*ii-1;
      go_back(E_order(ii)+NE)  =2*ii;
      end
      
      for ii=1:N_L_remove
      go_back(N_order(ii)+2 * NE)  = 2*ii-1 + 2 * E_L_remove ;
      go_back(N_order(ii)+2 * NE + N)  =2*ii   + 2 * E_L_remove ;
      end
      
      L_bc=2*E_L_remove+2*N_L_remove;
      
      for ii=1:E_L_internal
      go_back(E_order(ii+E_L_remove))= 2*ii-1 + L_bc ;
      go_back(E_order(ii+E_L_remove)+NE)  =2*ii   + L_bc ;
      end
      L_bc_E=2*NE + 2 *N_L_remove;
      for ii=1:N_L_internal
      go_back(N_order(ii+N_L_remove)+2 * NE)=2*ii-1 + L_bc_E ;
      go_back(N_order(ii+N_L_remove)+2 * NE+N )  =2*ii   + L_bc_E;
      end




 
      

    b1=b1(E_order);
    b2=b2(E_order);
    b3=b3(N_order);
    b4=b4(N_order);
        
    %compute the rhs
    b_lev(E_bc_loc_to_glob1)=b1;
    b_lev(E_bc_loc_to_glob2)=b2;
    b_lev(N_bc_loc_to_glob3)=b3;
    b_lev(N_bc_loc_to_glob4)=b4;
    

    
    A11_lev{lev}=A11_lev{lev}(E_order,E_order);
    A12_lev{lev}=A12_lev{lev}(E_order,E_order);
    A13_lev{lev}=A13_lev{lev}(E_order,N_order);
    A14_lev{lev}=A14_lev{lev}(E_order,N_order);
    
    A22_lev{lev}=A22_lev{lev}(E_order,E_order);
    A23_lev{lev}=A23_lev{lev}(E_order,N_order);     
    A24_lev{lev}=A24_lev{lev}(E_order,N_order);    

    A33_lev{lev}=A33_lev{lev}(N_order,N_order);
    A34_lev{lev}=A34_lev{lev}(N_order,N_order);
    A44_lev{lev}=A44_lev{lev}(N_order,N_order);  
    
    C11_lev{lev}= C11_lev{lev}(ND_order,ND_order);
    C12_lev{lev}= C12_lev{lev}(ND_order,ND_order);
    C13_lev{lev}= C13_lev{lev}(ND_order,N_order);
    C14_lev{lev}= C14_lev{lev}(ND_order,N_order);
    
    C22_lev{lev}= C22_lev{lev}(ND_order,ND_order);
    C23_lev{lev}= C23_lev{lev}(ND_order,N_order);
    C24_lev{lev}= C24_lev{lev}(ND_order,N_order);
    else
    RTCtoRTF{lev}=RTCtoRTF{lev}(:,E_order);
    P1CtoP1F{lev}=P1CtoP1F{lev}(:,N_order); 
    NDCtoNDF{lev}=NDCtoNDF{lev}(:,ND_order);
    
    if(lev-1>0)
    RTCtoRTF{lev-1}=RTCtoRTF{lev-1}(E_order,:);
    P1CtoP1F{lev-1}=P1CtoP1F{lev-1}(N_order,:); 
    NDCtoNDF{lev-1}=NDCtoNDF{lev-1}(ND_order,:);
    end
    
    A11_lev{lev}=RTCtoRTF{lev}'*A11_lev{lev+1}*RTCtoRTF{lev};
    A12_lev{lev}=RTCtoRTF{lev}'*A12_lev{lev+1}*RTCtoRTF{lev};
    A13_lev{lev}=RTCtoRTF{lev}'*A13_lev{lev+1}*P1CtoP1F{lev};   
    A14_lev{lev}=RTCtoRTF{lev}'*A14_lev{lev+1}*P1CtoP1F{lev};
    
    A22_lev{lev}=RTCtoRTF{lev}'*A22_lev{lev+1}*RTCtoRTF{lev};
    A23_lev{lev}=RTCtoRTF{lev}'*A23_lev{lev+1}*P1CtoP1F{lev};   
    A24_lev{lev}=RTCtoRTF{lev}'*A24_lev{lev+1}*P1CtoP1F{lev};    

    A33_lev{lev}=P1CtoP1F{lev}'*A33_lev{lev+1}*P1CtoP1F{lev}; 
    A34_lev{lev}=P1CtoP1F{lev}'*A34_lev{lev+1}*P1CtoP1F{lev}; 
    A44_lev{lev}=P1CtoP1F{lev}'*A44_lev{lev+1}*P1CtoP1F{lev}; 
    
    C11_lev{lev}=NDCtoNDF{lev}'*C11_lev{lev+1}*NDCtoNDF{lev};
    C12_lev{lev}=NDCtoNDF{lev}'*C12_lev{lev+1}*NDCtoNDF{lev};
    C13_lev{lev}=NDCtoNDF{lev}'*C13_lev{lev+1}*NDCtoNDF{lev};
    C14_lev{lev}=NDCtoNDF{lev}'*C14_lev{lev+1}*NDCtoNDF{lev};
    
    C22_lev{lev}=NDCtoNDF{lev}'*C22_lev{lev+1}*NDCtoNDF{lev};
    C23_lev{lev}=NDCtoNDF{lev}'*C23_lev{lev+1}*NDCtoNDF{lev};
    C24_lev{lev}=NDCtoNDF{lev}'*C24_lev{lev+1}*NDCtoNDF{lev};
    
    
    end
          

    M_lev{lev}(E_bc_loc_to_glob1,E_bc_loc_to_glob1)=A11_lev{lev};
    M_lev{lev}(E_bc_loc_to_glob1,E_bc_loc_to_glob2)=A12_lev{lev};
    M_lev{lev}(E_bc_loc_to_glob1,N_bc_loc_to_glob3)=A13_lev{lev};
    M_lev{lev}(E_bc_loc_to_glob1,N_bc_loc_to_glob4)=A14_lev{lev};

    M_lev{lev}(E_bc_loc_to_glob2,E_bc_loc_to_glob1)= (A12_lev{lev})';
    M_lev{lev}(E_bc_loc_to_glob2,E_bc_loc_to_glob2)=A22_lev{lev};
    M_lev{lev}(E_bc_loc_to_glob2,N_bc_loc_to_glob3)=A23_lev{lev};
    M_lev{lev}(E_bc_loc_to_glob2,N_bc_loc_to_glob4)=A24_lev{lev};
    
    M_lev{lev}(N_bc_loc_to_glob3,E_bc_loc_to_glob1)=(A13_lev{lev})';
    M_lev{lev}(N_bc_loc_to_glob3,E_bc_loc_to_glob2)=(A23_lev{lev})';    
    M_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob3)=A33_lev{lev};
    M_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob4)=A34_lev{lev};

    M_lev{lev}(N_bc_loc_to_glob4,E_bc_loc_to_glob1)=(A14_lev{lev})'; 
    M_lev{lev}(N_bc_loc_to_glob4,E_bc_loc_to_glob2)=(A24_lev{lev})';
    M_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob3)=(A34_lev{lev})';
    M_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob4)=A44_lev{lev};   
    
    ND_bc_loc_to_glob3=2*N+1:3*N;
    ND_bc_loc_to_glob4=3*N+1:4*N;
        
    C_lev{lev}(ND_bc_loc_to_glob1,ND_bc_loc_to_glob1)=C11_lev{lev};
    C_lev{lev}(ND_bc_loc_to_glob1,ND_bc_loc_to_glob2)=C12_lev{lev};
    C_lev{lev}(ND_bc_loc_to_glob1,ND_bc_loc_to_glob3)=C13_lev{lev};
    C_lev{lev}(ND_bc_loc_to_glob1,ND_bc_loc_to_glob4)=C14_lev{lev};

    C_lev{lev}(ND_bc_loc_to_glob2,ND_bc_loc_to_glob1)= (C12_lev{lev})';
    C_lev{lev}(ND_bc_loc_to_glob2,ND_bc_loc_to_glob2)=C22_lev{lev};
    C_lev{lev}(ND_bc_loc_to_glob2,ND_bc_loc_to_glob3)=C23_lev{lev};
    C_lev{lev}(ND_bc_loc_to_glob2,ND_bc_loc_to_glob4)=C24_lev{lev};
     
%     C_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob1)=(C13_lev{lev})';
%     C_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob2)=(C23_lev{lev})';    
%     C_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob3)=A33_lev{lev};
%     C_lev{lev}(N_bc_loc_to_glob3,N_bc_loc_to_glob4)=A34_lev{lev};
% 
%     C_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob1)=(C14_lev{lev})'; 
%     C_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob2)=(C24_lev{lev})';
%     C_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob3)=(A34_lev{lev})';
%     C_lev{lev}(N_bc_loc_to_glob4,N_bc_loc_to_glob4)=A44_lev{lev};
    
    
        
    
    
end

end