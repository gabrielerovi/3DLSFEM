function x = ArnoldVcycle(lev,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,Boundary_Node, Boundary_Edge,...
    x,b,mesh,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44, ...
    A11_lev,A12_lev,A13_lev,A14_lev,A22_lev,A23_lev,A24_lev,A33_lev,A34_lev,A44_lev,smoothing_steps,eta, ...
    RTCtoRTF,P1CtoP1F)

 % fine and coarse
 F=lev;
 C=lev-1;
 
 % N and NE on the fine mesh
 NE=mesh{F}.NE;
 N=mesh{F}.N;
 E_label=mesh{F}.E_label;
 N_label=mesh{F}.N_label;
 % fine level lenghts
 LFine=[mesh{F}.NE, 2*mesh{F}.NE, 2*mesh{F}.NE+mesh{F}.N, 2*mesh{F}.NE+2*mesh{F}.N ];
 
 
  E_remove=mesh{lev}.E_remove;
  N_remove=mesh{lev}.N_remove;    
  E_dirichlet=mesh{lev}.E_dirichlet;
  N_dirichlet=mesh{lev}.N_dirichlet;

  % create matrix at the rough level
  A=[A11_lev{lev}  A12_lev{lev}  A13_lev{lev}  A14_lev{lev};
     A12_lev{lev}' A22_lev{lev}  A23_lev{lev}  A24_lev{lev}; 
     A13_lev{lev}' A23_lev{lev}' A33_lev{lev}  A34_lev{lev};
     A14_lev{lev}' A24_lev{lev}' A34_lev{lev}' A44_lev{lev};];
 
  
  sigma1=x(1:LFine(1));
  sigma2=x(1+LFine(1):LFine(2));
  disp1= x(1+LFine(2):LFine(3));  
  disp2= x(1+LFine(3):LFine(4));  
  
contact=0;
for ii=E_remove
    
    A2x2=[A11_lev{lev}(ii,ii)   A12_lev{lev}(ii,ii);
          A12_lev{lev}(ii,ii)   A22_lev{lev}(ii,ii)];
   
    ww=[1:ii-1,ii+1:NE];
    b_tmp= A11_lev{lev}(ii,ww) * sigma1(ww) + A12_lev{lev}(ii,ww) * sigma2(ww) + A13_lev{lev}(ii,:) *disp1 + A14_lev{lev}(ii,:) * disp2;
    b2x2(1,1)=b(ii)   - b_tmp;
    b_tmp= A12_lev{lev}(ww,ii)' * sigma1(ww) + A22_lev{lev}(ii,ww) * sigma2(ww) + A23_lev{lev}(ii,:) *disp1 + A24_lev{lev}(ii,:) * disp2;
    b2x2(2,1)=b(ii + NE)   - b_tmp;

    U_xy=solve_2x2system(A2x2,b2x2);
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), contact);  
    % rescale the boundary condition
    % ATTENZIONE.   NON DEVI FARLO ALTROVE, PRIMA DI ADD BOUNDARY? 
    % DEVI CONFRONTARE GLI SFORZI, QUINDI NELLA FUNZIONE DEVI VALUTARE
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh);
    sigma1(ii)=tmp(1)/coeff;
    sigma2(ii)=tmp(2)/coeff; 
    A(ii,:)= 0;     A(ii,ii)= 1;
    jj=ii+NE;
    A(jj,:)= 0;     A(jj,jj)= 1;
end

NcontBC=0;   
type_of_dof=1;
for ii=N_remove
    
    A2x2=[A33_lev{lev}(ii,ii)   A34_lev{lev}(ii,ii);
          A44_lev{lev}(ii,ii)   A44_lev{lev}(ii,ii)];
   
    ww=[1:ii-1,ii+1:N];
    b_tmp= A33_lev{lev}(ii,ww) * disp1(ww) + A34_lev{lev}(ii,ww) * disp2(ww) + A13_lev{lev}(:,ii)' * sigma1 + A23_lev{lev}(:,ii)' * sigma2;
    b2x2(1,1)=b(ii+2*NE)   - b_tmp;
    b_tmp= A34_lev{lev}(ww,ii)' * disp1(ww) + A44_lev{lev}(ii,ww) * disp2(ww) + A14_lev{lev}(:,ii)' * sigma1 + A24_lev{lev}(:,ii)' * sigma2;
    b2x2(2,1)=b(ii+2*NE+N)   - b_tmp;

    U_xy=solve_2x2system(A2x2,b2x2);
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), contact);   
    disp1(ii)=tmp(1);
    disp2(ii)=tmp(2);
    jj=ii+2*NE;
    A(jj,:)= 0;     A(jj,jj)= 1;
    jj=jj+N;
    A(jj,:)= 0;     A(jj,jj)= 1;       
end














if(lev==1)
  % we impose homogeneous dirichlet for the correction where we impose (and
  % know) dirichlet bc for the original proble
% b(E_remove)=0;   b(E_remove2)=0;  b(N_remove1)=0;   b(N_remove2)=0;
 x = A\ b;
 
else
 % coarse level lengths
 LCoarse=[mesh{C}.NE, 2*mesh{C}.NE, 2*mesh{C}.NE+mesh{C}.N, 2*mesh{C}.NE+2*mesh{C}.N ];  
% presmoothing
top2bottom=1;
 x=ArnoldSmoother(top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},Boundary_Node{lev}, Boundary_Edge{lev},...
    x,b,mesh{lev},A11{lev},A12{lev},A13{lev},A14{lev},A22{lev},A23{lev},A24{lev},A33{lev},A34{lev},A44{lev},...
    A11_lev{lev},A12_lev{lev},A13_lev{lev},A14_lev{lev},A22_lev{lev},A23_lev{lev},A24_lev{lev},A33_lev{lev},A34_lev{lev},A44_lev{lev},...
smoothing_steps,eta);


resF= b - A* x;
resC(1:LCoarse(1), 1)=            RTCtoRTF{C}' * resF(1:LFine(1));
resC(1+LCoarse(1):LCoarse(2), 1)= RTCtoRTF{C}' * resF(1+LFine(1):LFine(2));
resC(1+LCoarse(2):LCoarse(3), 1)= P1CtoP1F{C}' * resF(1+LFine(2):LFine(3));
resC(1+LCoarse(3):LCoarse(4), 1)= P1CtoP1F{C}' * resF(1+LFine(3):LFine(4));


coorection=zeros(LCoarse(4),1);
% V-Cycle 
coorection = ArnoldVcycle(lev-1,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,Boundary_Node, Boundary_Edge,...
    coorection,resC,mesh,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44, ...
    A11_lev,A12_lev,A13_lev,A14_lev,A22_lev,A23_lev,A24_lev,A33_lev,A34_lev,A44_lev,smoothing_steps,eta, ...
    RTCtoRTF,P1CtoP1F);

x(1:LFine(1))          = x(1:LFine(1))          + RTCtoRTF{C} * coorection(1:LCoarse(1));
x(1+LFine(1):LFine(2)) = x(1+LFine(1):LFine(2)) + RTCtoRTF{C} * coorection(1+LCoarse(1):LCoarse(2));
x(1+LFine(2):LFine(3)) = x(1+LFine(2):LFine(3)) + P1CtoP1F{C} * coorection(1+LCoarse(2):LCoarse(3));
x(1+LFine(3):LFine(4)) = x(1+LFine(3):LFine(4)) + P1CtoP1F{C} * coorection(1+LCoarse(3):LCoarse(4));


% postsmoothing. Let us do a symmetric v-cycle!! yeah!
top2bottom=0;
 x=ArnoldSmoother(top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},Boundary_Node{lev}, Boundary_Edge{lev},...
    x,b,mesh{lev},A11{lev},A12{lev},A13{lev},A14{lev},A22{lev},A23{lev},A24{lev},A33{lev},A34{lev},A44{lev},...
    A11_lev{lev},A12_lev{lev},A13_lev{lev},A14_lev{lev},A22_lev{lev},A23_lev{lev},A24_lev{lev},A33_lev{lev},A34_lev{lev},A44_lev{lev},...
smoothing_steps,eta);

end


end


