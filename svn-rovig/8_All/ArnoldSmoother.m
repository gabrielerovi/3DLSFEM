


function x=ArnoldSmoother(top2bottom,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,Boundary_Node, Boundary_Edge,...
    x,b,mesh,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44, ...
    A11_lev,A12_lev,A13_lev,A14_lev,A22_lev,A23_lev,A24_lev,A33_lev,A34_lev,A44_lev,smoothing_steps,eta)
N=mesh.N;
NE=mesh.NE;
N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;

LGlob=[NE;2*NE;2*NE+N;2*NE+2*N ];


b_sigma1=b(1:LGlob(1));
b_sigma2=b(LGlob(1)+1:LGlob(2));
b_disp1 =b(LGlob(2)+1:LGlob(3));
b_disp2 =b(LGlob(3)+1:LGlob(4));

if(top2bottom)
    vertices=1:N;
else
    vertices=N:-1:1;
end


%  SMOOTHING-STEPS
for jj=1:smoothing_steps
 
x_old=x;    
sigma1=x(1:LGlob(1));
sigma2=x(LGlob(1)+1:LGlob(2));
disp1 =x(LGlob(2)+1:LGlob(3));
disp2 =x(LGlob(3)+1:LGlob(4));



EcontBC=0;   
type_of_dof=2;
% we have no contact right now
contact=0;
for ii=E_remove
    
    A2x2=[A11_lev(ii,ii)   A12_lev(ii,ii);
          A12_lev(ii,ii) A22_lev(ii,ii)];
   
    ww=[1:ii-1,ii+1:NE];
    b_tmp= A11_lev(ii,ww) * sigma1(ww) + A12_lev(ii,ww) * sigma2(ww) + A13_lev(ii,:) *disp1 + A14_lev(ii,:) * disp2;
    b2x2(1,1)=b_sigma1(ii)   - b_tmp;
    b_tmp= A12_lev(ww,ii)' * sigma1(ww) + A22_lev(ii,ww) * sigma2(ww) + A23_lev(ii,:) *disp1 + A24_lev(ii,:) * disp2;
    b2x2(2,1)=b_sigma2(ii)   - b_tmp;

    U_xy=solve_2x2system(A2x2,b2x2);
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), contact);  
    % rescale the boundary condition
    % ATTENZIONE.   NON DEVI FARLO ALTROVE, PRIMA DI ADD BOUNDARY? 
    % DEVI CONFRONTARE GLI SFORZI, QUINDI NELLA FUNZIONE DEVI VALUTARE
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh);
    sigma1(ii)=tmp(1)/coeff;
    sigma2(ii)=tmp(2)/coeff; 
       
end

NcontBC=0;   
type_of_dof=1;
for ii=N_remove
    
    A2x2=[A33_lev(ii,ii)   A34_lev(ii,ii);
       A44_lev(ii,ii)   A44_lev(ii,ii)];
   
    ww=[1:ii-1,ii+1:N];
    b_tmp= A33_lev(ii,ww) * disp1(ww) + A34_lev(ii,ww) * disp2(ww) + A13_lev(:,ii)' * sigma1 + A23_lev(:,ii)' * sigma2;
    b2x2(1,1)=b_disp1(ii)   - b_tmp;
    b_tmp= A34_lev(ww,ii)' * disp1(ww) + A44_lev(ii,ww) * disp2(ww) + A14_lev(:,ii)' * sigma1 + A24_lev(:,ii)' * sigma2;
    b2x2(2,1)=b_disp2(ii)   - b_tmp;

    U_xy=solve_2x2system(A2x2,b2x2);
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), contact);   
    disp1(ii)=tmp(1);
    disp2(ii)=tmp(2);
       
end

for nn=vertices

    NELoc=length(A11{nn}(1,:));
    NLoc=length(A33{nn}(1,:));
    
    % ##############################################################################
    % ATTENZIONE LI SETTI A ZERO, POI METTI LE BC, MA NON METTI I TERMINI
    % VOLUMETRICII!!! CRISTOOOOOOOOOOO00000OOOOIIIIAHAHAAHAH
    

    % ##############################################################################
    
    
    % border vertex - edge, Local - Global
    N_borderGlob=Boundary_Node{nn};
    E_borderGlob=Boundary_Edge{nn};
    N_borderLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_borderGlob,1)));
    E_borderLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_borderGlob,1)));
    
    % all vertex -edge, Local - Global
    N_Loc=1:NLoc;
    E_Loc=1:NELoc;   
    N_Glob=cell2mat(values(NmapLoc2Glob{nn},num2cell(N_Loc,1)));
    E_Glob=cell2mat(values(EmapLoc2Glob{nn},num2cell(E_Loc,1)));
    
    % internal vertex- edge, Local - Global
    N_Loc_internal=setdiff(N_Loc,N_borderLoc);
    N_Glob_internal=setdiff(N_Glob,N_borderGlob);
    E_Loc_internal=setdiff(E_Loc,E_borderLoc);
    E_Glob_internal=setdiff(E_Glob,E_borderGlob);    
    

    
    
    LLoc=[NELoc; 2*NELoc; 2*NELoc+NLoc;2*NELoc+2*NLoc;];
    
    a11=A11{nn}; a12=A12{nn}; a13=A13{nn}; a14=A14{nn};
    a21=a12';    a22=A22{nn}; a23=A23{nn}; a24=A24{nn}; 
    a31=a13';    a32=a23';    a33=A33{nn}; a34=A34{nn}; 
    a41=a14';    a42=a24';    a43=a34';    a44=A44{nn};

    % assign the global external forces to the rhs
    b_sigma1Loc=b_sigma1(E_Glob);
    b_sigma2Loc=b_sigma2(E_Glob);
    b_disp1Loc=b_disp1(N_Glob);
    b_disp2Loc=b_disp2(N_Glob);
    
    % then subtract the system computed wrt x_old
    b_sigma1Loc=b_sigma1Loc - a11 * sigma1(E_Glob) - a12 * sigma2(E_Glob)- a13 * disp1(N_Glob) - a14 * disp2(N_Glob);
    b_sigma2Loc=b_sigma2Loc - a21 * sigma1(E_Glob) - a22 * sigma2(E_Glob)- a23 * disp1(N_Glob) - a24 * disp2(N_Glob);
    b_disp1Loc=b_disp1Loc   - a31 * sigma1(E_Glob) - a32 * sigma2(E_Glob)- a33 * disp1(N_Glob) - a34 * disp2(N_Glob);
    b_disp2Loc=b_disp2Loc   - a41 * sigma1(E_Glob) - a42 * sigma2(E_Glob)- a43 * disp1(N_Glob) - a44 * disp2(N_Glob);
    

    
    cont1=0;
    for kk=E_borderLoc
        cont1=cont1+1;
        a11(kk,:)=0; a11(kk,kk)=1; 
        b_sigma1Loc(kk,1)= 0 ;    %%%%%%%%%%%%%%%%b_sigma1(E_borderGlob(cont1)); STO STUDIANDO IL PROBLEMA DIFFERENZA
        a22(kk,:)=0; a22(kk,kk)=1; 
        b_sigma2Loc(kk,1)=0; %%%%%%%%%%%%%%%5b_sigma2(E_borderGlob(cont1));
        a12(kk,:)=0; a13(kk,:)=0; a14(kk,:)=0; 
        a21(kk,:)=0; a23(kk,:)=0; a24(kk,:)=0;
        
    end
    cont1=0;
    for kk=N_borderLoc
        cont1=cont1+1;
        a33(kk,:)=0; a33(kk,kk)=1; 
        b_disp1Loc(kk,1)=   0;    %%%%%%%%%%%%b_disp1(N_borderGlob(cont1));
        a44(kk,:)=0; a44(kk,kk)=1; 
        b_disp2Loc(kk,1)=   0; %%%%%%%%%%%%%b_disp2(N_borderGlob(cont1));
        a31(kk,:)=0; a32(kk,:)=0; a34(kk,:)=0; 
        a41(kk,:)=0; a42(kk,:)=0; a43(kk,:)=0;   
    end
    
    A_Loc=[a11 a12 a13 a14;
           a21 a22 a23 a24;
           a31 a32 a33 a34;
           a41 a42 a43 a44;];
       
       %  A PARTE LE BC COME CALCOLI BLOC?
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    x_Loc=A_Loc\b_Loc;
    sigma1Loc=x_Loc(1:LLoc(1));
    sigma2Loc=x_Loc(LLoc(1)+1:LLoc(2));
    disp1Loc =x_Loc(LLoc(2)+1:LLoc(3));
    disp2Loc =x_Loc(LLoc(3)+1:LLoc(4));
    
%     sigma1(E_Glob_internal)=  eta* sigma1Loc(E_Loc_internal);
%     sigma2(E_Glob_internal)=  eta* sigma2Loc(E_Loc_internal);
%     disp1(N_Glob_internal) =  eta* disp1Loc(N_Loc_internal);
%     disp2(N_Glob_internal) =  eta* disp2Loc(N_Loc_internal);
coeff=1;
    sigma1(E_Glob_internal)=  coeff *(sigma1(E_Glob_internal)+ eta * sigma1Loc(E_Loc_internal) );
    sigma2(E_Glob_internal)=  coeff *(sigma2(E_Glob_internal)+ eta * sigma2Loc(E_Loc_internal) );
    disp1(N_Glob_internal) =  coeff *(disp1(N_Glob_internal)+ eta *disp1Loc(N_Loc_internal) );
    disp2(N_Glob_internal) =  coeff *(disp2(N_Glob_internal)+ eta * disp2Loc(N_Loc_internal) );
    
end
x=[sigma1;sigma2;disp1;disp2];
%norm(x_old-x);

end


end

