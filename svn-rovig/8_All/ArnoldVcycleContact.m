function x = ArnoldVcycleContact(graph,maxlev,L,maps,x,bnt,mesh,Ant,P, parameters,Constraint)

EmapGlob2Loc=maps.EmapGlob2Loc;
EmapLoc2Glob=maps.EmapLoc2Glob;
NmapGlob2Loc=maps.NmapGlob2Loc;
NmapLoc2Glob=maps.NmapLoc2Glob;
Patch_Boundary_Edge=maps.Patch_Boundary_Edge;
Patch_Boundary_Node=maps.Patch_Boundary_Node;
Patch_Edge=maps.Patch_Edge;
Patch_Node=maps.Patch_Node;

 % fine and coarse
 F=L;
 C=L-1;
 
 % N and NE on the fine mesh
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 
 % N and NE on the coarse mesh
 NEC=mesh{C}.NE;
 NC=mesh{C}.N;
 
 %  E_label=mesh{F}.E_label;
%  N_label=mesh{F}.N_label;
 % fine level lenghts
 LFine=[ NEF ; 2 * NEF; 2 * NEF + NF ; 2 * NEF + 2 * NF ];
 LCoarse=[NEC; 2 * NEC;2 * NEC + NC ; 2 * NEC + 2 * NC ];  

 
  E_remove=mesh{F}.E_remove;
  E_dirichlet=mesh{F}.E_dirichlet;

  N_remove=mesh{F}.N_remove;
  N_dirichlet=mesh{F}.N_dirichlet;
  
smoothing_steps=parameters.smoothing_steps;

    
E_remove=mesh{F}.E_remove;
N_remove=mesh{F}.N_remove;    
E_contact_tangent=mesh{F}.E_contact;  
% L=lev;
% A11_lev{L}=Ant{1,1}; 
% A12_lev{L}=Ant{1,2}; 
% A13_lev{L}=Ant{1,3}; 
% A14_lev{L}=Ant{1,4}; 
% 
% A21_lev{L}=Ant21; 
% A22_lev{L}=Ant22; 
% A23_lev{L}=Ant23; 
% A24_lev{L}=Ant24; 
% 
% A31_lev{L}=Ant31; 
% A32_lev{L}=Ant32;  
% A33_lev{L}=Ant33;  
% A34_lev{L}=Ant34;  
% 
% A41_lev{L}=Ant41;
% A42_lev{L}=Ant42; 
% A43_lev{L}=Ant43;  
% A44_lev{L}=Ant44;


% Ant has no bc for the moment. therefore we copy them in Anobc that will
% be used to make the assembly on the coarser levels for solving the linear
% problem.
% Ant will be instead used just for the finer level
for kk=1:4
    for mm=1:4
        Anobc{F,kk,mm}=Ant{kk,mm};
    end
end





for nn=1:2
    for mm=1:4
       Ant{nn,mm} (E_remove,:)=0; 
       Ant{nn+2,mm} (N_remove,:)=0;
    end
    
    for eee1=E_remove
    Ant{nn,nn}(eee1,eee1)=1;
    end
    
    for nnn1=N_remove
    Ant{nn+2,nn+2}(nnn1,nnn1)=1;
    end
    
end

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % here we fix everything also on GammaC
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for mm=1:4        
       Ant{1,mm} (Constraint.WorkingSetE,:)=0; 
    end
    for eee=Constraint.WorkingSetE
        Ant{1,1}(eee,eee)=1;
    end
    for mm=1:4        
       Ant{3,mm} (Constraint.WorkingSetN,:)=0; 
    end
    for nnn=Constraint.WorkingSetN
        Ant{3,3}(nnn,nnn)=1;
    end
    
    
A=[ Ant{1,1}, Ant{1,2}, Ant{1,3}, Ant{1,4};
    Ant{2,1}, Ant{2,2}, Ant{2,3}, Ant{2,4};
    Ant{3,1}, Ant{3,2}, Ant{3,3}, Ant{3,4};
    Ant{4,1}, Ant{4,2}, Ant{4,3}, Ant{4,4}];


% we assume to have at least two levels
% so at the moment we are not on coarser grid
is_on_coarser_grid=0;
top2bottom=1;    
%[x,Constraint]=ArnoldSmootherContact3(x,bnt,mesh,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)

[x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},Ant,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);

 
%  
[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L},1);
sol=M_Normal_Tangent*x;
print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');

% [x,WorkingSetNormal_N,WorkingSetNormal_E]=ArnoldSmootherContactNormalTangent...
%     (graph{F},top2bottom,EmapGlob2Loc{F},EmapLoc2Glob{F},NmapGlob2Loc{F},NmapLoc2Glob{F},...
%     Patch_Boundary_Edge{F},Patch_Boundary_Node{F}, Patch_Edge{F},Patch_Node{F},x,b,mesh{F},...
%     Ant{1,1}, Ant{1,2}, Ant{1,3}, Ant{1,4},...
%     Ant{2,1}, Ant{2,2}, Ant{2,3}, Ant{2,4},...
%     Ant{3,1}, Ant{3,2}, Ant{3,3}, Ant{3,4},...
%     Ant{4,1}, Ant{4,2}, Ant{4,3}, Ant{4,4},...
%     smoothing_steps,0);


[Acoarse,Pcoarse]=ArnoldAssembling(Anobc,P,[],[],mesh);

%[Acoarse,Pcoarse]=ArnoldAssembling(Anobc,P,Constraint.WorkingSetE,Constraint.WorkingSetN,mesh);

% for lev=1:C
%     Acoarse11{lev}=Anobc{lev,1,1};
%     Acoarse12{lev}=Anobc{lev,1,2};
%     Acoarse13{lev}=Anobc{lev,1,3};
%     Acoarse14{lev}=Anobc{lev,1,4};
% 
% 
%     Acoarse21{lev}=Anobc{lev,2,1};
%     Acoarse22{lev}=Anobc{lev,2,2};
%     Acoarse23{lev}=Anobc{lev,2,3};
%     Acoarse24{lev}=Anobc{lev,2,4};
% 
%     Acoarse31{lev}=Anobc{lev,3,1};
%     Acoarse32{lev}=Anobc{lev,3,2};
%     Acoarse33{lev}=Anobc{lev,3,3};
%     Acoarse34{lev}=Anobc{lev,3,4};
%     
%     Acoarse41{lev}=Anobc{lev,4,1};
%     Acoarse42{lev}=Anobc{lev,4,2};
%     Acoarse43{lev}=Anobc{lev,4,3};
%     Acoarse44{lev}=Anobc{lev,4,4};
% 
% end

% Acoarse=[Acoarse11{C},Acoarse12{C},Acoarse13{C},Acoarse14{C};
%          Acoarse21{C},Acoarse22{C},Acoarse23{C},Acoarse24{C};
%          Acoarse31{C},Acoarse32{C},Acoarse33{C},Acoarse34{C};
%          Acoarse41{C},Acoarse42{C},Acoarse43{C},Acoarse44{C}];

% compute fine residual on the fine level
resF= bnt - A * x;  
resF(Constraint.WorkingSetE)=0;
resF(Constraint.WorkingSetN+2*NEF)=0;

% then project it to the lower level
% now we want to solve a linear problem. all the degrees of freedom that
% are constrained on the fine levelwill not influence coarse levels
% for this reason, in normal direction we use RTCtoRTF_normal and P1CtoP1F_normal

resC=Pcoarse{C}'*resF;

E_removeC=mesh{C}.E_remove;
E_removeC1=E_removeC;
E_removeC2=E_removeC1+NEC;

N_removeC=mesh{C}.N_remove;
N_removeC1 = N_removeC + 2 * NEC;
N_removeC2 = N_removeC1+ NC;


removeCTOT=[E_removeC1,E_removeC2,N_removeC1,N_removeC2];

resC(removeCTOT,1)=0;

%%%%% qui fai Vcycle
corrC=sparse(2*NEC+2*NC,1);
LineSearch=0;

 corrC = ArnoldVcycleLinearElasticNormalTangent(graph,maxlev,C,maps,corrC,resC,mesh,...
    Acoarse,smoothing_steps,Pcoarse,LineSearch);


corrF=Pcoarse{C}*corrC;

E_removeF=mesh{F}.E_remove;
E_removeF1=E_removeF;
E_removeF2=E_removeF1+NEF;

N_removeF=mesh{F}.N_remove;
N_removeF1 = N_removeF + 2 * NEF;
N_removeF2 = N_removeF1+ NF;


removeCTOT=[E_removeF1,E_removeF2,N_removeF1,N_removeF2];

corrF(removeCTOT)=0;


% qui devi correggere corrF

% Check if correction satisfy constraint. If not, project it.
%corrF=ContactProjection(x,corrF,Constraint,mesh{F});

x=x+corrF;

sol=M_Normal_Tangent*x;
print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');





% resC(1:LCoarse(1),1)         =       RTCtoRTF_normal' * resF(1:LFine(1));
% resC(1+LCoarse(1): LCoarse(2),1)=    P.RTCtoRTF{C}' * resF(1+LFine(1):LFine(2));
% resC(1+LCoarse(2): LCoarse(3),1)=    P1CtoP1F_normal' * resF(1+LFine(2):LFine(3));
% resC(1+LCoarse(3): LCoarse(4),1)=    P.P1CtoP1F{C}' * resF(1+LFine(3):LFine(4));

% now we want to be sure that the bc on the coarser level are zero
% E_removeC=mesh{C}.E_remove;
% E_removeC1=E_removeC;
% E_removeC2=E_removeC1+NEC;
% E_removeCTOT=[E_removeC1,E_removeC2];
% 
% N_removeC=mesh{C}.N_remove;
% N_removeC1 = N_removeC + 2 * NEC;
% N_removeC2 = N_removeC1+ NC;
% 
% N_removeCTOT=[N_removeC1,N_removeC2];
% 
% removeCTOT=[E_removeCTOT,N_removeCTOT];
% 
% resC(removeCTOT,1)=0;
% 
% correction=zeros(LCoarse(end),1);


% Linear  V-Cycle 
% LineSearch=0;
% for mm=1:parameters.MGiter
%  correction = ArnoldVcycle3(graph,maxlev,C,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
%     Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,correction,resC,mesh,...
%     Acoarse11,Acoarse12,Acoarse13,Acoarse14,...
%     Acoarse21,Acoarse22,Acoarse23,Acoarse24,...
%     Acoarse31,Acoarse32,Acoarse33,Acoarse34,...
%     Acoarse41,Acoarse42,Acoarse43,Acoarse44,...
%     smoothing_steps,P.RTCtoRTF,P.P1CtoP1F,LineSearch);
% 
% norm(resC-Acoarse*correction)
% if(norm(resC-Acoarse*correction)<parameters.toll)
%     break
% end
% end
% 



% corrF1=RTCtoRTF_normal * correction(1:LCoarse(1));
% corrF2=P.RTCtoRTF{C} * correction(1+LCoarse(1):LCoarse(2));
% corrF3=P1CtoP1F_normal * correction(1+LCoarse(2):LCoarse(3));
% corrF4=P.P1CtoP1F{C} * correction(1+LCoarse(3):LCoarse(4));
% 
% E_removeF=mesh{F}.E_remove;
% N_removeF=mesh{F}.N_remove;
% 
% corrF1(E_removeF)=0;
% corrF2(E_removeF)=0;
% corrF3(N_removeF)=0;
% corrF4(N_removeF)=0;
% 
% if(LineSearch==true)
% corrF=[corrF1;corrF2;corrF3;corrF4];
% Ac=A*corrF;
% alpha = (resF'*corrF) /(Ac'*corrF);
% else
% alpha=1;
% end
% 
% 
% x(1:LFine(1))          = x(1:LFine(1))          + corrF1;
% x(1+LFine(1):LFine(2)) = x(1+LFine(1):LFine(2)) + corrF2;
% x(1+LFine(2):LFine(3)) = x(1+LFine(2):LFine(3)) + corrF3;
% x(1+LFine(3):LFine(4)) = x(1+LFine(3):LFine(4)) + corrF4;
% % postsmoothing. Let us do a symmetric v-cycle!! yeah!
% top2bottom=0;

is_on_coarser_grid=0;
%[x,Constraint]=ArnoldSmootherContact3(x,bnt,mesh,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)
 [x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},Ant,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)

sol=M_Normal_Tangent*x;
print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');

% end


end


