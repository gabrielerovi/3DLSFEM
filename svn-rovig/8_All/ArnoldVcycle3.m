function x = ArnoldVcycle3(graph,maxlev,lev,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,A12_lev,A13_lev,A14_lev,...
    A21_lev,A22_lev,A23_lev,A24_lev,...
    A31_lev,A32_lev,A33_lev,A34_lev,...
    A41_lev,A42_lev,A43_lev,A44_lev,...
    smoothing_steps,RTCtoRTF,P1CtoP1F,LineSearch)

 % fine and coarse
 F=lev;
 C=lev-1;
 
 % N and NE on the fine mesh
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 
 E_label=mesh{F}.E_label;
 N_label=mesh{F}.N_label;
 % fine level lenghts
 LFine=[ NEF ; 2 * NEF; 2 * NEF + NF ; 2 * NEF + 2 * NF ];
 
 
  E_remove=mesh{lev}.E_remove;
  E_dirichlet=mesh{lev}.E_dirichlet;

  N_remove=mesh{lev}.N_remove;
  N_dirichlet=mesh{lev}.N_dirichlet;
  
  % create matrix at the rough level
 A=[A11_lev{lev} A12_lev{lev} A13_lev{lev} A14_lev{lev};
    A21_lev{lev} A22_lev{lev} A23_lev{lev} A24_lev{lev};
    A31_lev{lev} A32_lev{lev} A33_lev{lev} A34_lev{lev};
    A41_lev{lev} A42_lev{lev} A43_lev{lev} A44_lev{lev};];
 

  
  if(maxlev==lev)
    is_on_coarser_grid=false;
else
    is_on_coarser_grid=true;
end



if(lev==1)

    x = A\ b;

else
 % coarse level lengths
 NEC=mesh{C}.NE;
 NC=mesh{C}.N;
 
 LCoarse=[NEC; 2 * NEC;2 * NEC + NC ; 2 * NEC + 2 * NC ];  
% presmoothing
top2bottom=1;    
x=ArnoldSmoother3(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
    A11_lev{lev},A12_lev{lev},A13_lev{lev},A14_lev{lev},...
    A21_lev{lev},A22_lev{lev},A23_lev{lev},A24_lev{lev},...
    A31_lev{lev},A32_lev{lev},A33_lev{lev},A34_lev{lev},...
    A41_lev{lev},A42_lev{lev},A43_lev{lev},A44_lev{lev},...
   smoothing_steps,is_on_coarser_grid);



resF= b - A* x;
resC(1:LCoarse(1),1)         =       RTCtoRTF{C}' * resF(1:LFine(1));
resC(1+LCoarse(1): LCoarse(2),1)=    RTCtoRTF{C}' * resF(1+LFine(1):LFine(2));
resC(1+LCoarse(2): LCoarse(3),1)=    P1CtoP1F{C}' * resF(1+LFine(2):LFine(3));
resC(1+LCoarse(3): LCoarse(4),1)=    P1CtoP1F{C}' * resF(1+LFine(3):LFine(4));


E_removeC=mesh{C}.E_remove;
E_removeC1=E_removeC;
E_removeC2=E_removeC1+NEC;
E_removeCTOT=[E_removeC1,E_removeC2];

N_removeC=mesh{C}.N_remove;
N_removeC1 = N_removeC + 2 * NEC;
N_removeC2 = N_removeC1+ NC;

N_removeCTOT=[N_removeC1,N_removeC2];

removeCTOT=[E_removeCTOT,N_removeCTOT];

resC(removeCTOT,1)=0;

correction=zeros(LCoarse(end),1);
% V-Cycle 
correction = ArnoldVcycle3(graph,maxlev,lev-1,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,correction,resC,mesh,...
    A11_lev,A12_lev,A13_lev,A14_lev,...
    A21_lev,A22_lev,A23_lev,A24_lev,...
    A31_lev,A32_lev,A33_lev,A34_lev,...
    A41_lev,A42_lev,A43_lev,A44_lev,...
    smoothing_steps,RTCtoRTF,P1CtoP1F,LineSearch);


corrF1=RTCtoRTF{C} * correction(1:LCoarse(1));
corrF2=RTCtoRTF{C} * correction(1+LCoarse(1):LCoarse(2));
corrF3=P1CtoP1F{C} * correction(1+LCoarse(2):LCoarse(3));
corrF4=P1CtoP1F{C} * correction(1+LCoarse(3):LCoarse(4));

E_removeF=mesh{F}.E_remove;
N_removeF=mesh{F}.N_remove;

corrF1(E_removeF)=0;
corrF2(E_removeF)=0;
corrF3(N_removeF)=0;
corrF4(N_removeF)=0;

if(LineSearch==true)
corrF=[corrF1;corrF2;corrF3;corrF4];
Ac=A*corrF;
alpha = (resF'*corrF) /(Ac'*corrF);
else
alpha=1;
end


x(1:LFine(1))          = x(1:LFine(1))          + corrF1;
x(1+LFine(1):LFine(2)) = x(1+LFine(1):LFine(2)) + corrF2;
x(1+LFine(2):LFine(3)) = x(1+LFine(2):LFine(3)) + corrF3;
x(1+LFine(3):LFine(4)) = x(1+LFine(3):LFine(4)) + corrF4;
% postsmoothing. Let us do a symmetric v-cycle!! yeah!
top2bottom=0;
x=ArnoldSmoother3(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
    A11_lev{lev},A12_lev{lev},A13_lev{lev},A14_lev{lev},...
    A21_lev{lev},A22_lev{lev},A23_lev{lev},A24_lev{lev},...
    A31_lev{lev},A32_lev{lev},A33_lev{lev},A34_lev{lev},...
    A41_lev{lev},A42_lev{lev},A43_lev{lev},A44_lev{lev},...
   smoothing_steps,is_on_coarser_grid);


end


end


