function x = ArnoldVcycle4(graph,maxlev,lev,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,A12_lev,A21_lev,A22_lev,...
    smoothing_steps,RTCtoRTF,P1CtoP1F,LineSearch)

 % fine and coarse
 F=lev;
 C=lev-1;
 
 % N and NE on the fine mesh
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 
 E_label=mesh{F}.E_label;
 % fine level lenghts
 LFine=[ NEF ; 2 * NEF; 2 * NEF + NF ; 2 * NEF + 2 * NF ];
 
 
  E_remove=mesh{lev}.E_remove;
  E_dirichlet=mesh{lev}.E_dirichlet;

  % create matrix at the rough level
 A=[A11_lev{lev} A12_lev{lev};
    A21_lev{lev} A22_lev{lev} ];
 
  % consider the stress and displacement unknowns
  sigma1=x(1:LFine(1));
  sigma2=x(1+LFine(1):LFine(2));

  
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
 
 LCoarse=[NEC; 2 * NEC];  
% presmoothing
top2bottom=1;    
x=ArnoldSmoother4(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
   A11_lev{lev},A12_lev{lev},A21_lev{lev},A22_lev{lev},...
   smoothing_steps,is_on_coarser_grid);



resF= b - A* x;
resC(1:LCoarse(1),1)         =    RTCtoRTF{C}' * resF(1:LFine(1));
resC(1+LCoarse(1): LCoarse(2),1)=    RTCtoRTF{C}' * resF(1+LFine(1):LFine(2));


E_removeC=mesh{C}.E_remove;
E_removeC=[E_removeC,E_removeC+mesh{C}.NE];
resC(E_removeC,1)=0;

correction=zeros(LCoarse(2),1);
% V-Cycle 
correction = ArnoldVcycle4(graph,maxlev,lev-1,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,correction,resC,mesh,...
    A11_lev,A12_lev,A21_lev,A22_lev,...
    smoothing_steps,RTCtoRTF,P1CtoP1F,LineSearch);


corrF1=RTCtoRTF{C} * correction(1:LCoarse(1));
corrF2=RTCtoRTF{C} * correction(1+LCoarse(1):LCoarse(2));

E_removeF=mesh{F}.E_remove;
corrF1(E_removeF)=0;
corrF2(E_removeF)=0;


d=doubleLineSearch (A11_lev{lev},A12_lev{lev},A21_lev{lev},A22_lev{lev},corrF1,corrF2,resF(1:mesh{lev}.NE), resF(mesh{lev}.NE+1:end));
  
  
x(1:LFine(1))          = x(1:LFine(1))          + d(1) * corrF1;
x(1+LFine(1):LFine(2)) = x(1+LFine(1):LFine(2)) + d(2) * corrF2;

% postsmoothing. Let us do a symmetric v-cycle!! yeah!
top2bottom=0;
x=ArnoldSmoother4(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
   A11_lev{lev},A12_lev{lev},A21_lev{lev},A22_lev{lev},...
   smoothing_steps,is_on_coarser_grid);


end


end


