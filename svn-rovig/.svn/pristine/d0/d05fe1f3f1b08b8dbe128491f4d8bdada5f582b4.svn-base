function x = ArnoldVcycle2(graph,maxlev,lev,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,smoothing_steps,eta,RTCtoRTF,P1CtoP1F,LineSearch)

 % fine and coarse
 F=lev;
 C=lev-1;
 
 % N and NE on the fine mesh
 NE=mesh{F}.NE;
 N=mesh{F}.N;
 E_label=mesh{F}.E_label;
 N_label=mesh{F}.N_label;
 % fine level lenghts
 LFine=[mesh{F}.NE];
 
 
  E_remove=mesh{lev}.E_remove;
  N_remove=mesh{lev}.N_remove;    
  E_dirichlet=mesh{lev}.E_dirichlet;
  N_dirichlet=mesh{lev}.N_dirichlet;

  % create matrix at the rough level
  A=[A11_lev{lev}];
 
  % consider the stress and displacement unknowns
  sigma1=x(1:LFine(1));

if(maxlev==lev)
    is_on_coarser_grid=false;
else
    is_on_coarser_grid=true;
end



if(lev==1)

    x = A\ b;

else
 % coarse level lengths
 LCoarse=[mesh{C}.NE ];  
% presmoothing
top2bottom=1;    
x=ArnoldSmoother2(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
   A11_lev{lev},smoothing_steps,eta,is_on_coarser_grid);



resF= b - A* x;

resC(1:LCoarse(1), 1)=            RTCtoRTF{C}' * resF(1:LFine(1));

E_removeC=mesh{C}.E_remove;
resC(E_removeC)=0;

correction=zeros(LCoarse(1),1);
% V-Cycle 
correction = ArnoldVcycle2(graph,maxlev,lev-1,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,correction,resC,mesh,...
    A11_lev,smoothing_steps,eta,RTCtoRTF,P1CtoP1F,LineSearch);


corrF1=RTCtoRTF{C} * correction(1:LCoarse(1));


E_removeF=mesh{F}.E_remove;N_removeF=mesh{F}.N_remove;
corrF1(E_removeF)=0;


if(LineSearch==true)
corrF=[corrF1];
Ac=A*corrF;
alpha = (resF'*corrF) /(Ac'*corrF);
else
alpha=1;
end

x(1:LFine(1))          = x(1:LFine(1))          + alpha * corrF1;



% postsmoothing. Let us do a symmetric v-cycle!! yeah!
top2bottom=0;
x=ArnoldSmoother2(graph{lev},top2bottom,EmapGlob2Loc{lev},EmapLoc2Glob{lev},NmapGlob2Loc{lev},NmapLoc2Glob{lev},...
    Patch_Boundary_Edge{lev},Patch_Boundary_Node{lev}, Patch_Edge{lev},Patch_Node{lev},x,b,mesh{lev},...
   A11_lev{lev},smoothing_steps,eta,is_on_coarser_grid);


end


end


