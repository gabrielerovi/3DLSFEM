function x = ArnoldVcycleContact1(graph,maxlev,L,maps,x,bnt,mesh,Ant,P, parameters)

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

A=[ Ant{F,1,1}, Ant{F,1,2}, Ant{F,1,3}, Ant{F,1,4};
    Ant{F,2,1}, Ant{F,2,2}, Ant{F,2,3}, Ant{F,2,4};
    Ant{F,3,1}, Ant{F,3,2}, Ant{F,3,3}, Ant{F,3,4};
    Ant{F,4,1}, Ant{F,4,2}, Ant{F,4,3}, Ant{F,4,4}];

if(F==1)
   x=A\bnt; 
else
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
E_contact_tangentF=mesh{F}.E_contact;  

removeF=[E_remove,E_remove+NEF,N_remove+2*NEF,N_remove+NF+2*NEF ];
removeF=[removeF,E_contact_tangentF,E_contact_tangentF+NEF];
removeF=unique(removeF);
is_on_coarser_grid=0;
top2bottom=1;    



 for ii=1:4
     for jj=1:4
        Antlev{ii,jj}=   Ant{F,ii,jj};
     end
 end
 
[x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},Antlev,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);

 
% compute fine residual on the fine level
resF= bnt - A * x;  

resF(removeF)=0;


resC=P{C}'*resF;

E_removeC=mesh{C}.E_remove;
E_removeC1=E_removeC;
E_removeC2=E_removeC1+NEC;

N_removeC=mesh{C}.N_remove;
N_removeC1 = N_removeC + 2 * NEC;
N_removeC2 = N_removeC1+ NC;

E_contact_tangentC=mesh{C}.E_contact;  

removeC=[E_removeC1,E_removeC2,N_removeC1,N_removeC2];

removeC=[removeC,E_contact_tangentC,E_contact_tangentC+NEC];
removeC=unique(removeC);

resC(removeC,1)=0;

%%%%% qui fai Vcycle
corrC=sparse(2*NEC+2*NC,1);
LineSearch=0;

corrC = ArnoldVcycleContact1(graph,maxlev,C,maps,corrC,resC,mesh,Ant,P, parameters);

corrF=P{C}*corrC;

corrF(removeF)=0;

x=x+corrF;


is_on_coarser_grid=0;
[x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},Antlev,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);

end


end


