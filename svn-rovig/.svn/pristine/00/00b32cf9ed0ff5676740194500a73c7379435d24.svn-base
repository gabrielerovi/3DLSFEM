function x = ArnoldVcycleLinearElasticNormalTangent(graph,maxlev,lev,maps,x,bnt,mesh,...
    Ant,smoothing_steps,Pnt,LineSearch)

 % fine and coarse
 F=lev;
 C=lev-1;
 
 % N and NE on the fine mesh
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 
 E_label=mesh{F}.E_label;
 N_label=mesh{F}.N_label;
 % fine level lenghts
 LF=[1; NEF ; 2 * NEF; 2 * NEF + NF ; 2 * NEF + 2 * NF ];
 
 
  E_remove=mesh{lev}.E_remove;
  E_dirichlet=mesh{lev}.E_dirichlet;

  N_remove=mesh{lev}.N_remove;
  N_dirichlet=mesh{lev}.N_dirichlet;
  
  % create matrix at the rough level
for mm=1:4
    for nn=1:4
        AntF{mm,nn}=Ant{F}(LF(mm):LF(mm+1),LF(nn):LF(nn+1));
        A(LF(mm):LF(mm+1),LF(nn):LF(nn+1))=AntF{mm,nn};
    end
end


  
  if(maxlev==lev)
    is_on_coarser_grid=false;
else
    is_on_coarser_grid=true;
end




if(lev==1)

    x = A\ bnt;

else
 % coarse level lengths
 NEC=mesh{C}.NE;
 NC=mesh{C}.N;
 
 LC=[NEC; 2 * NEC;2 * NEC + NC ; 2 * NEC + 2 * NC ];  
% presmoothing
top2bottom=1;    







[x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},AntF,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);




resF= bnt - A* x;
resC   =  Pnt{C}' * resF;


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

corrC=sparse(LC(end),1);
% V-Cycle 
%corrC = ArnoldVcycleLinearElasticNormalTangent(graph,maxlev,lev-1,maps,corrC,resC,mesh,Ant,smoothing_steps,Pnt,LineSearch);


corrF   =  Pnt{C} * corrC;



E_removeF=mesh{F}.E_remove;
E_removeF1=E_removeF;
E_removeF2=E_removeF1+NEF;

N_removeF=mesh{F}.N_remove;
N_removeF1 = N_removeF + 2 * NEF;
N_removeF2 = N_removeF1+ NF;


removeCTOT=[E_removeF1,E_removeF2,N_removeF1,N_removeF2];


corrF(removeCTOT)=0;










if(LineSearch==true)
corrF=[corrF1;corrF2;corrF3;corrF4];
Ac=A*corrF;
alpha = (resF'*corrF) /(Ac'*corrF);
else
alpha=1;
end


x=x+ alpha * corrF;

% postsmoothing. Let us do a symmetric v-cycle!! yeah!
top2bottom=0;
[x]=ArnoldSmootherContact5(F,x,bnt,mesh{F},AntF,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);



end


end


