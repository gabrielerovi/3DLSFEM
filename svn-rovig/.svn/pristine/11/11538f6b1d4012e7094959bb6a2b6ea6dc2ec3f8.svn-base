clear all
close all
clc
%echo on

% data

parameters =parameters();


[mesh,h]=mesh(parameters);

for ii=1:length(mesh)
vv{ii}=0;
end
% print_mesh_L(mesh,vv);
householder=1;
[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L},householder);

[maps,graph]=graph_and_maps(mesh);
L=length(mesh);
maps=ArnoldMonotoneConstraint(mesh,maps,parameters);

AArnoldLocal=ArnoldLocalAssembling(L,mesh,maps,parameters);

[A,b,AFine,P,Ant,bnt,Antbc,AFinenobc,AFinenobcnt,bnt_lev,Complementarity]=create_system(parameters,mesh,h);

Pnt=[];
if(parameters.L>1)
 Pnt=ArnoldPnt(P,mesh);
end



% for ii=1:2
%  
%      for jj=1:4
%           Ant{ii,jj}(mesh{L}.E_remove,:)=0;
%      end
% 
%      for jj=1:4
%      Ant{ii+2,jj}(mesh{L}.N_remove,:)=0;
%      end
%  
%  
%      for jj=mesh{L}.E_remove
%           Ant{ii,ii}(jj,jj)=1;
%      end
%  
%      for jj=mesh{L}.N_remove
%      Ant{ii+2,ii+2}(jj,jj)=1;
%      end
%  
%  end
 

 
%  EcontBC=0;
%  type_of_dof=2;
%  for ii=mesh{L}.E_remove
%  U_xy=[0;0];
%  EcontBC=EcontBC+1;
%  tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,mesh{L}.E_label(EcontBC), 0);
%  coeff = RT_dirichlet_coeff(mesh{L}.E_remove(EcontBC), mesh{L});
%  jj1=ii;
%  jj2=ii+mesh{L}.NE;
%  b(jj1)=tmp(1)/coeff;
%  b(jj2)=tmp(2)/coeff;
%  end
%  
%  NcontBC=0;
%  type_of_dof=1;
%  for ii=mesh{L}.N_remove
%  U_xy=[0;0];
%  NcontBC=NcontBC+1;
%  tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,mesh{L}.N_label(NcontBC), 0);
%  jj1=ii+2*mesh{L}.NE;
%  jj2=ii+2*mesh{L}.NE+mesh{L}.N;
%  b(jj1)=tmp(1);
%  b(jj2)=tmp(2);
%  end
 
 
 
 
 
Constraint = CheckConstraintsNT(mesh,maps,parameters);

activeset=Constraint.WorkingSetE;

big_value=10^10;

c=big_value*ones(2*mesh{L}.NE+2*mesh{L}.N,1);


c(mesh{L}.E_contact)=Constraint.RhsE2;
c(mesh{L}.N_contact+2*mesh{L}.NE)=Constraint.RhsN2;

A=[Antbc{L,1,1} Antbc{L,1,2} Antbc{L,1,3} Antbc{L,1,4};];
for ii=2:4
A=[A;Antbc{L,ii,1} Antbc{L,ii,2} Antbc{L,ii,3} Antbc{L,ii,4};];
end

for ii=1:4
    for jj=1:4
        Aprova{ii,jj}=Antbc{L,ii,jj};
    end
end

% [y,activeset]=activesetResidual(A,bnt,c,activeset);
% sol=M_Normal_Tangent*y;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% p = ContactPressure(mesh,sol,parameters)
% resy=bnt-A*y;
% resy(activeset)=0;
% resy=norm(resy);
% resy=ones(100,1)*resy;

B=speye(length(bnt));
WorkingSet=[];
A=A+Complementarity{L}.Ant_complementarity;
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
bnt=bnt+Complementarity{L}.bnt_complementarity;
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER
% IT IS ALREADY HERE DO NOT ADD IT LATER

[z,lambda,WorkingSet] = ArnoldActiveset2(A,B,bnt,c,WorkingSet);
sol=M_Normal_Tangent*z;
p = ContactPressure(mesh,sol,parameters);
ContactDisplacement(mesh,sol,parameters)

resz=bnt-A*z;
resz(WorkingSet)=0;
resz=norm(resz);
resz=ones(100,1)*resz;

% [w,activesetw]=activesetResidual(A,bnt,c,activesetw);
% sol=M_Normal_Tangent*w;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% 
% 
% B=speye(length(bnt));
% WorkingSet=[];
% [y,lambda,WorkingSet] = ArnoldActiveset2(A,B,bnt,c,WorkingSet);
% resy=bnt-A*y;
% resy(WorkingSet)=0;
% resy=norm(resy);
% resy=ones(100,1)*resy;
% sol=M_Normal_Tangent*y;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% p = ContactPressure(mesh,sol,parameters)



% devi ruotare il sistema di riferimento localeeee

removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE, 2*mesh{L}.NE + mesh{L}.N_remove,2*mesh{L}.NE + mesh{L}.N + mesh{L}.N_remove];
x=bnt;
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes);
x(remove)=0;


is_on_coarser_grid=false;
top2bottom=1;
smoothing_steps=5000;
WorkingSetContact{1}=sparse(mesh{L}.NE*2+mesh{L}.N*2,1);
[x,WorkingSetContact{1}] = ArnoldNestedIteration(Ant,bnt_lev,c,mesh,maps,Pnt);
toll=10^(-12);


% WorkingSet=activesetw';
% mesh{L}.RemoveNT=unique([mesh{L}.RemoveNT,WorkingSet']);
% ConstraintTot=WorkingSet-mesh{L}.NE*2;
% ConstraintNodes=ConstraintTot(ConstraintTot>0);
% ConstraintEdge=WorkingSet(WorkingSet<=mesh{L}.NE*2);
% bnt(ConstraintNodes+mesh{L}.NE*2)=gap_function(mesh{L}.node(ConstraintNodes,1),mesh{L}.node(ConstraintNodes,2));
% bnt(ConstraintEdge)=0;
% for lev=L-1:-1:1
% WorkingSet=WorkingSet(find(WorkingSet<=mesh{lev+1}.NE*2+mesh{lev}.N));
% WorkingSet=WorkingSet-mesh{lev+1}.NE*2+mesh{lev}.NE*2;
% mesh{lev}.RemoveNT=unique([mesh{lev}.RemoveNT,WorkingSet']);
% end

% 


% x=z;
% 
% 
% WorkingSetContact{1}=sparse(2*mesh{L}.NE+2*mesh{L}.N,1);
% WorkingSetContact{1}(WorkingSet)=1;


for mm=1:50

 [x,WorkingSetContact{mm+1},residual(mm),corrF(mm)] = ArnoldVcycleContact3(Complementarity{L}.Ant_complementarity,h,AFinenobcnt,bnt,x,c,WorkingSetContact{mm},parameters,mesh,maps,graph,Pnt,L,L);

% [x] = ArnoldVcycle7(AFinenobcnt,bnt,x,parameters,mesh,maps,graph,Pnt,L,L);
% res=bnt-AFinenobcnt*x;
% res(mesh{L}.RemoveNT)=0;
% residual(mm)=norm(res);
% [y,WorkingSet]=ArnoldSmootherContact7(WorkingSet,A,x,b,c,mesh{L},maps.Patch_Internal_All{L},smoothing_steps,graph{L},top2bottom) 
%[z,WorkingSet]=ArnoldSmootherContact6(L,WorkingSet,x,bnt,mesh,AFinenobcnt,Aprova,c,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);

%[x,WorkingSet,residual(mm)] = ArnoldVcycleContact2(graph,L,L,maps,x,bnt,mesh,Ant,Antbc,Pnt, parameters,c,WorkingSet);


% [x,WorkingSet]=ArnoldSmootherContact6(L,WorkingSet,x,bnt,mesh,Aprova,c,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);
% sol=M_Normal_Tangent*x;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% p = ContactPressure(mesh,sol,parameters)
mm;
% norm(y-x);
energy(mm)=0.5*x'*AFinenobcnt*x-bnt'*x;
norma(mm)=norm(x-z);
[mm,energy(mm),residual(mm),norma(mm)]

if(residual(mm)<toll)
break
end
end

figure
plot(log10(residual))
hold on
plot(log10(resz))
top2bottom=1;
lev=L;
is_on_coarser_grid=0;
smoothing_steps=10;
% figure
%     x=ArnoldSmoother3(graph{lev},top2bottom,maps.EmapGlob2Loc{lev},maps.EmapLoc2Glob{lev},maps.NmapGlob2Loc{lev},maps.NmapLoc2Glob{lev},...
%     maps.Patch_Boundary_Edge{lev}, maps.Patch_Boundary_Node{lev},maps.Patch_Edge{lev}, maps.Patch_Node{lev},x,b,mesh{lev},...
%     A{lev,1,1},A{lev,1,2},A{lev,1,3},A{lev,1,4},...
%     A{lev,2,1},A{lev,2,2},A{lev,2,3},A{lev,2,4},...
%     A{lev,3,1},A{lev,3,2},A{lev,3,3},A{lev,3,4},...
%     A{lev,4,1},A{lev,4,2},A{lev,4,3},A{lev,4,4},...
%     smoothing_steps,is_on_coarser_grid)


















%  for mm=1:5
%      mm
% x = ArnoldVcycleContact1(graph,L,L,maps,x,bnt,mesh,Antbc,Pnt, parameters);
%  sol=M_Normal_Tangent*x;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
%  end
%  end
%  sol=M_Normal_Tangent*x;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');


% [solnt,Constraint]=ArnoldSmootherContact3(x,bnt,mesh,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)
% sol=M_Normal_Tangent*solnt;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');

 for ii=1:4
     for jj=1:4
        Antlev{ii,jj}=   Antbc{L,ii,jj};
     end
 end
 
[solnt,is_the_right_solution,Constraint]=activeset3(Antlev,bnt,Constraint,mesh{L});
sol=M_Normal_Tangent*solnt;
print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
wn{1}=Constraint.WorkingSetN;
we{1}=Constraint.WorkingSetE;
[solnt,is_the_right_solution,Constraint]=activeset3(Ant,bnt,Constraint,mesh{L});
sol=M_Normal_Tangent*solnt;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
wn{2}=Constraint.WorkingSetN;
we{2}=Constraint.WorkingSetE;
[solnt,is_the_right_solution,Constraint]=activeset3(Ant,bnt,Constraint,mesh{L});
sol=M_Normal_Tangent*solnt;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
wn{3}=Constraint.WorkingSetN;
we{3}=Constraint.WorkingSetE;
[solnt,is_the_right_solution,Constraint]=activeset3(Ant,bnt,Constraint,mesh{L});
sol=M_Normal_Tangent*solnt;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
wn{4}=Constraint.WorkingSetN;
we{4}=Constraint.WorkingSetE;
[solnt,is_the_right_solution,Constraint]=activeset3(Ant,bnt,Constraint,mesh{L});
sol=M_Normal_Tangent*solnt;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');





% L=parameters.L;
% N=mesh{L}.N;
% NE=mesh{L}.NE;
% [dirichletN,n_and_or_tN,bool_bcN]= boundary_value_bool(1);


[AFine,b,ST,STT,CT,B,BT,C] = ActiveSetBmatrixCvector(mesh,parameters,AFine,b);



%%% qui proviamo a forzare le bc
% N=mesh{parameters.L}.N;
% NE=mesh{parameters.L}.NE;
% Aforzata=AFine;
% for ii=1:N/2
%     Aforzata(2*NE+N+ii,:)=0;
%     Aforzata(2*NE+N+ii,2*NE+N+ii)=1;
%     b(2*NE+N+ii,1)=0;
% end
% sol=Aforzata\b;
%  print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% 
%  S=B(1:N/2,:)*AFine^(-1)*BT(:,1:N/2);
%  lambda=S\(B(1:N/2,:)*AFine^(-1)*(b-AFine*sol));


removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
x=b;
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes);
x(remove)=0;

[sol1,lambda,WorkingSet] = activeset2(AFine,STT,BT,b,ST,CT,B,C,x);

%[sol1,lambda,Working Set] = activeset(AFine,STT,BT,b,ST,CT,B,C,x);

sol=sol1;
 print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');


  [AFinent,bnt,STnt,STTnt,CTnt,Bnt,BTnt,Cnt] = ActiveSetBmatrixCvectorNormalTangent(mesh,parameters,M_Normal_Tangent,M_Normal_TangentT,AFine,b);
removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
bnt=M_Normal_Tangent*b;
xnt=bnt;
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes);
xnt(remove)=0; 
[solnt,lambdant,WorkingSetnt] = activeset(AFinent,STTnt,BTnt,bnt,STnt,CTnt,Bnt,Cnt,xnt);
 sol=M_Normal_Tangent*solnt;
 print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');

 
removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes); 
x=bnt;
x(remove)=0; 

    Nconstraint=zeros(mesh{L}.N,1);
    Econstraint=zeros(mesh{L}.NE,1);
    
    
% for lev=1:L
%    Ant11{lev}=Ant{lev,1,1};
%    Ant12{lev}=Ant{lev,1,2};
%    Ant13{lev}=Ant{lev,1,3};
%    Ant14{lev}=Ant{lev,1,4};
%    Ant21{lev}=Ant{lev,2,1};
%    Ant22{lev}=Ant{lev,2,2};
%    Ant23{lev}=Ant{lev,2,3};
%    Ant24{lev}=Ant{lev,2,4};
%    
%    Ant31{lev}=Ant{lev,3,1};
%    Ant32{lev}=Ant{lev,3,2};
%    Ant33{lev}=Ant{lev,3,3};
%    Ant34{lev}=Ant{lev,3,4};  
%    
%    Ant41{lev}=Ant{lev,4,1};
%    Ant42{lev}=Ant{lev,4,2};
%    Ant43{lev}=Ant{lev,4,3};
%    Ant44{lev}=Ant{lev,4,4}; 
%    Ass{lev}=[A11{lev} A12{lev};
%              A21{lev} A22{lev}];
%    Asu{lev}=[A13{lev} A14{lev};
%              A23{lev} A24{lev}];
%    Aus{lev}=[A31{lev} A32{lev};
%              A41{lev} A42{lev}];    
%    Auu{lev}=[A33{lev} A34{lev};
%              A43{lev} A44{lev}];      
% end



for mm=1:parameters.MGiter
    
    xold=x;
    lev=L;
    is_on_coarser_grid=0;
    
    top2bottom=0;
    % ArnoldSmootherContact
%     x=ArnoldSmootherContact(graph{lev},top2bottom,maps.EmapGlob2Loc{lev},maps.EmapLoc2Glob{lev},maps.NmapGlob2Loc{lev},maps.NmapLoc2Glob{lev},...
%     maps.Patch_Boundary_Edge{lev}, maps.Patch_Boundary_Node{lev},maps.Patch_Edge{lev}, maps.Patch_Node{lev},x,b,mesh{lev},...
%     A11{lev},A12{lev},A13{lev},A14{lev},...
%     A21{lev},A22{lev},A23{lev},A24{lev},...
%     A31{lev},A32{lev},A33{lev},A34{lev},...
%     A41{lev},A42{lev},A43{lev},A44{lev},...
%     parameters.smoothing_steps,is_on_coarser_grid)

NE=mesh{L}.NE;
N=mesh{L}.N;
LL=[NE,2*NE,2*NE+N,2*NE+2*N ];
A11=AFinent(1:LL(1),1:LL(1) );
A12=AFinent(1:LL(1),LL(1)+1:LL(2) );
A13=AFinent(1:LL(1),LL(2)+1:LL(3) );
A14=AFinent(1:LL(1),LL(3)+1:LL(4) );

A21=AFinent(LL(1)+1:LL(2),1:LL(1) );
A22=AFinent(LL(1)+1:LL(2),LL(1)+1:LL(2)  );
A23=AFinent(LL(1)+1:LL(2),LL(2)+1:LL(3) );
A24=AFinent(LL(1)+1:LL(2),LL(3)+1:LL(4)  );

A31=AFinent(LL(2)+1:LL(3),1:LL(1) );
A32=AFinent(LL(2)+1:LL(3),LL(1)+1:LL(2)  );
A33=AFinent(LL(2)+1:LL(3),LL(2)+1:LL(3) );
A34=AFinent(LL(2)+1:LL(3),LL(3)+1:LL(4)  );

A41=AFinent(LL(3)+1:LL(4) ,1:LL(1) );
A42=AFinent(LL(3)+1:LL(4) ,LL(1)+1:LL(2)  );
A43=AFinent(LL(3)+1:LL(4) ,LL(2)+1:LL(3) );
A44=AFinent(LL(3)+1:LL(4) ,LL(3)+1:LL(4)  );



% [x,WorkingSetNormal_N,WorkingSetNormal_E]=ArnoldSmootherContactNormalTangent(graph{lev},top2bottom,maps.EmapGlob2Loc{lev},maps.EmapLoc2Glob{lev},maps.NmapGlob2Loc{lev},maps.NmapLoc2Glob{lev},...
%     maps.Patch_Boundary_Edge{lev}, maps.Patch_Boundary_Node{lev},maps.Patch_Edge{lev}, maps.Patch_Node{lev},x,b,mesh{lev},...
%     A11,A12,A13,A14,...
%     A21,A22,A23,A24,...
%     A31,A32,A33,A34,...
%     A41,A42,A43,A44,...
%     parameters.smoothing_steps,is_on_coarser_grid);

    x = ArnoldVcycleContact(graph,L,L,maps.EmapGlob2Loc,maps.EmapLoc2Glob,maps.NmapGlob2Loc,maps.NmapLoc2Glob,...
    maps.Patch_Boundary_Edge, maps.Patch_Boundary_Node,maps.Patch_Edge, maps.Patch_Node,x,b,mesh,...
    Ant,Pnt, parameters);




contact_integral=ContactIntegral(mesh,x,parameters);
  % print_displacement_solution(mesh,x(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',x(1+2*mesh{L}.NE+mesh{L}.N:end)');
  y=x;

  [J]=LSfunctional(x',parameters.force1,parameters.force2,parameters.alpha,parameters.beta,parameters.C_eq,parameters.C_const,parameters.C_asym,mesh,parameters.qrule);
  residual(mm)=norm(AFine*y-b);
  if(mm>1)
      rate(mm-1)=residual(mm)/residual(mm-1);
  end
  energy(mm)=sum(J);
%   norm(sol-M_Normal_Tangent*x)
end

%[WorkingSetNormal_E,WorkingSetNormal_N]=WorkingSetsNormalTangent(x,mesh{L})


x=M_Normal_Tangent*x;
print_displacement_solution(mesh,x(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',x(1+2*mesh{L}.NE+mesh{L}.N:end)');
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
% Ant=AFine*M_Normal_Tangent;
% bnt=b;
 Bnt=B*M_Normal_TangentT;
 BTnt=Bnt';
 STnt=ST*M_Normal_TangentT;
 STTnt=STnt'; 
 Cnt=C;
% BTnt=Bnt';
% STnt=ST*M_Normal_Tangent;
% STTnt=STnt';
[solnt,lambdant] = activeset(Ant,STTnt,BTnt,bnt,STnt,CTnt,Bnt,Cnt,xnt);
 sol=M_Normal_Tangent*solnt;
 print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% 
% v=[15,16];

% for ii=v
% AFine(ii,:)=0;
% AFine(ii,ii)=1;
% b(ii)=0;
% end


% sol=AFine\b;


% z=[1 2 3 4 5 6 7 8 9 10 15 16 17 18 11 12 13 14];
% AFine(:,z);
% % for ii=v
% % tryA(ii,:)=0;
% % tryA(ii,ii)=1;
% % b(ii)=0;
% % end
% % trysol=tryA\b;
% % sol=M_Normal_Tangent*trysol;
% % w=[13,14,17,18];
% % for ii=w
% % tryA(ii,:)=0;
% % tryA(ii,ii)=1;
% % b(ii)=0;
% % end
% % b([17,18])=-0.05;
% % for nn=1:mesh{L}.N
% %     if(mesh{L}.N_bc(nn)>0 && dirichletN(mesh{L}.N_bc(nn),3)==1)
% %     tryA(2*NE+nn,:)=0;
% %     tryA(2*NE+nn,2*NE+nn)=1;
% %     b(2*NE+nn)=0;
% %     end
% % end


% trysol=tryA\b;
% sol=M_Normal_Tangent*trysol;
%sol=full(AFine\b);
% print_bdf(mesh)
% meshes_write(mesh)

% for i=1:mesh{parameters.L}.N
%     mesh{parameters.L}.normal_node{i}=[0;-1];
% end


% [B,BT,C] = ActiveSetBmatrixCvectorDisplacementFormulation(mesh,parameters)
% M=[AFine,-BT;B, sparse(length(C),length(C))];
% F=[b;C];
% sol1=M\F;
% sol=sol1(1:mesh{parameters.L}.N*2);
% lambda=sol1(1+mesh{parameters.L}.N*2:end);
% print_displacement_solution(mesh,sol(1:mesh{L}.N)',sol(1+mesh{L}.N:end)');

removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
x=b;
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes);
x(remove)=0; 

% ST, STT, CT equality constraints and its transpose (for frictionless) and
% the corresponding rhs
% B, BT, C inequality constraints and their transpose and the corresponding rhs




[AFine,b,ST,STT,CT,B,BT,C] = ActiveSetBmatrixCvector(mesh,parameters,AFine,b);

Ant=AFine*M_Normal_Tangent;
STnt=ST*M_Normal_Tangent;
STTnt=STnt';
Bnt= B*M_Normal_Tangent;%([mesh{parameters.L}.NE*2+mesh{parameters.L}.N_contact, mesh{parameters.L}.E_contact],[mesh{parameters.L}.NE*2+mesh{parameters.L}.N_contact, mesh{parameters.L}.E_contact]);
BTnt=Bnt';
[solnt,lambdant] = activeset(Ant,STTnt,BTnt,b,STnt,CT,Bnt,C,x);
sol=M_Normal_Tangent*solnt;
print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');

BT=[];
C=[];
[sol,lambda] = activeset(AFine,STT,BT,b,ST,CT,B,C,x);
%  
%  
%  removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
% x=b;
% tmp=1:(mesh{L}.NE*2);
% remove=setdiff(tmp,removes);
% x(remove)=0;
% M=[A,-BT;B, sparse(length(C),length(C))];
% F=[b;C];
% M=[AFine -B';B zeros(length(C),length(C))];
% F=[b;C];
% sol1=M\F;
% sol=sol1(1:mesh{L}.NE*2+mesh{L}.N*2);
% lambda=sol1(1+mesh{L}.NE*2+mesh{L}.N*2:end);
% [sol,lambda]=activeset(A,BT,b,B,C,x);
% 
% [sol,lambda]=activeset(AFine,B',b,B,C,x);




%print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% print_displacement_solution(mesh,sol(1:mesh{L}.N)',sol(1+mesh{L}.N:end)');


% for nn=1:mesh{1}.N
%   if(mesh{1}.node(nn,2)==0)
%       mesh{1}.N_bc(nn)=1;
%   end
% end
% 
% for ii=1:2
% mesh{1}.normal_node{ii}=[0;-1];
% end


% parameters.input_name='LSelasticityAsymmetric';
% [A,b,AFine,P]=create_system(parameters,mesh,h);
% solasym=full(AFine\b);
% norm(solasym(1+mesh{L}.NE*2)-soldisp)
% norm(solasym(1+mesh{L}.NE*2)-soldisp)/mesh{L}.N
% max(abs(solasym(1+mesh{L}.NE*2)-soldisp))
% 
%  sol=solasym(1+mesh{L}.NE*2)-soldisp;
% print_displacement_solution(mesh,sol(1:mesh{L}.N)',sol(1+mesh{L}.N:end)');
% if(strcmp(parameters.input_name,'DispElasticity'))
% print_displacement_solution(mesh,sol(1:mesh{L}.N)',sol(1+mesh{L}.N:end)');
% else
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
% end
% x=zeros(mesh{L}.NE*parameters.E_components + mesh{L}.N*parameters.N_components,1);






removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
x=b;
tmp=1:(mesh{L}.NE*2);
remove=setdiff(tmp,removes);
x(remove)=0; 
if(strcmp(parameters.input_name,'LSstressblock'))
    remove=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
    bs=b(1:mesh{L}.NE*2); 
    bu=[];
for lev=1:L
   A11{lev}=A{lev,1,1};
   A12{lev}=A{lev,1,2};
   A21{lev}=A{lev,2,1};
   A22{lev}=A{lev,2,2};
end

elseif(strcmp(parameters.input_name,'DispElasticity'))
    remove=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE];
    bu=b; 
    bs=[];
for lev=1:L
   A11{lev}=A{lev,1,1};
   A12{lev}=A{lev,1,2};
   A21{lev}=A{lev,2,1};
   A22{lev}=A{lev,2,2};
   
   Auu{lev}=[ A11{lev} A12{lev}; 
              A21{lev} A22{lev}];
          
end
elseif(strcmp(parameters.input_name,'LSelasticity')||strcmp(parameters.input_name,'LSelasticityAsymmetric'))
    remove=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE, mesh{L}.N_remove + 2 * mesh{L}.NE, mesh{L}.N_remove + 2 * mesh{L}.NE+  mesh{L}.N ];
    bs=b(1:mesh{L}.NE*2); 
    bu=b(1+mesh{L}.NE*2:end);
for lev=1:L
   A11{lev}=A{lev,1,1};
   A12{lev}=A{lev,1,2};
   A13{lev}=A{lev,1,3};
   A14{lev}=A{lev,1,4};
   A21{lev}=A{lev,2,1};
   A22{lev}=A{lev,2,2};
   A23{lev}=A{lev,2,3};
   A24{lev}=A{lev,2,4};
   
   A31{lev}=A{lev,3,1};
   A32{lev}=A{lev,3,2};
   A33{lev}=A{lev,3,3};
   A34{lev}=A{lev,3,4};  
   
   A41{lev}=A{lev,4,1};
   A42{lev}=A{lev,4,2};
   A43{lev}=A{lev,4,3};
   A44{lev}=A{lev,4,4}; 
   Ass{lev}=[A11{lev} A12{lev};
             A21{lev} A22{lev}];
   Asu{lev}=[A13{lev} A14{lev};
             A23{lev} A24{lev}];
   Aus{lev}=[A31{lev} A32{lev};
             A41{lev} A42{lev}];    
   Auu{lev}=[A33{lev} A34{lev};
             A43{lev} A44{lev}];      
end
end





% DISPLACEMENT FORMULATION

if(strcmp(parameters.input_name,'DispElasticity'))
x=b;
remove=[mesh{L}.N_dirichlet;mesh{L}.N_dirichlet+mesh{L}.N];
tmp=1:(mesh{L}.N*2);
remove=setdiff(tmp,remove);
x(remove)=0;

disp=x;

resuu = bu -  Auu{L} * disp ;
    
for mm=1:parameters.MGiter

%     qui cambia ATTENZIONE. 
     lev=L;
     top2bottom=1;
     is_on_coarser_grid=false;

    normresuu_old=norm(resuu);
    correctionuu=zeros(2*mesh{L}.N,1); 
           
    
    correctionuu =v_cycle(L,Auu,resuu,correctionuu,P.P1CtoP1Fuu,parameters.smoothing_steps,parameters.GS_is_symmetric,mesh);
    
    disp=disp+correctionuu;
    
    
   % [disp]=gauss_seidel_contact_displacement(Auu{L},b,disp,parameters.smoothing_steps,mesh)
%print_displacement_solution(mesh,disp(1:mesh{L}.N)',disp(1+mesh{L}.N:end)')

    resuu = bu - Auu{L} * disp ; 
    
    normresuu_new=norm(resuu);
    
    
    %rate(mm)=normresuu_new/normresuu_old
    residual(mm)=normresuu_new;
    
 if(normresuu_new<parameters.toll)
     break;
 end

end

x=disp;
print_displacement_solution(mesh,x(1:mesh{L}.N)',x(1+mesh{L}.N:end)')

end















x=b;
tmp=1:(mesh{L}.N*2+mesh{L}.NE*2);
remove=setdiff(tmp,remove);
x(remove)=0;
if(parameters.contact==0)
stress=x(1:mesh{L}.NE*2);
disp=x(1+mesh{L}.NE*2:end);

    resss = bs - (Ass{L} * stress  +  Asu{L} * disp );
    resuu = bu - (Aus{L} * stress  +  Auu{L} * disp );

    contss=0;
    cont=0;
for mm=1:parameters.MGiter

%     qui cambia ATTENZIONE. 
     lev=L;
     top2bottom=1;
     is_on_coarser_grid=false;

% QUI FACCIAMO STAGGERED CON LA CONDENSAZIONE DELLE BC AL BORDO
%    correction1= ass\resss;
%    stress=stress+correction1;
%    resss = bss - (ass * stress  +  asu * disp );
%    resuu = buu - (aus * stress  +  auu * disp );
%     
%    correction2= auu\resuu; 
%    disp=disp+correction2;
%    resss = bss - (ass * stress  +  asu * disp );
%    resuu = buu - (aus * stress  +  auu * disp );
%    res=norm([resss;resuu]);

  
%     stress=stress+correction(1:mesh{L}.NE*2);
%     disp=disp+correction(1+mesh{L}.NE*2:end);
%     resss = bs - (Ass{L} * stress  +  Asu{L} * disp );
%     resuu = bu - (Aus{L} * stress  +  Auu{L} * disp );
    
    


  
   
%     residualss=parameters.toll+1;
%     mgss=1;
%     correctionss=zeros(mesh{L}.NE*2,1);
%     while(residualss>parameters.toll && mgss<parameters.MGiterSS)
%     correctionss = ArnoldVcycleStress(graph,L,lev,maps,correctionss,resss,mesh,...
%     A11,A12,A21,A22,...
%     parameters.smoothing_steps,P.RTCtoRTF,P.P1CtoP1F,parameters.LineSearch);
%     residualss=norm(resss-Ass{L}*correctionss)
%     mgss=mgss+1;
%     contss=contss+1;
%     cont=cont+1;
%     end
%     stress=stress+correctionss;
% 
%     resss = bs - (Ass{L} * stress  +  Asu{L} * disp );
%     resuu = bu - (Aus{L} * stress  +  Auu{L} * disp );
% 
%  
%     residualuu=parameters.toll+1;
%     mguu=1;
%     correctionuu=zeros(2*mesh{L}.N,1); 
%     
%     while(residualuu>parameters.toll && mguu<parameters.MGiterUU)
%     residualuuold=residualuu;
%     correctionuu =v_cycle(L,Auu,resuu,correctionuu,P.P1CtoP1Fuu,parameters.smoothing_steps,parameters.GS_is_symmetric,mesh);
%     residualuu=norm(resuu-Auu{L}*correctionuu) 
%     mguu=mguu+1;
%     residualuuuuuu(mguu)=residualuu;
%     cont=cont+1;
%     end
% %     correctionuu=Auu{L}\resuu;
%     disp=disp+correctionuu;
%     resss = bs - (Ass{L} * stress  +  Asu{L} * disp );
%     resuu = bu - (Aus{L} * stress  +  Auu{L} * disp );   
%     residual(mm)=norm([resss;resuu]);
%     if(mm>1)
%     rate(mm-1)=residual(mm)/residual(mm-1);
%     end
%  if(residual(mm)<parameters.toll)
%      break;
%  end
    
     lev=L;
     top2bottom=1;
     is_on_coarser_grid=false;
LineSearch=false;

    x = ArnoldVcycle3(graph,L,L,maps.EmapGlob2Loc,maps.EmapLoc2Glob,maps.NmapGlob2Loc,maps.NmapLoc2Glob,...
    maps.Patch_Boundary_Edge, maps.Patch_Boundary_Node,maps.Patch_Edge, maps.Patch_Node,x,b,mesh,...
    A11,A12,A13,A14,...
    A21,A22,A23,A24,...
    A31,A32,A33,A34,...
    A41,A42,A43,A44,...
   parameters.smoothing_steps,P.RTCtoRTF,P.P1CtoP1F,parameters.LineSearch);
res=norm(b-AFine*x)
%res=norm([resss;resuu])   
residual(mm,1)=res;
if(mm>1)
    rate(mm-1)=residual(mm,1)/residual(mm-1,1);
end
if(res<parameters.toll)
    break
end
if(mm>6 && res>resold)
    break;
end
resold=res;

 
% J_discretized=LSfunctional(x',force1,force2,alpha,beta, coeff_equilibrium,coeff_constitutive,coeff_symmetry,mesh,qrule);
% 
% J_approx(:,mm)=J_discretized;
% J_diff(:,mm)=abs(J_discretized-J_exact)./J_exact;

end




print_displacement_solution(mesh,x(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',x(1+2*mesh{L}.NE+mesh{L}.N:end)');
% print_displacement_solution(mesh,sol(1:mesh{L}.N)',sol(1+mesh{L}.N:end)');


end




if(parameters.contact==1)
x=b;
x(remove)=0; 

    Nconstraint=zeros(mesh{L}.N,1);
    Econstraint=zeros(mesh{L}.NE,1);
for mm=1:parameters.MGiter
    
    xold=x;
    lev=L;
    is_on_coarser_grid=0;
    
    top2bottom=0;
    % ArnoldSmootherContact
%     x=ArnoldSmootherContact(graph{lev},top2bottom,maps.EmapGlob2Loc{lev},maps.EmapLoc2Glob{lev},maps.NmapGlob2Loc{lev},maps.NmapLoc2Glob{lev},...
%     maps.Patch_Boundary_Edge{lev}, maps.Patch_Boundary_Node{lev},maps.Patch_Edge{lev}, maps.Patch_Node{lev},x,b,mesh{lev},...
%     A11{lev},A12{lev},A13{lev},A14{lev},...
%     A21{lev},A22{lev},A23{lev},A24{lev},...
%     A31{lev},A32{lev},A33{lev},A34{lev},...
%     A41{lev},A42{lev},A43{lev},A44{lev},...
%     parameters.smoothing_steps,is_on_coarser_grid)

NE=mesh{L}.NE;
N=mesh{L}.N;
LL=[NE,2*NE,2*NE+N,2*NE+2*N ];
A11=Ant(1:LL(1),1:LL(1) );
A12=Ant(1:LL(1),1:LL(2) );
A13=Ant(1:LL(1),1:LL(3) );
A14=Ant(1:LL(1),1:LL(4) );

A21=Ant(1:LL(2),1:LL(1) );
A22=Ant(1:LL(2),1:LL(2) );
A23=Ant(1:LL(2),1:LL(3) );
A24=Ant(1:LL(2),1:LL(4) );

A31=Ant(1:LL(3),1:LL(1) );
A32=Ant(1:LL(3),1:LL(2) );
A33=Ant(1:LL(3),1:LL(3) );
A34=Ant(1:LL(3),1:LL(4) );

A41=Ant(1:LL(4),1:LL(1) );
A42=Ant(1:LL(4),1:LL(2) );
A43=Ant(1:LL(4),1:LL(3) );
A44=Ant(1:LL(4),1:LL(4) );



x=ArnoldSmootherContactNormalTangent(graph{lev},top2bottom,maps.EmapGlob2Loc{lev},maps.EmapLoc2Glob{lev},maps.NmapGlob2Loc{lev},maps.NmapLoc2Glob{lev},...
    maps.Patch_Boundary_Edge{lev}, maps.Patch_Boundary_Node{lev},maps.Patch_Edge{lev}, maps.Patch_Node{lev},x,b,mesh{lev},...
    A11,A12,A13,A14,...
    A21,A22,A23,A24,...
    A31,A32,A33,A34,...
    A41,A42,A43,A44,...
    parameters.smoothing_steps,is_on_coarser_grid)


   contact_integral=ContactIntegral(mesh,x,parameters);
  % print_displacement_solution(mesh,x(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',x(1+2*mesh{L}.NE+mesh{L}.N:end)');
  y=x;

  [J]=LSfunctional(x',parameters.force1,parameters.force2,parameters.alpha,parameters.beta,parameters.C_eq,parameters.C_const,parameters.C_asym,mesh,parameters.qrule);
  residual(mm)=norm(AFine*y-b);
  if(mm>1)
      rate(mm-1)=residual(mm)/residual(mm-1);
  end
  energy(mm)=sum(J);
  norm(xold-x)
end
print_displacement_solution(mesh,x(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',x(1+2*mesh{L}.NE+mesh{L}.N:end)');
p = ContactPressure(mesh,x,parameters)
figure
plot(energy)
% y=x;
% y=[y;zeros(2*mesh{L}.N,1)];
% [J]=LSfunctional(y,force1,force2,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,mesh,qrule)



p = ContactPressure(mesh,sol,parameters)

%print_vector_solution(mesh,x)
%     component=4;
%     if(component==3)
%     vector=1:mesh{L}.N;
%     vector_extra=1+mesh{L}.N:2*mesh{L}.N;
%     resuu(vector)=resuu(vector)-A34{L}*disp(vector_extra);
%     
%     else
%     vector_extra=1:mesh{L}.N;
%     vector=1+mesh{L}.N:2*mesh{L}.N;
%     resuu(vector)=resuu(vector)-A43{L}*disp(vector_extra);
%     end
%     
%     resuu(vector)=ones(length(resuu(vector)),1);
%     P.P1CtoP1F{1}(abs(P.P1CtoP1F{1})<parameters.toll)=0
%     
%     
%     while(residualuu>parameters.toll && mguu<parameters.MGiterUU)
%     residualuuold=residualuu;
%     if(component==3)
%     correctionuu(vector) =v_cycle(L,A33,resuu(vector),correctionuu(vector),P.P1CtoP1F,parameters.smoothing_steps,parameters.GS_is_symmetric,mesh);
%     else
%     correctionuu(vector) =v_cycle(L,A44,resuu(vector),correctionuu(vector),P.P1CtoP1F,parameters.smoothing_steps,parameters.GS_is_symmetric,mesh);
%     end
%     if(component==3)
%     residualuu=norm(resuu(vector)-A33{L}*correctionuu(vector));
%     else
%     residualuu=norm(resuu(vector)-A44{L}*correctionuu(vector));
%     end
%     residualuu
%         mguu=mguu+1;
%     rate(mguu)=residualuu/residualuuold
%     residual222(mguu)=residualuu
%     end
end
