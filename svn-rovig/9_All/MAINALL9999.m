clear all
close all
clc
%echo on

% data

parameters =parameters();


[mesh,h]=mesh(parameters);
L=length(mesh);
    
[maps,graph]=graph_and_maps(mesh);
maps=ArnoldMonotoneConstraint(mesh,maps,parameters);

for lev=1:length(mesh)
%     subplot(2,length(mesh)/2,lev);
    vv{lev}=0;
% plotsquare(lev) = showmesh(mesh{lev}.node,mesh{lev}.elem)
end
print_mesh_L(mesh,vv);


[nodes_colors,colors_cardinality]=UTOPIAcoloring(mesh,maps,0);
print_mesh_colors(mesh,nodes_colors);
vertices_per_processor= UTOPIA_processsors(mesh,parameters,6);
[AFinenobcnt,Complementarity,bnt,bnt1,Pnt]=create_system_contact(parameters,mesh,h,maps.Patch_Node_Monotone);

 Constraint = CheckConstraintsNT(mesh,parameters);


% activeset=Constraint.WorkingSetE;
big_value=10^10;
c=big_value*ones(2*mesh{L}.NE+2*mesh{L}.N,1);
c(mesh{L}.E_contact)=Constraint.RhsE2;
c(mesh{L}.N_contact+2*mesh{L}.NE)=Constraint.RhsN2;

% B=speye(length(bnt));
% WorkingSet=[];

% IT IS ALREADY HERE DO NOT ADD IT LATER



% devi ruotare il sistema di riferimento localeeee

% removes=[mesh{L}.E_remove,mesh{L}.E_remove+mesh{L}.NE, 2*mesh{L}.NE + mesh{L}.N_remove,2*mesh{L}.NE + mesh{L}.N + mesh{L}.N_remove];
% x=bnt;
% tmp=1:(mesh{L}.NE*2);
% remove=setdiff(tmp,removes);
% x(remove)=0;


for lev=1:L
    mesh_parameters{lev}.RemoveNT=mesh{lev}.RemoveNT;
    mesh_parameters{lev}.NE=mesh{lev}.NE;
    mesh_parameters{lev}.N=mesh{lev}.N;
    mesh_parameters{lev}.E_contact=mesh{lev}.E_contact;
    mesh_parameters{lev}.N_contact=mesh{lev}.N_contact;
    mesh_parameters{lev}.RemoveNT=mesh{lev}.RemoveNT;
    mesh_parameters{lev}.Patch_Internal_All=maps.Patch_Internal_All{lev};
    if(lev<L)
    mesh_parameters{lev}.Patch_Edge_Monotone=maps.Patch_Edge_Monotone{lev};
    mesh_parameters{lev}.Patch_Node_Monotone=maps.Patch_Node_Monotone{lev};
    end
end


% clearvars mesh
% clearvars maps

[x,WorkingSet] = ArnoldNestedIteration(AFinenobcnt,Complementarity,h,bnt1,c,mesh_parameters,Pnt);


WorkingSetContact{1}=WorkingSet;%sparse(mesh{L}.NE*2+mesh{L}.N*2,1);
% load('x.mat');
% load('activesetlambda1e50mu1MIN.mat')
% load('isitlinearlambda1e50mu1MIN.mat')
% WorkingSetContact{1}=sparse(length(WorkingSetContact{1}),1);
% WorkingSetContact{1}(activeset)=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRY UTOPIA SMOOTHER PARALLEL VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                BEGIN                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% given the processor pp, compute the vertices_per_processor
num_procs=6;
for lev=1:L
    Nlev=mesh{lev}.N;
    % if Nlev=8, num_procs=3, min_num_nodes=2 and
    % processors_with_more_nodes=2 (they have 3 nodes)
    min_num_nodes=floor(Nlev/num_procs);
    processors_with_more_nodes=Nlev-min_num_nodes*num_procs;
    
    for pp=1:processors_with_more_nodes
        vertices_per_processor{lev}{pp}= (1+(pp-1)*(min_num_nodes+1) ): pp*(min_num_nodes+1) ;
    end
   % add_me_if_previous_cycle_not_empty=sign(processors_with_more_nodes);
    for pp=processors_with_more_nodes+1:num_procs
        vertices_per_processor{lev}{pp}= (1+processors_with_more_nodes+(pp-1)*(min_num_nodes) ): (processors_with_more_nodes+pp*(min_num_nodes)) ;
    end
end






for mm=1:parameters.max_iter

%ArnoldVcycleContact(AFinenobcnt,Complementarity,bnt,x,h,c,WorkingSetContact{mm},parameters.smoothing_steps,mesh_parameters,graph,Pnt,L,L);
[x,residual(mm)] = ArnoldVcycleUTOPIA(AFinenobcnt,bnt,x,vertices_per_processor,parameters.smoothing_steps,mesh_parameters,graph,Pnt,L,L)
energy(mm)=0.5*x'*AFinenobcnt*x-bnt'*x;
[mm,log10(residual(mm)),energy(mm)]

if(residual(mm)<parameters.toll_loop)
break
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRY UTOPIA SMOOTHER PARALLEL VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 END                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



isitlinear=[];
for mm=1:parameters.max_iter

% allvars = whos;
% memused = sum([allvars.bytes])
 [x,WorkingSetContact{mm+1},residual(mm)] = ArnoldVcycleContact(AFinenobcnt,Complementarity,bnt,x,h,c,WorkingSetContact{mm},parameters.smoothing_steps,mesh_parameters,graph,Pnt,L,L);
%[x,WorkingSetContact{mm+1},residual(mm)] = ArnoldVcycleContactNoCoarseConstraint(AFinenobcnt,Complementarity,bnt,x,h,c,WorkingSetContact{mm},parameters.smoothing_steps,mesh_parameters,graph,Pnt,L,L);

energy(mm)=0.5*x'*AFinenobcnt*x-bnt'*x;
[mm,log10(residual(mm)),energy(mm)]
if(~isequal(WorkingSetContact{mm},WorkingSetContact{mm+1}))
   isitlinear(end+1)=0
else
    isitlinear(end+1)=1
end
activeset=find(WorkingSetContact{mm}==1);
save('activesetlambda1mu1MEAN.mat','activeset');
save('isitlinearlambda1mu1MEAN.mat','isitlinear')
if(residual(mm)<parameters.toll_loop)
break
end

end


% string1='SquareResidualLambda1e50Mu1C1F';
% string2=num2str(parameters.FINE);
% string3='.mat';
% string=strcat(string1,string2,string3);
% residual=log10(residual);
% save(string,'residual');
householder=1;
[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L},householder);
sol=M_Normal_Tangent*x;
% print_displacement_solution(mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');
p = ContactPressure(mesh,sol,parameters);
ContactDisplacement(mesh,sol,parameters);
string0='';
string1='SignoriniNoConstraint';
string1=strcat(string1,string0,'Residuallambda');
string2=num2str(parameters.lambda);
string3='mu';
string4=num2str(parameters.mu);
string5='Ceq';
string6=num2str(parameters.C_eq);
string7='Casym';
string8=num2str(parameters.C_asym);
string9='Csmoothing';
string10=num2str(parameters.smoothing_steps);
string11='.mat';
stringext=strcat(string2,string3,string4,string5,string6,string7,string8,string9,string10,string11);
string=strcat(string1,stringext);
save(string,'residual');


string1='SignoriniNoConstraint';
string1=strcat(string1,string0,'Solutionlambda');
string=strcat(string1,stringext);
save(string,'sol');




% figure
% plot(log10(residual))
% hold on
% plot(log10(resz))
% top2bottom=1;
% lev=L;
% is_on_coarser_grid=0;
% smoothing_steps=10;
% figure