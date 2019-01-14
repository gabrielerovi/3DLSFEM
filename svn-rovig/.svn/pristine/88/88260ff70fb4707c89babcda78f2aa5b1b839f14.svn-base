clear all
close all
clc
%echo on

% data

parameters =parameters_lambdainfsmoothing5();


[mesh,h]=mesh(parameters);
L=length(mesh);



dim=parameters.dim;
Dofs=3*(mesh{1}.N+mesh{1}.NF)
Dofs=3*(mesh{L}.N+mesh{L}.NF)
fprintf('lambda= %f\n',parameters.lambda)
fprintf('smoothing= %f\n',parameters.smoothing_steps)
% for ii=1:length(mesh)
% vv{ii}=0;
% end
% print_mesh_L(mesh,vv);

[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L});


[maps,graph]=graph_and_maps(mesh);
maps=ArnoldMonotoneConstraint(mesh,maps,parameters);
Constraint = CheckConstraintsNT(mesh,parameters);






[AFinenobcnt,Complementarity,bnt,bnt1,Pnt]=create_system_contact(parameters,mesh,h,maps.Patch_Node_Monotone);



% Aprova=AFinenobcnt;
% Aprova(mesh{L}.Remove,:)=0;
% for ii=mesh{L}.Remove
% Aprova(ii,ii)=1;
% end
% % remove the tangent components of sigma n
% Aprova(mesh{L}.F_contact+mesh{L}.NF,:)=0;
% Aprova(mesh{L}.F_contact+mesh{L}.NF*2,:)=0;
% for ii=mesh{L}.F_contact+mesh{L}.NF
% Aprova(ii,ii)=1;
% end
% for ii=mesh{L}.F_contact+2*mesh{L}.NF
% Aprova(ii,ii)=1;
% end
% % remove normal components of displacement 
% Aprova(mesh{L}.N_contact+mesh{L}.NF*3,:)=0;
% for ii=mesh{L}.N_contact+mesh{L}.NF*3
% Aprova(ii,ii)=1;
% end
% x=Aprova\bnt;
% ContactPressure(mesh,x,parameters);
% xnt=M_Normal_Tangent'*x;
% print_displacement_solution3D(mesh,xnt);
% 
% 
% 
% 
% 
% Aprova=M_Normal_Tangent'*AFinenobcnt*M_Normal_Tangent;
% 
% Aprova(mesh{L}.Remove,:)=0;
% for ii=mesh{L}.Remove
% Aprova(ii,ii)=1;
% end
% 
% % remove x,y components of sigma n
% Aprova(mesh{L}.F_contact,:)=0;
% Aprova(mesh{L}.F_contact+mesh{L}.NF,:)=0;
% for ii=mesh{L}.F_contact
% Aprova(ii,ii)=1;
% end
% for ii=mesh{L}.F_contact+mesh{L}.NF
% Aprova(ii,ii)=1;
% end
% 
% % remove z components of displacement 
% Aprova(mesh{L}.N_contact+mesh{L}.NF*3+mesh{L}.N*2,:)=0;
% for ii=mesh{L}.N_contact+mesh{L}.NF*3+mesh{L}.N*2
% Aprova(ii,ii)=1;
% end
% x=Aprova\bnt;
% print_displacement_solution3D(mesh,x);
% xnt=M_Normal_Tangent'*x
% ContactPressure(mesh,xnt,parameters);
% 
% 
% Dofs=3*(mesh{L}.N+mesh{L}.NF)

%  AFinenobcnt=add_boundary_bc_system(AFinenobcnt,mesh{L}.RemoveNT);
%  
%  sol=AFinenobcnt\bnt;
%  [M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L});
%  sol=M_Normal_Tangent*sol;
%  print_displacement_solution3D(mesh,sol)
big_value=10^10;
c=big_value*ones(dim * (mesh{L}.NF+mesh{L}.N),1);
c(mesh{L}.F_contact)=Constraint.RhsF;
c(mesh{L}.N_contact+dim*mesh{L}.NF)=Constraint.RhsN;
WorkingSet_Loc=[];
B=speye(length(bnt));


% Aprova=AFinenobcnt+Complementarity;
% Aprova(mesh{L}.RemoveNT,:)=0;
% for ii=mesh{L}.RemoveNT
% Aprova(ii,ii)=1;
% end
%  C11=find(c<10);
%  C=B(C11,:);
%  gamma=1;
%  C=gamma*C'*C;
%  eigA=eig(full(Aprova))
% eigC=eig(full(Aprova+C))
% [x,lambda,WorkingSet] = ArnoldActiveset(Aprova,B,bnt,c,WorkingSet_Loc);
% ContactPressure(mesh,x,parameters);
% xnt=M_Normal_Tangent'*x;
% print_displacement_solution3D(mesh,xnt);
% IT IS ALREADY HERE DO NOT ADD IT LATER



% WorkingSetContact{1}=sparse(mesh{L}.NE*2+mesh{L}.N*2,1);
for lev=1:L
    mesh_parameters{lev}.RemoveNT=mesh{lev}.RemoveNT;
    mesh_parameters{lev}.NF=mesh{lev}.NF;
    mesh_parameters{lev}.N=mesh{lev}.N;
    mesh_parameters{lev}.F_contact=mesh{lev}.F_contact;
    mesh_parameters{lev}.N_contact=mesh{lev}.N_contact;
    mesh_parameters{lev}.Remove=mesh{lev}.Remove;
    mesh_parameters{lev}.RemoveNT=mesh{lev}.RemoveNT;
    mesh_parameters{lev}.Patch_Internal_All=maps.Patch_Internal_All{lev};
    if(lev<L)
    mesh_parameters{lev}.Patch_Face_Monotone=maps.Patch_Face_Monotone{lev};
    mesh_parameters{lev}.Patch_Node_Monotone=maps.Patch_Node_Monotone{lev};
    end
end


% clearvars mesh
% clearvars maps


%%% guarda meglio ArnoldCoarseConstraint3D
[x,WorkingSetContact{1}] = ArnoldNestedIteration3D(AFinenobcnt,Complementarity,h,bnt1,c,mesh_parameters,Pnt);
Dofs=3*(mesh{L}.N+mesh{L}.NF)

load('x/xlambdainfsmoothing5.mat');
isitlinear=[];
ffname{1}='x';
ffname{2}='sol';
ffname{3}='residuals';
ffname{4}='activeset';
ffname{5}='isitlinear';

for mm=1:parameters.max_iter

% allvars = whos;
% memused = sum([allvars.bytes])
[x,WorkingSetContact{mm+1},residual(mm)] = ArnoldVcycleContact3D(AFinenobcnt,Complementarity,bnt,x,h,c,WorkingSetContact{mm},parameters.smoothing_steps,mesh_parameters,graph,Pnt,L,L);

mm
log10(residual)
if(~isequal(WorkingSetContact{mm},WorkingSetContact{mm+1}))
   isitlinear(end+1)=0;
else
    isitlinear(end+1)=1;
end
activeset=find(WorkingSetContact{mm}==1);

save( fullfile( pwd,ffname{4},'activesetlambdainfsmoothing5.mat'), 'activeset')
save( fullfile( pwd,ffname{5},'isitlinearlambdainfsmoothing5.mat'), 'isitlinear')
save( fullfile( pwd,ffname{1},'xlambdainfsmoothing5.mat'), 'x')
if(residual(mm)<parameters.toll_loop)
break
end

end





% ContactPressure(mesh,x,parameters);
[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{parameters.L});
sol=M_Normal_Tangent*x;
% print_displacement_solution3D(mesh,sol)
% ContactDisplacement(mesh,sol,parameters);
% savevtkdisplacement(mesh{L},sol)






stringtype='CubeContact';
string0='';
string1=stringtype;
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
string11='Dofs';
string12=num2str(3*(mesh{L}.N+mesh{L}.NF));
string13='C1F';
string14=num2str(L);
string15='.mat';
stringext=strcat(string2,string3,string4,string5,string6,string7,string8,string9,string10,string11,string12,string13,string14,string15);
string=strcat(string1,stringext);
save( fullfile( pwd,ffname{3},'residuallambdainfsmoothing5.mat'), 'residual')


string1=stringtype;
string1=strcat(string1,string0,'sol/Solutionlambda');
string=strcat(string1,stringext);
save( fullfile( pwd,ffname{2},'sollambdainfsmoothing5.mat'), 'sol')


string1=stringtype;
string1=strcat(string1,string0,'x/SolutionNTlambda');
string=strcat(string1,stringext);
save( fullfile( pwd,ffname{1},'xlambdainfsmoothing5.mat'), 'x')








% string1='CubeResidualLambda1e50Mu1C1F';
% string2=num2str(parameters.FINE);
% string3='.mat';
% string=strcat(string1,string2,string3);
% residual=log10(residual);
% save(string,'residual');
% [erruL2,errgrad,errsigmaL2,errdivsigmaL2]=error_displacement(mesh,sol,parameters)


% print_displacement_solution3D(mesh,sol)
% figure
% plot(log10(residual))
% hold on
% plot(log10(resz))
% top2bottom=1;
% lev=L;
% is_on_coarser_grid=0;
% smoothing_steps=10;
% figure
