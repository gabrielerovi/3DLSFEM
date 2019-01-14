
function [A,b,AFine,P,Ant,bnt,Antbc,AFinenobc,AFinenobcnt,bnt_lev,C]=create_system(parameters,mesh,h)
qrule=parameters.qrule;
L=length(mesh);
toll=10^-12;
% create projections
P.NDtoRT=P1_to_RT(mesh);
P.RTCtoRTF= RTCtoRTF2D(mesh,qrule);
P.P1CtoP1F =P1CtoP1F2D(mesh,qrule);
P.NDCtoNDF=P.P1CtoP1F;
P.P1CtoP1Fuu=[];
P.RTCtoRTFss=[];
Antbc=[];
if(L>1)
for lev=1:L-1
% here we remove errors in the computation of P1 projection    
P.P1CtoP1F{lev}(abs(P.P1CtoP1F{lev})<toll)=0;

[m,n]=size(P.RTCtoRTF{lev});
P.RTCtoRTFss{lev}= [P.RTCtoRTF{lev} zeros(m,n);
                    zeros(m,n)       P.RTCtoRTF{lev}];
[m,n]=size(P.P1CtoP1F{lev});
P.P1CtoP1Fuu{lev}= [P.P1CtoP1F{lev} zeros(m,n);
                    zeros(m,n)      P.P1CtoP1F{lev}];
end
end
N_components=parameters.N_components;
E_components=parameters.E_components;

dim_problem=N_components+E_components;
A_eq=cell(L,dim_problem,dim_problem);
A_const_asym=cell(L,dim_problem,dim_problem);
A=cell(L,dim_problem,dim_problem);


C_eq=parameters.C_eq;
C_const=parameters.C_const;
C_asym=parameters.C_asym;
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;
contact=parameters.contact;

if(parameters.h_dependence==false)
    h=ones(L,1);
end


switch N_components
    case 0
        N_remove=[];
        N_label=[];
    otherwise
        N_remove=[mesh{L}.N_remove];
        N_label=mesh{L}.N_label;
end

switch E_components
    case 0
        E_remove=[];
        E_label=[];
    otherwise
        E_remove=[mesh{L}.E_remove];
        E_label=mesh{L}.E_label;
end




 Anobc=[];   
 bnt_lev=[];
 C=[];
if(strcmp(parameters.input_name,'LSstressblock'))
[A,b,AFine,P]=create_system_LSstressblock(parameters,mesh,h,P);
elseif(strcmp(parameters.input_name,'LSelasticity')||strcmp(parameters.input_name,'LSelasticityAsymmetric'))
[A,b,AFine,P, Ant,bnt,Antbc,AFinenobc,AFinenobcnt,bnt_lev,C]=create_system_LSelasticity(parameters,mesh,h,P);    
elseif(strcmp(parameters.input_name,'DispElasticity'))
[A,b,AFine,P ]=create_system_DispElasticity(parameters,mesh,h,P);  
end



end



