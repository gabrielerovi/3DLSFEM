


function [Anobc_tot,Pcoarse]=ArnoldAssembling(Anobc,P,WorkingSetNormal_E,WorkingSetNormal_N,mesh)

% from the finest level L to L-1, we have to remove the dofs relative to
% constrained nodes, that are the nodes on gammaC such that:
% n_obst'sigma n == 0
% u n_obst == g
% sigma n - n_obst'sigma n == 0 (frictionless case)
% The frictionless constrained is true for every edge dof on gammaC
% since we have the first edge GammaC dofs that are the normal components
% and the last ones that are the tangent dofs, it is sufficient to homogenous dirichlet conditions
% on these.

Pcoarse=P;
L=length(mesh);

NE=mesh{L}.NE;
N=mesh{L}.N;


% P is a projection operator, P: Coarse -> Fine

% 1) P: RTC -> RTF (on GammaC, the first component is the normal one)
% WorkingSetNormal_E=find(WorkingSetNormal_E>0);
% WorkingSetNormal_N=find(WorkingSetNormal_N>0);

% 1)edge-normal component,  2) edge-tangent component, 
% 3) node-normal component, 4) node-tangent component,
WorkingSetNormal_TOT=[WorkingSetNormal_E, WorkingSetNormal_N+2*NE];

% we put to zero all the normal edges and nodes that are constrained
Pcoarse{L-1} (WorkingSetNormal_TOT,:)=0;

%%%% . Compute A linear elastic on level L-1


 lev=L;
 Anobc_tot{lev}=[Anobc{lev,1,1} Anobc{lev,1,2} Anobc{lev,1,3} Anobc{lev,1,4};
                 Anobc{lev,2,1} Anobc{lev,2,2} Anobc{lev,2,3} Anobc{lev,2,4};
                 Anobc{lev,3,1} Anobc{lev,3,2} Anobc{lev,3,3} Anobc{lev,3,4};
                 Anobc{lev,4,1} Anobc{lev,4,2} Anobc{lev,4,3} Anobc{lev,4,4};];
   
 
for lev=L-1:-1:1
 Anobc_tot{lev}=Pcoarse{lev}'* Anobc_tot{lev+1} * Pcoarse{lev};
end


for lev=L-1:-1:1
E_remove=mesh{lev}.E_remove;
N_remove=mesh{lev}.N_remove;    
E_contact_tangent=mesh{lev}.E_contact;  
NE=mesh{lev}.NE;
N=mesh{lev}.N;
remove_bc=[E_remove, E_remove+NE, N_remove +2 *NE,N_remove + N + 2 *NE];

for rrr=remove_bc
 Anobc_tot{lev}(rrr,:)=0;
 Anobc_tot{lev}(rrr,rrr)=1;
end

remove_tangent_stress=E_contact_tangent+NE;
for rrr=remove_tangent_stress
 Anobc_tot{lev}(rrr,:)=0;
 Anobc_tot{lev}(rrr,rrr)=1;
end


end


% Anobc{lev-1,1,1} = RTCtoRTF_normal' * Anobc{lev,1,1} * RTCtoRTF_normal;
% Anobc{lev-1,1,2} = RTCtoRTF_normal' * Anobc{lev,1,2} * P.RTCtoRTF{L-1};
% Anobc{lev-1,1,3} = RTCtoRTF_normal' * Anobc{lev,1,3} * P1CtoP1F_normal;
% Anobc{lev-1,1,4} = RTCtoRTF_normal' * Anobc{lev,1,4} * P.P1CtoP1F{L-1};
% 
% Anobc{lev-1,2,1} = P.RTCtoRTF{L-1}' * Anobc{lev,2,1} * RTCtoRTF_normal;
% Anobc{lev-1,2,2} = P.RTCtoRTF{L-1}' * Anobc{lev,2,2} * P.RTCtoRTF{L-1};
% Anobc{lev-1,2,3} = P.RTCtoRTF{L-1}' * Anobc{lev,2,3} * P1CtoP1F_normal;
% Anobc{lev-1,2,4} = P.RTCtoRTF{L-1}' * Anobc{lev,2,4} * P.P1CtoP1F{L-1};
% 
% Anobc{lev-1,3,1} = P1CtoP1F_normal' * Anobc{lev,3,1} * RTCtoRTF_normal;
% Anobc{lev-1,3,2} = P1CtoP1F_normal' * Anobc{lev,3,2} * P.RTCtoRTF{L-1};
% Anobc{lev-1,3,3} = P1CtoP1F_normal' * Anobc{lev,3,3} * P1CtoP1F_normal;
% Anobc{lev-1,3,4} = P1CtoP1F_normal' * Anobc{lev,3,4} * P.P1CtoP1F{L-1};
% 
% Anobc{lev-1,4,1} = P.P1CtoP1F{L-1}' * Anobc{lev,4,1} * RTCtoRTF_normal;
% Anobc{lev-1,4,2} = P.P1CtoP1F{L-1}' * Anobc{lev,4,2} * P.RTCtoRTF{L-1};
% Anobc{lev-1,4,3} = P.P1CtoP1F{L-1}' * Anobc{lev,4,3} * P1CtoP1F_normal;
% Anobc{lev-1,4,4} = P.P1CtoP1F{L-1}' * Anobc{lev,4,4} * P.P1CtoP1F{L-1};
% 
% 
% if(L>2)
% for lev=L-1:-1:2
% 
% Anobc{lev-1,1,1} = P.RTCtoRTF{lev-1}' * Anobc{lev,1,1} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,1,2} = P.RTCtoRTF{lev-1}' * Anobc{lev,1,2} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,1,3} = P.RTCtoRTF{lev-1}' * Anobc{lev,1,3} * P.P1CtoP1F{lev-1};
% Anobc{lev-1,1,4} = P.RTCtoRTF{lev-1}' * Anobc{lev,1,4} * P.P1CtoP1F{lev-1};
% 
% Anobc{lev-1,2,1} = P.RTCtoRTF{lev-1}' * Anobc{lev,2,1} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,2,2} = P.RTCtoRTF{lev-1}' * Anobc{lev,2,2} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,2,3} = P.RTCtoRTF{lev-1}' * Anobc{lev,2,3} * P.P1CtoP1F{lev-1};
% Anobc{lev-1,2,4} = P.RTCtoRTF{lev-1}' * Anobc{lev,2,4} * P.P1CtoP1F{lev-1};
% 
% Anobc{lev-1,3,1} = P.P1CtoP1F{lev-1}' * Anobc{lev,3,1} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,3,2} = P.P1CtoP1F{lev-1}' * Anobc{lev,3,2} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,3,3} = P.P1CtoP1F{lev-1}' * Anobc{lev,3,3} * P.P1CtoP1F{lev-1};
% Anobc{lev-1,3,4} = P.P1CtoP1F{lev-1}' * Anobc{lev,3,4} * P.P1CtoP1F{lev-1};
% 
% Anobc{lev-1,4,1} = P.P1CtoP1F{lev-1}' * Anobc{lev,4,1} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,4,2} = P.P1CtoP1F{lev-1}' * Anobc{lev,4,2} * P.RTCtoRTF{lev-1};
% Anobc{lev-1,4,3} = P.P1CtoP1F{lev-1}' * Anobc{lev,4,3} * P.P1CtoP1F{lev-1};
% Anobc{lev-1,4,4} = P.P1CtoP1F{lev-1}' * Anobc{lev,4,4} * P.P1CtoP1F{lev-1};
% 
% end
% end


% for lev=L-1:-1:1
%     
% E_remove=mesh{lev}.E_remove;
% N_remove=mesh{lev}.N_remove;    
% E_contact_tangent=mesh{lev}.E_contact;  
% 
% % remove the stress tangent component (frictionless case)
% Anobc{lev,2,1} (E_contact_tangent,:)=0;
% Anobc{lev,2,2} (E_contact_tangent,:)=0;
% Anobc{lev,2,3} (E_contact_tangent,:)=0;
% Anobc{lev,2,4} (E_contact_tangent,:)=0;
% 
% for eee=E_contact_tangent
%     Anobc{lev,2,2}(eee,eee)=1;
% end
%     
% Anobc{lev,1,1} (E_remove,:)=0; 
% Anobc{lev,1,2} (E_remove,:)=0; 
% Anobc{lev,1,3} (E_remove,:)=0; 
% Anobc{lev,1,4} (E_remove,:)=0; 
% 
% Anobc{lev,2,1} (E_remove,:)=0; 
% Anobc{lev,2,2} (E_remove,:)=0; 
% Anobc{lev,2,3} (E_remove,:)=0; 
% Anobc{lev,2,4} (E_remove,:)=0; 
% 
% Anobc{lev,3,1} (N_remove,:)=0; 
% Anobc{lev,3,2} (N_remove,:)=0;  
% Anobc{lev,3,3} (N_remove,:)=0;  
% Anobc{lev,3,4} (N_remove,:)=0;  
% 
% Anobc{lev,4,1} (N_remove,:)=0; 
% Anobc{lev,4,2} (N_remove,:)=0; 
% Anobc{lev,4,3} (N_remove,:)=0;  
% Anobc{lev,4,4} (N_remove,:)=0; 
% 
% 
%     
% for eee=E_remove
%     Anobc{lev,1,1}(eee,eee)=1;
%     Anobc{lev,2,2}(eee,eee)=1;
% end
% 
% for nnn=N_remove
%     Anobc{lev,3,3}(nnn,nnn)=1;
%     Anobc{lev,4,4}(nnn,nnn)=1;
% end


% end



% for lev=L-1:-1:1
% 
% AnobcFine{lev}=[Anobc{lev,1,1} Anobc{lev,1,2} Anobc{lev,1,3} Anobc{lev,1,4};
%        Anobc{lev,2,1} Anobc{lev,2,2} Anobc{lev,2,3} Anobc{lev,2,4};
%        Anobc{lev,3,1} Anobc{lev,3,2} Anobc{lev,3,3} Anobc{lev,3,4};
%        Anobc{lev,4,1} Anobc{lev,4,2} Anobc{lev,4,3} Anobc{lev,4,4};];
%    [min(eig(abs(full(AnobcFine{lev})))),max(eig(abs(full(AnobcFine{lev}))))]
%    
% end  

end