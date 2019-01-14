function x =v_cycle(lev,A_lev,b,x,P1CtoP1F,smoothing_steps,is_symmetric,mesh)



if(lev==1)
    x=A_lev{1}\b;
else
    
%     %pre-smoothing
%     perm=symrcm(A_lev{lev});
%     perm_inv(perm)=1:length(perm);
%     x=gauss_seidel_symmetric(A_lev{lev}(perm,perm),b(perm),x(perm),smoothing_steps,is_symmetric) ;
%     x=x(perm_inv);
%     residualF=b-A_lev{lev}*x;
%     residualC=P1CtoP1F{lev-1}'*residualF;
%     correctionC=zeros(length(residualC),1);
%     % v-cycle
%     correctionC =v_cycle(lev-1,A_lev,residualC,correctionC,P1CtoP1F,smoothing_steps);
%     correctionF=P1CtoP1F{lev-1}*correctionC;
%     x=x+correctionF;
%     %post-smoothing
%     perm=symrcm(A_lev{lev});
%     perm_inv(perm)=1:length(perm);
%     x=gauss_seidel_symmetric(A_lev{lev}(perm,perm),b(perm),x(perm),smoothing_steps,is_symmetric) ;
%     x=x(perm_inv);    
    
    
    
    
C=lev-1;
F=lev;
NC=mesh{C}.N;
NF=mesh{F}.N;
A=A_lev{F};

    %pre-smoothing
    
    
    [x]=gauss_seidel(A_lev{lev},b,x,smoothing_steps)
    
%     smoothing_steps=100;
%     
%     perm=symrcm(A_lev{lev});
%     perm_inv(perm)=1:length(perm);
%     x=gauss_seidel_symmetric(A_lev{lev}(perm,perm),b(perm),x(perm),smoothing_steps,is_symmetric) ;
%     x=x(perm_inv);
    residualF=b-A*x;
    residualC=P1CtoP1F{C}'*residualF;
    
    N_removeC=mesh{C}.N_remove;
    N_removeC=[N_removeC,N_removeC+NC];
    
    residualC(N_removeC)=0;
    
    
    
    
    correctionC=zeros(length(residualC),1);
    % v-cycle
    correctionC =v_cycle(lev-1,A_lev,residualC,correctionC,P1CtoP1F,smoothing_steps,is_symmetric,mesh);
    correctionF=P1CtoP1F{lev-1}*correctionC;
    
    N_removeF=mesh{F}.N_remove;
    N_removeF=[N_removeF,N_removeF+NF];
    correctionF(N_removeF)=0;
    
    
    N=[length(A(:,1))/2; length(A(:,1))];
    %N=length(A(:,1))
    d=1.0*ones(2,1);
    %[d]=doubleLineSearch (A(1:N(1),1:N(1)),A(1:N(1),N(1)+1:N(2)),A(N(1)+1:N(2),1:N(1)),A(N(1)+1:N(2),N(1)+1:N(2)),correctionF(1:N(1)), correctionF(N(1)+1:N(2)),residualF(1:N(1)), residualF(N(1)+1:N(2)))
    x(1:N(1))=x(1:N(1))+ d(1) *correctionF(1:N(1));
    x(1+N(1):end)=x(1+N(1):end)+ d(2) *correctionF(1+N(1):end);
    %x=x+ alpha * correctionF;
    %post-smoothing
    [x]=gauss_seidel(A_lev{lev},b,x,smoothing_steps)
%     perm=symrcm(A_lev{lev});
%     perm_inv(perm)=1:length(perm);
%     x=gauss_seidel_symmetric(A_lev{lev}(perm,perm),b(perm),x(perm),smoothing_steps,is_symmetric) ;
%     x=x(perm_inv);
end

end
