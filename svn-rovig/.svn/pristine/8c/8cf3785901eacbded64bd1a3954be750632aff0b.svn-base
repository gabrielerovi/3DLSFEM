function x =v_cycle_hybrid(components,lev,A_lev,b,x,C_lev,RTCtoRTF,NDtoRT,smoothing_steps,mesh)
if(lev==1)
    x=A_lev{1}\b;
else
    %pre-smoothing
    x=gauss_seidel_hybrid(components,A_lev{lev},b,x,NDtoRT{lev},C_lev{lev},smoothing_steps,mesh{lev})   ;
    residualF=b-A_lev{lev}*x;
    residualC=RTCtoRTF{lev-1}'*residualF;
    correctionC=zeros(length(residualC),1);
    % v-cycle
    correctionC =v_cycle_hybrid(components,lev-1,A_lev,residualC,correctionC,C_lev,RTCtoRTF,NDtoRT,smoothing_steps,mesh);
    correctionF=RTCtoRTF{lev-1}*correctionC;
    
    Ac=A_lev{lev} *correctionF;
    alpha= (residualF'*correctionF)/(Ac'*correctionF);
    x=x+alpha*correctionF;
    
%    x=x+correctionF;
    
    
    %post-smoothing
    x=gauss_seidel_hybrid(components,A_lev{lev},b,x,NDtoRT{lev},C_lev{lev},smoothing_steps,mesh{lev})   ;
    
end

end
