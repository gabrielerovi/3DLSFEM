%% 
syms uxx uxy uyx uyy
syms vxx vxy vyx vyy
syms sxx sxy syx syy
syms dx_sxx dy_syy
syms txx txy tyx tyy
syms dx_txx dy_tyy
syms phi1x phi1y phi2x phi2y beta alpha
syms psi1x psi1y psi2x psi2y
syms C1 C2


sigma=[ sxx, sxy;syx,syy];
tau=[txx,txy;tyx tyy];
divsigma= [dx_sxx; dy_syy];
divtau= [dx_txx;dy_tyy];
Asigma=beta *sigma + alpha * (sxx+syy) *[1 0; 0 1];
Atau=  beta * tau + alpha * (txx+tyy) *[1 0; 0 1];

 epsilonu= [uxx   0.5*(uxy+uyx);       0.5*(uxy+uyx) uyy];
 epsilonv= [vxx   0.5*(vxy+vyx);       0.5*(vxy+vyx) vyy];

Constitutive1=Asigma -epsilonu;
Constitutive2=Atau - epsilonv;

F=  C2 * sum(sum(Constitutive1 .*Constitutive2))
F= F + C1 * sum(divsigma' * divtau);


F1=diff(F,txx)* txx +diff(F,txy)* txy +diff(F,dx_txx)* dx_txx;
F2=diff(F,tyx)* tyx +diff(F,tyy)* tyy +diff(F,dy_tyy)* dy_tyy;
F3=diff(F,vxx)*vxx+diff(F,vxy)*vxy;
F4=diff(F,vyx)*vyx+diff(F,vyy)*vyy;

F11=diff(F1,sxx)*sxx+diff(F1,sxy)*sxy+diff(F1,dx_sxx)*dx_sxx;
F12=diff(F1,syx)*syx+diff(F1,syy)*syy+diff(F1,dy_syy)*dy_syy;
F13=diff(F1,uxx)*uxx+diff(F1,uxy)*uxy;
F14=diff(F1,uyx)*uyx+diff(F1,uyy)*uyy;


F22=diff(F2,syx)*syx+diff(F2,syy)*syy+diff(F2,dy_syy)*dy_syy;
F23=diff(F2,uxx)*uxx+diff(F2,uxy)*uxy;
F24=diff(F2,uyx)*uyx+diff(F2,uyy)*uyy;

F33=diff(F3,uxx)*uxx+diff(F3,uxy)*uxy;
F34=diff(F3,uyx)*uyx+diff(F3,uyy)*uyy;
F44=diff(F4,uyx)*uyx+diff(F4,uyy)*uyy;


%% epsilonu=[uxx           0.5*(uxy+uyx);
syms uxx uxy uyx uyy
syms phi1x phi1y phi2x phi2y beta alpha
syms psi1x psi1y psi2x psi2y

 epsilonu= [uxx   0.5*(uxy+uyx);       0.5*(uxy+uyx) uyy];
      
Aphi=[phi1y, - phi1x;
            phi2y, - phi2x;];
Apsi=[psi1y, - psi1x;
      psi2y, - psi2x;];

 Aphi=beta*Aphi + alpha * (phi1y-phi2x)*[1 0; 0 1];
 Apsi=beta*Apsi + alpha * (psi1y-psi2x)*[1 0; 0 1];
 
D=sum(sum(Aphi.*Apsi))
C=-sum(sum(Aphi.*epsilonu))

C1=diff(C,phi1x)*phi1x+diff(C,phi1y)*phi1y;
C2=diff(C,phi2x)*phi2x+diff(C,phi2y)*phi2y;

C13=diff(C1,uxx)*uxx+diff(C1,uxy)*uxy;
C14=diff(C1,uyx)*uyx+diff(C1,uyy)*uyy;

C23=diff(C2,uxx)*uxx+diff(C2,uxy)*uxy;
C24=diff(C2,uyx)*uyx+diff(C2,uyy)*uyy;

D1=diff(D,psi1x)*psi1x+diff(D,psi1y)*psi1y;
D2=diff(D,psi2x)*psi2x+diff(D,psi2y)*psi2y;

D11=diff(D1,phi1x)*phi1x+diff(D1,phi1y)*phi1y;
D12=diff(D1,phi2x)*phi2x+diff(D1,phi2y)*phi2y;
D22=diff(D2,phi2x)*phi2x+diff(D2,phi2y)*phi2y;