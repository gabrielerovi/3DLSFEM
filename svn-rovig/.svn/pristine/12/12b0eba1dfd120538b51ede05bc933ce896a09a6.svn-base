function Schur= SchurComplement(Asu,Ass,Auu,bool,mesh)

% we consider N_remove for the displacement in x an y,
% so we add N_remove+N to N_remove 
N_remove=cat(2,mesh.N_remove,mesh.N_remove+mesh.N)';

% approximate Ass=diag(Ass)
if(bool==false)
DAss=diag(Ass).^(-1); 
Assunew=DAss.*Asu;
Schur=Asu'* Assunew;
else
% use Ass
Assminus1=Ass^(-1);
Schur=Asu'* Assminus1*Asu;
end

Schur=Auu- Schur;
Schur(N_remove,:)=[];
Schur(:,N_remove)=[];
end