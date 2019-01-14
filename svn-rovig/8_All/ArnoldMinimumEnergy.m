function sol = ArnoldMinimumEnergy (A,b,Constraint,ContactDofs)

% ContactDofs is a vector containing all the dofs that can be in contact
% we then build all the possible combinations

toll=10^(-16);


cont=0;
for ii=1:length(ContactDofs)
tmp=combnk(ContactDofs,ii);

for jj=1:length(tmp(:,1))
    cont=cont+1;
    combo{cont}=tmp(jj,:);
    
end

end

combo{end+1}=[];



for ii=1:length(combo)
    
    Atmp=A;
    btmp=b;
    Atmp(combo{ii},:)=0;
    for jj=combo{ii}
       Atmp(jj,jj)=1;
       btmp(jj)=Constraint(jj);
    end
    
%     Atmp=Atmp'*Atmp;
%     btmp=Atmp'*btmp;
    c=Atmp\btmp;
    
    % we are considering the corrections, so Constraint=C-x_k
    
    if(sum(find(c-Constraint>toll))>0)
        J(ii)=10^10;
    else
        v=c;
        v(combo{ii})=0;
        J(ii)=0.5 * v'*A * v - b'*v;
        C_combo(ii,:)=c;
    end
        
end

[minJ,minJpos]=min(J);
sol=C_combo( minJpos,:);


if(isrow(sol))
    sol=sol';
end

end