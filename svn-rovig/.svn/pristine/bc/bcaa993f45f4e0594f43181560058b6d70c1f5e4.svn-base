function [A_condensation,b_condensation]=add_bc_condensation(A,b,Remove)



b_condensation=sparse(length(b),1);
b_condensation(Remove)=b(Remove);
b_condensation=A*b;
b_condensation(Remove)=0;
b_condensation=b-b_condensation;
A_condensation=A;
A_condensation(Remove,:)=0;
A_condensation(:,Remove)=0;
for ii=Remove
A_condensation(ii,ii)=1;
end



end