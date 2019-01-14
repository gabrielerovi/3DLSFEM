function [A,b]=add_boundary_bc_system_symmetric(A,b,Remove)

A(Remove,:)=0;
bbc=zeros(length(b),1);
bbc(Remove)=b(Remove);
Aprova=A;
Atimesbbc=A*bbc;
Atimesbbc(Remove)=0;
A(:,Remove)=0;
for ii=Remove
    A(ii,ii)=1;
end

b=b-Atimesbbc;
    
end