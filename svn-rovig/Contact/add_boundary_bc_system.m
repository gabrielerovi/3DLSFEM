function A=add_boundary_bc_system(A,Remove)

A(Remove,:)=0;
for ii=Remove
    A(ii,ii)=1;
end
    

end