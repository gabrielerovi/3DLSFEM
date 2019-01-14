function x=gauss_seidel_const_bc(components,A,b,x,smoothing_steps,mesh)

Lengths=mesh.ND_lengths;
n=length(b);
no_boundary=0;
% impose potential = 0 fon a fixed boundary for each potential component
% (two in 2D)
for i=1:components*Lengths(1)
    x(i)=0;
end
x_old=x;
if(n>1)
for j=1:smoothing_steps
x_old=x;


 S=components*Lengths(1);
 
 for j=2:length(Lengths)-1
 i=S+1;
 x(i)=1/A(i,i) * (b(i) - x(1:i-1)'*A(i,1:i-1)' -x_old(i+1:n)'*A(i,i+1:n)' );
 x(i+1)=1/A(i+1,i+1) * (b(i+1) - x(1:(i+1)-1)'*A(i+1,1:(i+1)-1)' -x_old((i+1)+1:n)'*A(i+1,(i+1)+1:n)' );
 for i=S+3:2:S+components*Lengths(j)
     % impose potential phi1=const, potential phi2=const
     x(i)=x(S+1);
     x(i+1)=x(S+2);
 end
 S=S+components*Lengths(j);
 
 end
 
 if(S==0)
     x(1)=1/A(1,1) * (b(1) -x_old(1:n)'*A(1,1:n)' );
     S=S+1;
     no_boundary=-1;
 end
 for i=S+1:S+components*Lengths(end)+no_boundary
 x(i)=1/A(i,i) * (b(i) - x(1:i-1)'*A(i,1:i-1)' -x_old(i+1:n)'*A(i,i+1:n)' );
end
end
else
    x(1)=1/A(1,1) * (b(1) -x_old(1:n)'*A(1,1:n)' );    
end

end
