function[x]=gauss_seidel(A,b,x,smoothing_steps)

n=length(b);
if(n>1)
for j=1:smoothing_steps
x_old=x;
x(1)=1/A(1,1) *(b(1) -x_old(2:n)'*A(1,2:n)' );
for i=2:n
  
    x(i)=1/A(i,i) * (b(i) - x(1:i-1)'*A(i,1:i-1)' -x_old(i+1:n)'*A(i,i+1:n)' );


    energy(i)=0.5*x'*A*x-b'*x;

end
    hold on
    plot(energy)
end
else
    x(1)=b(1)/A(1,1);
end

end