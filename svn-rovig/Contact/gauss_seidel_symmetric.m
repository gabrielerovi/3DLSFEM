function[x]=gauss_seidel_symmetric(A,b,x,smoothing_steps,is_symmetric)

    
n=length(b);
if(n>1)
    
if(is_symmetric==false)
for j=1:smoothing_steps
x_old=x;
x(1)=1/A(1,1) *(b(1) -x_old(2:n)'*A(1,2:n)' );
for i=2:n
  
    if(A(i,i)==0)
        erroreeeeee=1
    end
    x(i)=1/A(i,i) * (b(i) - x(1:i-1)'*A(i,1:i-1)' -x_old(i+1:n)'*A(i,i+1:n)' );


end

end

else
    
    smoothing_steps=ceil(smoothing_steps/2);
 for j=1:smoothing_steps
x_old=x;
x(1)=1/A(1,1) *(b(1) -x_old(2:n)'*A(1,2:n)' );
for i=2:n
  
    if(A(i,i)==0)
        erroreeeeee=1
    end
    x(i)=1/A(i,i) * (b(i) - x(1:i-1)'*A(i,1:i-1)' -x_old(i+1:n)'*A(i,i+1:n)' );


end
x_old=x;
x(n)=1/A(n,n) *(b(n) -x_old(1:n-1)'*A(n,1:n-1)' );
for i=n-1:1
  
    if(A(i,i)==0)
        erroreeeeee=1
    end
    x(i)=1/A(i,i) * (b(i) - x_old(1:i-1)'*A(i,1:i-1)' -x(i+1:n)'*A(i,i+1:n)' );


end

end   
end



else
    x(1)=b(1)/A(1,1);
end

end