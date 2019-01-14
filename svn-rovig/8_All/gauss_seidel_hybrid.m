function[x]=gauss_seidel_hybrid(components,A,b,x,NDtoRT,C,smoothing_steps,mesh)
x=gauss_seidel(A,b,x,smoothing_steps);

residualND=b-A*x;

residualND=NDtoRT'*residualND;

err=gauss_seidel(C,residualND,zeros(length(C(:,1)),1),smoothing_steps);
%err=gauss_seidel_const_bc(components,C,residualND,zeros(length(C(:,1)),1),smoothing_steps,mesh);

correction=NDtoRT*err;

x=x+correction;

end
