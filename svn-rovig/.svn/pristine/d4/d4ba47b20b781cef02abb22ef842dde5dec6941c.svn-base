function[x]=gauss_seidel_hybridREORDER(components,A,b,x,NDtoRT,C,smoothing_steps,mesh)

reorder=symrcm(A);
inversereorder(reorder) = 1:length(reorder);

Areorder=A(reorder,reorder);
breorder=b(reorder);
xreorder=x(reorder);
xreorder=gauss_seidel(Areorder,breorder,xreorder,smoothing_steps);

x=xreorder(inversereorder);
residualND=b-A*x;

residualND=NDtoRT'*residualND;

reorderC=symrcm(C);
inversereorderC(reorderC) = 1:length(reorderC);
Creorder=C(reorderC,reorderC);
residualNDreorderC=residualND(reorderC);


errreorderC=gauss_seidel(Creorder,residualNDreorderC,zeros(length(C(:,1)),1),smoothing_steps);
%err=gauss_seidel_const_bc(components,C,residualND,zeros(length(C(:,1)),1),smoothing_steps,mesh);

err=errreorderC(inversereorderC);
correction=NDtoRT*err;

x=x+correction;

end
