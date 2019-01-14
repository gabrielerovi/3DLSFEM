function [Constraint] = CheckConstraintsNT(mesh,parameters)


gap=parameters.gap;
L=size(mesh);
L=L(1);
grid=mesh{L};
N=grid.N;
NE=grid.NE;
node=grid.node;
E_contact=grid.E_contact;
N_contact=grid.N_contact;
contE=0;
RhsN2=[];
RhsE2=[];

for ee=E_contact        
contE=contE+1;
RhsE2(contE,1)=0; 
end


contN=0;

for nn=N_contact        
contN=contN+1;
vertices=node(nn,[1,2]);
RhsN2(contN,1)=gap_function(vertices(1),vertices(2) );
end

Constraint.RhsE2=RhsE2;
Constraint.RhsN2=RhsN2;

  
end