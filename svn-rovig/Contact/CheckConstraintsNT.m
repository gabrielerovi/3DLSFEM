function [Constraint] = CheckConstraintsNT(mesh,parameters)


dim=parameters.dim;
L=size(mesh);
L=L(1);
grid=mesh{L};
N=grid.N;
NF=grid.NF;
node=grid.node;
F_contact=grid.F_contact;
N_contact=grid.N_contact;
contF=0;
RhsN=[];
RhsF=[];

for ff=F_contact        
contF=contF+1;
RhsF(contF,1)=0; 
end


contN=0;

for nn=N_contact        
contN=contN+1;
vertices=node(nn,[1:dim]);
RhsN(contN,1)=gap_function(vertices(1),vertices(2),vertices(3) );
end

Constraint.RhsF=RhsF;
Constraint.RhsN=RhsN;

  
end