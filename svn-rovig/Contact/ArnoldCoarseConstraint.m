function ConstraintC= ArnoldCoarseConstraint(ConstraintF,mesh_parameters)

NEC=mesh_parameters.NEC;
NEF=mesh_parameters.NEF;
NC=mesh_parameters.NC;
NF=mesh_parameters.NF;
EC_contact=mesh_parameters.EC_contact;
% now we want to build coarse gap function
% in particular the constraint for the pressure is always negative
% for the constraint on the displacement in normal direction is take care
% with proper Monotone Restriction
ConstraintC=sparse(2*NEC+2*NC,1)+10^10;
ConstraintC(EC_contact)=0;
ConstraintNF=ConstraintF(2*NEF+1:2*NEF+NF);
for nC=1:NC
    tmp=10^10;
    for nb=1:length(mesh_parameters.Patch_Node_Monotone{nC})
    points=mesh_parameters.Patch_Node_Monotone{nC}{nb};
    c_ext=ConstraintNF(points(1));
    c_midpoint=ConstraintNF(points(2));
    
    %tmp=max(tmp, 2 * c_midpoint - c_ext);
    
    %tmp= min(tmp,  max(c_midpoint,2 * c_midpoint - c_ext ) );
    tmp=min(tmp,min(c_midpoint,c_ext));
    %tmp=min(tmp,c_midpoint);
    end
    
    ConstraintC(2*NEC+nC)=min(ConstraintNF(nC), tmp);
end

ConstraintEF=ConstraintF(1:2*NEF);
for eC=1:NEC
    points=mesh_parameters.Patch_Edge_Monotone{eC};    
    ConstraintC(eC)=min(ConstraintEF(points) );
end


end