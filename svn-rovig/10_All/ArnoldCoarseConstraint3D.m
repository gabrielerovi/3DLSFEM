function ConstraintC= ArnoldCoarseConstraint3D(ConstraintF,mesh_parameters)

NFC=mesh_parameters.NFC;
NFF=mesh_parameters.NFF;
NC=mesh_parameters.NC;
NF=mesh_parameters.NF;
FC_contact=mesh_parameters.FC_contact;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Now we want to build coarse gap function
%%%%% in particular the constraint for the pressure is always  pos
%%%%% for the constraint on the displacement in normal direction is take care
%%%%% with proper Monotone Restriction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConstraintC=sparse( 3 * (NFC+NC),1)+10^10;
ConstraintC(FC_contact)=0;
ConstraintNF=ConstraintF(3*NFF+1:3*NFF+NF);
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
    
    ConstraintC(3*NFC+nC)=min(ConstraintNF(nC), tmp);
end

ConstraintEF=ConstraintF(1:3*NFF);
for fC=1:NFC
    points=mesh_parameters.Patch_Face_Monotone{fC};    
    ConstraintC(fC)=min(ConstraintEF(points) );
end


end