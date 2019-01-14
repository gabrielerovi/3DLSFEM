function c=ContactProjection(x,c,Constraint,mesh)

% x is the current solution
% c is the correction
 y = x + c;
 
NE=mesh.NE;


cont=0;
for ee=mesh.E_contact
    cont=cont+1;
    tmp=Constraint.CheckConstraintE2(cont,:)*y;
    if( tmp > Constraint.RhsE2(cont))
         c(ee)=Constraint.RhsE2(cont)-tmp;
    end
end

cont=0;
for nn=mesh.N_contact
    cont=cont+1;
    tmp=Constraint.CheckConstraintN2(cont,:)*y;
    if( tmp > Constraint.RhsN2(cont))
         c(nn+2*NE)=Constraint.RhsN2(cont)-tmp;
    end

end





end