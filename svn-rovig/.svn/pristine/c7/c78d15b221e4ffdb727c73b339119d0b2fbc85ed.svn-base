function [x_Loc,Econstraint]=ArnoldActiveSet(mesh,parameters,A_Loc,b_Loc,x_Loc,E_contactLoc,E_contactGlob,E_InternalGlob,NELoc,edge_in_contact,Econstraint)

if(edge_in_contact==0)
x_Loc=A_Loc\b_Loc;
else
correction_Loc=1; 
contatore=0;
alpha=0.5;
while(norm(correction_Loc)>0.000000001 )
    % compute A(x+correction) = b, A correction = b - A correction
    Ax_Loc= A_Loc * x_Loc;
    rhs_Loc = b_Loc - Ax_Loc;


B_Loc=[];
C_Loc=[];
D_Loc=[];

    if(parameters.frictionless==1)
    % Constraint on the tangent component of the vector, 
    % (sigma n )_tangent =sigma n - (n'sigma n) n =0
    B_loc=sparse(edge_in_contact,length(A_Loc(:,1)) );
    for ee=1:edge_in_contact
        
        ee_loc=E_contactLoc(ee);
        dof=E_InternalGlob(ee_loc);
        normal_edge=mesh.normal_edge{dof};

        if(abs(normal_edge(1))<0.00001)
        B_Loc(ee, [ee_loc,ee_loc+NELoc]) = [(1-normal_edge(1)^2),         -normal_edge(1)*normal_edge(2)];        
        else
        B_Loc(ee, [ee_loc,ee_loc+NELoc]) = [-normal_edge(1)*normal_edge(2)          (1-normal_edge(2)^2)];            
        end
        
        
    end
    
    end
  
    EconstraintLoc=[];
    EfreeLoc=[];
    cont_constraint=0;
    cont_free=0;
    for ee=1:edge_in_contact
        ee_loc=E_contactLoc(ee);
        ee_glob=E_contactGlob(ee);
        % this edge is already constrained, so (n' sigma n ) = 0
        if(Econstraint(ee_glob)==1)
        EconstraintLoc=[Econstraint;ee_loc];
        normal_edge=mesh.normal_edge{ee_glob};
        cont_constraint=cont_constraint+1;
        C_Loc(cont_constraint, :)=sparse(1,length(A_Loc));
        C_Loc(cont_constraint, [ee_loc,ee_loc+NELoc])=[normal_edge(1),normal_edge(2)];
        else
        cont_free=cont_free+1;
        EfreeLoc = [ EfreeLoc; E_contactLoc(ee)];
        D_Loc(cont_free, :)=sparse(1,length(A_Loc));
        D_Loc(cont_free, [ee_loc,ee_loc+NELoc])=[normal_edge(1),normal_edge(2)];
        end
    end
    EfreeLoc=setdiff(E_contactLoc,EconstraintLoc);
    
                                                
                                                
    if(Nconstraint())
       N_Loc=sparse(1,length(A_Loc));
       normal_node=mesh.normal_node{nn_glob};

       N_Loc(1,[2 * NELoc + nn, 2 * NELoc + NLoc + nn])=[normal_node(1), normal_node(2)];
    end
                                                
                                                
    M_Loc=[A_Loc,  -B_Loc',                 -C_Loc';
           B_Loc,  sparse(edge_in_contact,  edge_in_contact+cont_constraint );
           C_Loc,  sparse(cont_constraint,  edge_in_contact+cont_constraint )];
       
    f_Loc=[rhs_Loc;sparse(edge_in_contact+cont_constraint,1)];
    
    lagrange_correction_Loc=M_Loc\f_Loc;
    correction_Loc=lagrange_correction_Loc(1:length(x_Loc));
    
    alpha=1;
    blocking_constraint=[];
    for ee=1:length(EfreeLoc)
        alpha_ee= -( D_Loc(ee,:) * x_Loc ) / ( D_Loc(ee,:)*correction_Loc );
        
        if(alpha_ee<alpha&&alpha_ee>=0)
            alpha=alpha_ee;
            blocking_constraint=ee;
        end
    end
    
    
    if(NconstraintFree==1)
        alpha_nn= (gap-N_Loc(1,:) * x_Loc ) / ( D_Loc(ee,:)*correction_Loc );
        if(alpha_nn<alpha&&alpha_nn>=0)
        alpha=alpha_nn;
        blocking_constraint=nn;
        end
    end
                                                
                                                
    % update the constraint that has been reached
    Econstraint(E_contactGlob(blocking_constraint))=1;
    x_Loc = x_Loc + alpha * correction_Loc;

    contatore=contatore+1;
    [contatore,alpha,norm(correction_Loc)]
end

end

end
